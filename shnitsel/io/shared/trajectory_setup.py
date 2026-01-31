from dataclasses import asdict, dataclass
from itertools import combinations, permutations
import math
from typing import Dict, List, Literal, Tuple, TypeVar

import pandas as pd
import xarray as xr

import numpy as np

from shnitsel.data.atom_helpers import (
    get_atom_number_from_symbol,
    get_symbol_from_atom_number,
)
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.io.shared.helpers import LoadingParameters
from shnitsel.io.shared.variable_flagging import (
    is_variable_assigned,
    mark_dataset_for_cleanup,
    mark_variable_assigned,
)

TrajType = TypeVar("TrajType", bound=xr.Dataset)


@dataclass
class RequiredTrajectorySettings:
    """Helper class to make sure all required settings are applied and to
    avoid typos when assigning required settings to consistent keys in trajectories.

    """

    t_max: float
    delta_t: float
    max_ts: int
    completed: bool
    input_format: Literal["sharc", "newtonx", "ase", "pyrai2md"]
    input_type: Literal["static", "dynamic", "unknown"]
    input_format_version: str

    num_singlets: int
    num_doublets: int
    num_triplets: int


@dataclass
class OptionalTrajectorySettings:
    """Helper class to make sure consistent keys are used in trajectory settings
    for optional configuration information.
    """

    has_forces: bool | Literal["all", "active_only"] | None = None
    trajectory_id: int | None = None
    is_multi_trajectory: bool | None = None
    trajectory_input_path: str | None = None

    theory_basis_set: str | None = None
    est_level: str | None = None

    # To keep track of input settings we do not explicitly use anywhere else.
    misc_input_settings: Dict | None = None


def assign_required_settings(
    dataset: TrajType, settings: RequiredTrajectorySettings
) -> TrajType:
    """
    Function to assign all required settings to the dataset.
    Some of the settings may be applied to a coordinate variable of zero-dimensional shape instead.

    Just a handy tool so all values are assigned because all fields in `settings` should be assigned upon its creation

    Parameters
    ----------
    dataset : xr.Dataset | ShnitselDataset
        The dataset/trajectory to write the required settings into
    settings : RequiredTrajectorySettings
        The fully assigned settings object containing all keys and values to be assigned.
    """
    kv_dict = dict(asdict(settings))

    if 'time' in dataset.coords:
        time_unit = dataset.time.attrs.get('units', 'fs')
    else:
        time_unit = 'fs'

    drop_keys = []
    new_coords = {}
    for var, var_type, attrs in [
        (
            'delta_t',
            float,
            {
                'long_name': 'The simulation time step between frames',
                'unitdim': 'time',
                'units': time_unit,
            },
        ),
        (
            'max_ts',
            int,
            {
                'long_name': 'The maximum time step found in the simulation data during import.'
            },
        ),
        (
            't_max',
            float,
            {
                'long_name': 'The maximum configured time for a simulation to run',
                'unitdim': 'time',
                'units': time_unit,
            },
        ),
    ]:
        if var in dataset.coords:
            if var_type(dataset[var]) == var_type(kv_dict[var]):
                # Correct value already set, no need to move things around
                continue
            else:
                drop_keys.append(var)

        new_coords[var] = xr.DataArray(
            var_type(kv_dict[var]), dims=(), name=var, attrs=attrs
        ).astype(var_type)
        mark_variable_assigned(new_coords[var])

    # Do not also set the metadata as attrs if we set it as variable.
    # Keeping it all consistent would be a mess.
    del kv_dict['t_max']
    del kv_dict['delta_t']
    del kv_dict['max_ts']

    res = dataset.drop_vars(drop_keys).assign_attrs(kv_dict).assign_coords(new_coords)

    # for key in new_coords.keys():
    #     mark_variable_assigned(res.coords[key])

    return res if isinstance(dataset, xr.Dataset) else wrap_dataset(res)


def assign_optional_settings(
    dataset: TrajType, settings: OptionalTrajectorySettings
) -> TrajType:
    """
    Function to assign all assigned optional settings to a dataset.

    Just a handy tool so we can be sure the settings are assigned with the correct keys.

    Parameters
    ----------
    dataset : xr.Dataset | Trajectory | Frames
        The dataset/trajectory to write the optional settings into
    settings : OptionalTrajectorySettings
        The dataclass object that has all optional setting keys with optional values. Only assigned settings (not None) will be inserted.
    """
    kv_dict = {k: v for k, v in asdict(settings).items() if v is not None}
    res = dataset.assign_attrs(kv_dict)
    return res if isinstance(dataset, xr.Dataset) else wrap_dataset(res)


def get_statecomb_coordinate(states: xr.DataArray) -> xr.Coordinates:
    """Helper function to create a statecombination coordinate if it is missing, based on the states registered.

    Parameters
    ----------
    states : xr.DataArray
        The state coordinate

    Returns
    -------
    xr.Coordinates
        The new coordinate having all non-ordered state combinations
    """
    return xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(
            combinations(states.values, 2), names=["statecomb_from", "statecomb_to"]
        ),
        dim="statecomb",
    )


# TODO: FIXME: Consider caching here. May speed up ICONDS in large ICOND datasets
def create_initial_dataset(
    num_time_steps: int,
    num_states: int,
    num_atoms: int,
    format_name: Literal["sharc", "newtonx", "ase", "pyrai2md"],
    loading_parameters: LoadingParameters | None,
    **kwargs,
) -> Tuple[xr.Dataset, dict]:
    """
    Function to initialize an `xr.Dataset` with appropriate variables and coordinates to acommodate loaded data.

    All arguments are used to accurately size the dimensions of the dataset or assign.
    Also returns the default attributes associated with the variables in the dataset for later assignment of certain values like "time" which is not initialized with values.


    Parameters
    ----------
    num_time_steps : int
        The number of expected time steps in this trajectory. Set to 0 to not create a "time" dimension.
    num_states : int
        The number of states within the datasets.
    num_atoms : int
        _description_
    format_name : Literal["sharc", "newtonx", "ase", "pyrai2md"]
        The number of atoms within the datasets. Set to 0 to remove all observables tied to an "atom" index.
    loading_parameters : LoadingParameters | None
        Optional parameters to influence import and object creation. Used to set some input units and default handlers.

    Returns
    -------
    Tuple[xr.Dataset, dict]
        First An xarray Dataset with appropriately sized DataArrays and coordinates also including default attributes for all variables.
        Second the key-value dict, where the key is the name of standard variables/coordinates in Shnitsel terminology and the value is the dict of default attributes associated with this variable in this format.
    """
    from shnitsel.units.defaults import get_default_input_attributes

    # This is the list of observables/variables we currently support.
    template = {
        "energy": ["time", "state"],
        "e_kin": ["time"],
        "velocities": ["time", "atom", "direction"],
        "forces": ["time", "state", "atom", "direction"],
        "atXYZ": ["time", "atom", "direction"],
        "nacs": ["time", "statecomb", "atom", "direction"],
        # "dip_all": ["time", "state", "state2", "direction"],
        "dip_perm": ["time", "state", "direction"],
        "dip_trans": ["time", "statecomb", "direction"],
        # SOCs need to be a state x state (without repetition) matrix as it is not necessarily hermitean
        # "socs": ["time", "statecomb"],
        "socs": ["time", "full_statecomb"],
        "state_names": ["state"],
        "state_types": ["state"],
        "state_magnetic_number": ["state"],
        "state_degeneracy_group": ["state"],
        "state_charges": ["state"],
        "astate": ["time"],
        "sdiag": ["time"],
        "phases": ["time", "state"],
        "atNames": ["atom"],
        "atNums": ["atom"],
    }

    template_default_values = {
        "energy": np.nan,
        "e_kin": np.nan,
        "velocities": np.nan,
        "forces": np.nan,
        "atXYZ": np.nan,
        "nacs": np.nan,
        # "dip_all": np.nan,
        "dip_perm": np.nan,
        "dip_trans": np.nan,
        "socs": np.nan + 0j,
        "state_names": "",
        "state_types": 0,
        "state_magnetic_number": np.nan,
        "state_degeneracy_group": 0,
        "state_charges": 0,
        "astate": -1,
        "sdiag": -1,
        "phases": np.nan,
        "atNames": "",
        "atNums": -1,
    }

    default_float_type = np.dtypes.Float32DType
    default_string_type = "U8"

    template_default_dtypes = {
        "energy": default_float_type,
        "e_kin": default_float_type,
        "velocities": default_float_type,
        "forces": default_float_type,
        "atXYZ": default_float_type,
        "nacs": default_float_type,
        # "dip_all": np.float32,
        "dip_perm": default_float_type,
        "dip_trans": default_float_type,
        "socs": np.complex128,
        "state_names": default_string_type,
        "state_types": np.dtypes.Int8DType,
        "state_magnetic_number": default_float_type,
        "state_degeneracy_group": np.dtypes.Int8DType,
        "state_charges": default_float_type,
        "astate": np.dtypes.Int16DType,
        "sdiag": np.dtypes.Int16DType,
        "phases": default_float_type,
        "atNames": default_string_type,
        "atNums": np.dtypes.Int8DType,
    }

    dim_lengths = {
        "time": num_time_steps,
        "state": num_states,
        "state2": num_states,
        "atom": num_atoms,
        "direction": 3,
        "statecomb": math.comb(num_states, 2),
        "full_statecomb": math.perm(num_states, 2),
    }

    coords: dict | xr.Dataset = {
        "state": (states := np.arange(1, num_states + 1)),
        "state2": states,
        "atom": np.arange(num_atoms),
        "direction": ["x", "y", "z"],
    }

    def template_purge_dim(template_dict: Dict[str, List[str]], dim: str):
        """Helper function to remove all variables dependent only on a certain dimension or remove the dimension from other variables' index list.

        Args:
            template_dict (Dict[str, List[str]]): The current state of the template dictionary
            dim (str): The dimension key to purge.
        """
        obsolete_keys = []
        for k, v in template_dict.items():
            if dim in v:
                if len(v) == 1:
                    obsolete_keys.append(k)
                else:
                    template_dict[k].remove(dim)

        for key in obsolete_keys:
            del template_dict[key]

    if num_time_steps == 0:
        template_purge_dim(template, "time")
        del dim_lengths["time"]
        # del coords["time"]

    if num_states == 0:
        # On the other hand, we don't worry about not knowing nstates,
        # because energy is always written.
        pass

    if num_atoms == 0:
        template_purge_dim(template, "atom")
        del dim_lengths["atom"]
        del coords["atom"]

    coords = (
        xr.Coordinates.from_pandas_multiindex(
            pd.MultiIndex.from_tuples(combinations(states, 2), names=["from", "to"]),
            dim="statecomb",
        )
        .merge(
            xr.Coordinates.from_pandas_multiindex(
                pd.MultiIndex.from_tuples(
                    permutations(states, 2),
                    names=["full_statecomb_from", "full_statecomb_to"],
                ),
                dim="full_statecomb",
            )
        )
        .merge(coords)
    )

    default_format_attributes = get_default_input_attributes(
        format_name, loading_parameters
    )
    # logging.debug(default_format_attributes)

    datavars = {
        varname: (
            dims,
            (
                x
                if (x := kwargs.get(varname)) is not None
                else np.full(
                    [dim_lengths[d] for d in dims],
                    fill_value=template_default_values[varname],
                    dtype=template_default_dtypes[varname],
                )
            ),
            (
                default_format_attributes[varname]
                if varname in default_format_attributes
                else {}
            ),
        )
        for varname, dims in template.items()
    }

    res_dataset = xr.Dataset(datavars, coords)

    # Try and set some default attributes on the coordinates for the dataset
    for coord_name in res_dataset.coords:
        if coord_name in default_format_attributes:
            res_dataset[coord_name].attrs.update(
                default_format_attributes[str(coord_name)]
            )

    res_dataset.attrs["input_format"] = format_name

    # Assign some of the variables as coordinates
    isolated_keys = list(
        set(
            [
                "atNames",
                "atNums",
                "statecomb",
                "state_names",
                "state_charges",
                "state_magnetic_number",
                "state_degeneracy_group",
                "state_types",
            ]
        ).intersection(template.keys())
    )

    # Prevent dimension labels from being lost.
    if "state" in res_dataset:
        mark_variable_assigned(res_dataset.state)
    if "atom" in res_dataset:
        mark_variable_assigned(res_dataset.atom)
    mark_variable_assigned(res_dataset.direction)
    if "statecomb" in res_dataset:
        mark_variable_assigned(res_dataset.statecomb)
        mark_variable_assigned(res_dataset["from"])
        mark_variable_assigned(res_dataset["to"])
    if "state_charges" in res_dataset:
        mark_variable_assigned(res_dataset.state_charges)

    mark_dataset_for_cleanup(res_dataset)

    res_dataset = res_dataset.set_coords(isolated_keys)
    res_dataset.attrs["_shnitsel_setup_for_cleanup"] = True

    return res_dataset, default_format_attributes


def fill_missing_dataset_variables(ds: xr.Dataset) -> xr.Dataset:
    """Helper function to fill in missing variables in datasets that we need for initialization.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to complete with state variables and atom variables

    Returns
    -------
    xr.Dataset
        The completed dataset with default variables applied.
    """
    from shnitsel.units.defaults import get_default_input_attributes

    template = {
        "state_names": ["state"],
        "state_types": ["state"],
        "state_magnetic_number": ["state"],
        "state_degeneracy_group": ["state"],
        "state_charges": ["state"],
        "atNames": ["atom"],
        "atNums": ["atom"],
    }

    template_default_values = {
        "state_names": "",
        "state_types": 0,
        "state_magnetic_number": np.nan,
        "state_degeneracy_group": 0,
        "state_charges": 0,
        "atNames": "",
        "atNums": -1,
    }

    default_float_type = np.dtypes.Float32DType
    default_string_type = "U8"

    template_default_dtypes = {
        "state_names": default_string_type,
        "state_types": np.dtypes.Int8DType,
        "state_magnetic_number": default_float_type,
        "state_degeneracy_group": np.dtypes.Int8DType,
        "state_charges": default_float_type,
        "atNames": default_string_type,
        "atNums": np.dtypes.Int8DType,
    }

    num_states = ds.sizes["state"] if "state" in ds.sizes else 0
    num_atoms = ds.sizes["atom"] if "atom" in ds.sizes else 0

    dim_lengths = {
        "state": num_states,
        "atom": num_atoms,
    }

    isolated_keys = list(
        set(
            [
                "atNames",
                "atNums",
                "state_names",
                "state_charges",
                "state_magnetic_number",
                "state_degeneracy_group",
                "state_types",
            ]
        ).intersection(template.keys())
    )

    default_format_attributes = get_default_input_attributes("shnitsel", None)

    for var in ds.data_vars:
        mark_variable_assigned(ds[var])
    for var in ds.coords:
        mark_variable_assigned(ds[var])

    var_data = {
        varname: (
            dims,
            (
                np.full(
                    [dim_lengths[d] for d in dims],
                    fill_value=template_default_values[varname],
                    dtype=template_default_dtypes[varname],
                )
            ),
            (
                default_format_attributes[varname]
                if varname in default_format_attributes
                else {}
            ),
        )
        for varname, dims in template.items()
        if varname not in ds
    }
    if len(var_data) == 0:
        return ds

    res_dataset = ds.assign(var_data)

    if "atNums" in var_data and "atNames" not in var_data:
        nums = [get_atom_number_from_symbol(x) for x in ds.atNames.values]
        res_dataset = res_dataset.assign(
            atNums=(res_dataset.atNums.dims, nums, res_dataset.atNums.attrs)
        )
        mark_variable_assigned(res_dataset["atNums"])

    if "atNames" in var_data and "atNums" not in var_data:
        names = [get_symbol_from_atom_number(x) for x in ds.atNums.values]
        res_dataset = res_dataset.assign(
            atNames=(res_dataset.atNames.dims, names, res_dataset.atNames.attrs)
        )
        mark_variable_assigned(res_dataset["atNames"])

    res_dataset = res_dataset.set_coords(isolated_keys)

    if var_data:
        mark_dataset_for_cleanup(res_dataset)

    return res_dataset
