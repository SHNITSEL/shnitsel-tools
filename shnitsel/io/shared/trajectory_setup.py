from dataclasses import asdict, dataclass
from itertools import combinations
import logging
import math
from typing import Dict, List, Literal, Tuple

import pandas as pd
import xarray as xr

import numpy as np

from shnitsel.io.helpers import LoadingParameters
from shnitsel.io.shared.variable_flagging import mark_variable_assigned


@dataclass
class RequiredTrajectorySettings:
    t_max: float
    delta_t: float
    max_ts: int
    completed: bool
    input_format: Literal["sharc", "newtonx", "ase", "pyrai2md"]
    input_type: Literal['static', 'dynamic']
    input_format_version: str

    num_singlets: int
    num_doublets: int
    num_triplets: int


@dataclass
class OptionalTrajectorySettings:
    has_forces: bool | Literal['all', 'active_only'] | None = None
    trajid: int | None = None
    is_multi_trajectory: bool | None = None
    trajectory_input_path: str | None = None

    theory_basis_set: str | None = None
    est_level: str | None = None

    # To keep track of input settings we do not explicitly use anywhere else.
    misc_input_settings: Dict | None = None


def assign_required_settings(
    dataset: xr.Dataset | xr.DataTree, settings: RequiredTrajectorySettings
) -> None:
    """Function to assign all required settings to the dataset.

    Just a handy tool so all values are assigned because all fields in `settings` should be assigned upon its creation

    Args:
        dataset (xr.Dataset): The dataset to write the required settings into
        settings (RequiredTrajectorySettings): The fully assigned settings object containing all keys and values to be assigned.
    """
    dataset.attrs.update(asdict(settings))


def assign_optional_settings(
    dataset: xr.Dataset | xr.DataTree, settings: OptionalTrajectorySettings
) -> None:
    """Function to assign all assigned optional settings to a dataset.

    Just a handy tool so we can be sure the settings are assigned with the correct keys.

    Args:
        dataset (xr.Dataset): The dataset to write the optional settings into
        settings (OptionalTrajectorySettings): The dataclass object that has all optional setting keys with optional values. Only assigned settings (not None) will be inserted.
    """
    kv_dict = asdict(settings)
    for k, v in kv_dict.items():
        if v is not None:
            dataset.attrs[k] = v


# TODO: FIXME: Consider caching here. May speed up ICONDS in large ICOND datasets
def create_initial_dataset(
    num_time_steps: int,
    num_states: int,
    num_atoms: int,
    format_name: Literal["sharc", "newtonx", "ase", "pyrai2md"],
    loading_parameters: LoadingParameters | None,
    **kwargs,
) -> Tuple[xr.Dataset, Dict]:
    """Function to initialize an `xr.Dataset` with appropriate variables and coordinates to acommodate loaded data.

    All arguments are used to accurately size the dimensions of the dataset or assign.
    Also returns the default attributes associated with the variables in the dataset for later assignment of certain values like "time" which is not initialized with values.

    Args:
        num_time_steps (int): The number of expected time steps in this trajectory. Set to 0 to not create a "time" dimension.
        num_states (int): The number of states within the datasets.
        num_atoms (int): The number of atoms within the datasets. Set to 0 to remove all observables tied to an "atom" index.

    Returns:
        xr.Dataset: An xarray Dataset with appropriately sized DataArrays and coordinates also including default attributes for all variables.
        Dict: The key-value dict, where the key is the name of standard variables/coordinates in Shnitsel terminology and the value is the dict of default attributes associated with this variable in this format.
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
        # TODO: FIXME: Check the correct dimensions for socs
        "socs": ["time", "statecomb", "atom", "direction"],
        "state_names": ["state"],
        "state_types": ["state"],
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
        "socs": np.nan,
        "state_names": "",
        "state_types": 0,
        "state_charges": 0,
        "astate": -1,
        "sdiag": -1,
        "phases": np.nan,
        "atNames": "",
        "atNums": -1,
    }

    dim_lengths = {
        "time": num_time_steps,
        "state": num_states,
        "state2": num_states,
        "atom": num_atoms,
        "direction": 3,
        "statecomb": math.comb(num_states, 2),
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

    coords = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(states, 2), names=["from", "to"]),
        dim="statecomb",
    ).merge(coords)

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
    isolated_keys = [
        "atNames",
        "atNums",
        "statecomb",
        "state_names",
        "state_charges",
        "state_types",
    ]

    # Prevent dimension labels from being lost.
    mark_variable_assigned(res_dataset.state)
    mark_variable_assigned(res_dataset.atom)
    mark_variable_assigned(res_dataset.direction)
    mark_variable_assigned(res_dataset.statecomb)
    mark_variable_assigned(res_dataset["from"])
    mark_variable_assigned(res_dataset["to"])
    mark_variable_assigned(res_dataset.state_charges)

    res_dataset = res_dataset.set_coords(isolated_keys)

    return res_dataset, default_format_attributes
