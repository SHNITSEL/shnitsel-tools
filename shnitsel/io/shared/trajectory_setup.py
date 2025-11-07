
from dataclasses import asdict, dataclass
from itertools import combinations
import math
from typing import Dict, List, Literal

import pandas as pd
import xarray as xr

import numpy as np

from shnitsel.io.helpers import LoadingParameters
from shnitsel.units.definitions import get_default_input_attributes


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
    has_forces: bool | None = None
    trajid: int | None = None
    is_multi_trajectory: bool | None = None
    trajectory_input_path: str | None = None


def assign_required_settings(dataset: xr.Dataset, settings: RequiredTrajectorySettings) -> None:
    """Function to assign all required settings to the dataset.

    Just a handy tool so all values are assigned because all fields in `settings` should be assigned upon its creation

    Args:
        dataset (xr.Dataset): The dataset to write the required settings into
        settings (RequiredTrajectorySettings): The fully assigned settings object containing all keys and values to be assigned.
    """
    dataset.attrs.update(asdict(settings))


def assign_optional_settings(dataset: xr.Dataset, settings: OptionalTrajectorySettings) -> None:
    """Function to assign all assigned optional settings to a dataset.

    Just a handy tool so we can be sure the settings are assigned with the correct keys.

    Args:
        dataset (xr.Dataset): The dataset to write the optional settings into
        settings (OptionalTrajectorySettings): The dataclass object that has all optional setting keys with optional values. Only assigned settings (not None) will be inserted.
    """
    for k, v in asdict(settings):
        if v is not None:
            dataset.attrs[k] = v


def create_initial_dataset(
    num_time_steps: int,
    num_states: int,
    num_atoms: int,
    format_name: Literal["sharc", "newtonx", "ase", "pyrai2md"],
    loading_parameters: LoadingParameters | None,
    **kwargs,
) -> xr.Dataset:
    """Function to initialize an `xr.Dataset` with appropriate variables and coordinates to acommodate loaded data.

    All arguments are used to accurately size the dimensions of the dataset or assign 

    Args:
        num_time_steps (int): The number of expected time steps in this trajectory. Set to 0 to not create a "time" dimension.
        num_states (int): The number of states within the datasets. 
        num_atoms (int): The number of atoms within the datasets. Set to 0 to remove all observables tied to an "atom" index.

    Returns:
        xr.Dataset: An xarray Dataset with appropriately sized DataArrays and coordinates also including default attributes for all variables.
    """
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
        "astate": -1,
        "sstate": -1,
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
                    obsolete_keys.append(v)
                else:
                    template_dict[k].remove(dim)

        for key in obsolete_keys:
            del template_dict[key]

    if num_time_steps == 0:
        template_purge_dim(template, "time")
        del dim_lengths["time"]
        del coords["time"]

    if num_states == 0:
        # On the other hand, we don't worry about not knowing nstates,
        # because energy is always written.
        pass

    if num_atoms == 0:
        template_purge_dim(template, "atom")
        del dim_lengths["atom"]
        del coords["atom"]

    coords = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(
            states, 2), names=["from", "to"]),
        dim="statecomb",
    ).merge(coords)

    default_format_attributes = get_default_input_attributes(
        format_name, loading_parameters)

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
        "state_types",
    ]

    res_dataset = res_dataset.set_coords(isolated_keys)

    return res_dataset
