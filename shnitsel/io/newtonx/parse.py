from io import TextIOWrapper
import pathlib
from typing import NamedTuple, Tuple
import numpy as np
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults
from shnitsel.units.definitions import get_default_input_attributes
import xarray as xr
from itertools import combinations
import pandas as pd
import logging
import os
import re
import math

from shnitsel.io.helpers import LoadingParameters, PathOptionsType, make_uniform_path

from ..xyz import parse_xyz


# TODO: FIXME: use loading_parameters to configure units and state names
def parse_newtonx(
    traj_path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read a NewtonX trajectory directory into a Dataset with standard shnitsel annotations and units

    Args:
        pathlist (PathOptionsType): Path to the NewtonX trajectory output
        loading_parameters (LoadingParameters | None, optional): Parameter settings for e.g. standard units or state names.

    Returns:
        xr.Dataset: The Dataset object containing all of the loaded data in default shnitsel units
    """
    path_obj = make_uniform_path(traj_path)

    # TODO: FIXME: use loading_parameters to configure state names
    with open(path_obj / "RESULTS" / "nx.log") as f:
        settings = parse_settings_from_nx_log(f)
        trajectory = create_initial_dataset(
            settings.num_states, settings.num_atoms, loading_parameters
        )

        default_attributes = get_default_input_attributes("newtonx", loading_parameters)
        # Add time dimension

        isolated_keys = [
            "atNames",
            "atNums",
            "state",
            "statecomb",
            "state_names",
            "state_types",
        ]

        trajectory = trajectory.set_coords(isolated_keys)

        trajectory = trajectory.expand_dims({"time": settings.num_steps}, axis=0)
        trajectory.attrs["delta_t"] = settings.delta_t
        # TODO: FIXME: Make sure this is actually a common naming scheme
        trajectory.attrs["real_tmax"] = settings.t_max
        trajectory.attrs["max_ts"] = settings.num_steps
        trajectory.attrs["completed"] = settings.completed
        # Create "active state" variable
        trajectory = trajectory.assign(
            {
                "astate": xr.DataArray(
                    np.zeros((settings.num_steps)),
                    dims=["time"],
                    name="astate",
                    attrs=default_attributes["astate"],
                )
            }
        )

        # trajectory = trajectory.assign_coords(time=("time", [0.0]))

        # Read several datasets into trajectory and get first indicator of actual performed steps
        actual_steps = parse_nx_log_data(f, trajectory)

        if actual_steps < settings.num_steps:
            # Filter
            trajectory.attrs["completed"] = False
            trajectory = trajectory.isel(time=slice(0, actual_steps))
        elif actual_steps > settings.num_steps:
            raise ValueError(
                f"Trajectory data at {path_obj} contained data for {actual_steps} frames instead of initially denoted {settings.num_steps} frames. Cannot allocate space after the fact."
            )

        # Assign times to trajectory
        trajectory.assign_coords(time=np.arange(0, actual_steps) * settings.delta_t)

    # nsteps = single_traj.sizes['ts']

    with open(os.path.join(traj_path, "RESULTS", "dyn.xyz")) as f:
        atNames, atNums, atXYZ = parse_xyz(f)

    # if (
    #     not single_traj.attrs["completed"]
    #     and atXYZ.shape[0] == single_traj.sizes["time"] + 1
    # ):
    #     logging.info("Geometry file contains time after error. Truncating.")
    #     atXYZ = atXYZ[:-1]

    trajectory.atNames.values = atNames
    trajectory.atNums.values = atNums
    trajectory.atXYZ.values = atXYZ

    trajectory["time"].attrs.update(default_attributes[str("time")])
    return convert_all_units_to_shnitsel_defaults(trajectory)


class NewtonXSettingsResult(NamedTuple):
    t_max: float
    delta_t: float
    num_steps: int
    num_atoms: int
    num_states: int
    completed: bool


def parse_settings_from_nx_log(f) -> NewtonXSettingsResult:
    completed = True

    # hack to deal with restarts obscuring real tmax
    real_tmax = float(0)
    for line in f:
        stripline = line.strip()
        if stripline.startswith("FINISHING"):
            real_tmax = max(real_tmax, float(stripline.split()[4]))
        elif stripline.startswith("xx:"):
            splitline = stripline.split()
            if splitline[1] == "::ERROR::":
                real_tmax = max(real_tmax, float(splitline[7]))
                completed = False

    logging.debug(f"found real_tmax: {real_tmax}")
    f.seek(0)

    # skip to settings
    for line in f:
        if line.strip().startswith("Initial parameters:"):
            break

    settings = {}
    for line in f:
        # blank line marks end of settings
        if line.strip() == "":
            break

        key, val = re.split(" *= *", line.strip())
        if "." in val:
            val = float(val)
        else:
            val = int(val)
        settings[key] = val

    delta_t = settings["dt"]
    max_ts = int(real_tmax / delta_t)
    nsteps = max_ts + 1
    nstates = settings["nstat"]
    natoms = settings["Nat"]
    assert isinstance(nstates, int)
    assert isinstance(natoms, int)

    return NewtonXSettingsResult(real_tmax, delta_t, nsteps, natoms, nstates, completed)


def create_initial_dataset(
    nstates: int,
    natoms: int,
    loading_parameters: LoadingParameters | None,
    **kwargs,
) -> xr.Dataset:
    """Function to initialize an `xr.Dataset` with appropriate variables and coordinates to acommodate loaded mewtonx data.

    Args:
        nstates (int): Number of states within the datasets
        natoms (int): The number of atoms within the datasets

    Returns:
        xr.Dataset: An xarray Dataset with appropriately sized DataArrays and coordinates also including default attributes for all variables.
    """
    template = {
        "energy": ["state"],
        # "dip_all": ["state", "state2", "direction"],
        "dip_perm": ["state", "direction"],
        "dip_trans": ["statecomb", "direction"],
        # Has only active state forces
        "forces": ["atom", "direction"],  # ["state", "atom", "direction"],
        # 'has_forces': ['placeholder'],
        # 'has_forces': [],
        "phases": ["state"],
        "nacs": ["statecomb", "atom", "direction"],
        "atXYZ": ["atom", "direction"],
        "atNames": ["atom"],
        "atNums": ["atom"],
        "state_names": ["state"],
        "state_types": ["state"],
    }

    template_default_values = {
        "energy": np.nan,
        # "dip_all": np.nan,
        "dip_perm": np.nan,
        "dip_trans": np.nan,
        "forces": np.nan,
        "phases": np.nan,
        "nacs": np.nan,
        "atXYZ": np.nan,
        "atNames": "",
        "atNums": -1,
        "state_names": "",
        "state_types": 0,
    }

    if natoms == 0:
        # No atoms expected, so let us delete all variables with atom dimension
        del template["forces"]
        del template["nacs"]
        del template["atXYZ"]
        del template["atNames"]
        del template["atNums"]

    if nstates == 0:
        # On the other hand, we don't worry about not knowing nstates,
        # because energy is always written.
        pass

    dim_lengths = {
        "state": nstates,
        "state2": nstates,
        "atom": natoms,
        "direction": 3,
        "statecomb": math.comb(nstates, 2),
    }

    coords: dict | xr.Dataset = {
        "state": (states := np.arange(1, nstates + 1)),
        "state2": states,
        "atom": np.arange(natoms),
        "direction": ["x", "y", "z"],
    }

    coords = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(states, 2), names=["from", "to"]),
        dim="statecomb",
    ).merge(coords)

    # attrs = {
    #    'atXYZ': {'long_name': "positions", 'units': 'Bohr', 'unitdim': 'Length'},
    #    'energy': {'units': 'hartree', 'unitdim': 'Energy'},
    #    'e_kin': {'units': 'hartree', 'unitdim': 'Energy'},
    #    'dip_perm': {'long_name': "permanent dipoles", 'units': 'au'},
    #    'dip_trans': {'long_name': "transition dipoles", 'units': 'au'},
    #    'sdiag': {'long_name': 'active state (diag)'},
    #    'astate': {'long_name': 'active state (MCH)'},
    #    'forces': {'units': 'hartree/bohr', 'unitdim': 'Force'},
    #    'nacs': {'long_name': "nonadiabatic couplings", 'units': 'au'},
    # }
    default_attributes = get_default_input_attributes("newtonx", loading_parameters)

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
            (default_attributes[varname] if varname in default_attributes else {}),
        )
        for varname, dims in template.items()
    }

    res_dataset = xr.Dataset(datavars, coords)

    # Try and set some default attributes on the coordinates for the dataset
    for coord_name in res_dataset.coords:
        if coord_name in default_attributes:
            res_dataset[coord_name].attrs.update(default_attributes[str(coord_name)])

    res_dataset.attrs["input_format"] = "newtonx"
    res_dataset.attrs["input_type"] = "dynamic"
    res_dataset.attrs["has_forces"] = "active_only"

    return res_dataset


def parse_nx_log_data(f: TextIOWrapper, dataset: xr.Dataset) -> int:
    # f should be after the setting

    natoms = dataset.sizes["atom"]
    nstates = dataset.sizes["state"]
    nstatecomb = dataset.sizes["statecomb"]
    ntimesteps = dataset.sizes["time"]

    actual_max_ts: int = -1

    delta_t = dataset.attrs["delta_t"]

    ts: int = 0
    time: float = 0

    tmp_astate = np.zeros((ntimesteps,))
    tmp_forces = np.zeros((natoms, 3))
    tmp_energy = np.zeros((nstates,))
    tmp_nacs = np.zeros((nstatecomb, natoms, 3))
    tmp_full_forces = np.zeros((ntimesteps, natoms, 3))
    tmp_full_energy = np.zeros(
        (
            ntimesteps,
            nstates,
        )
    )
    tmp_full_nacs = np.zeros((ntimesteps, nstatecomb, natoms, 3))

    # See page 101, section 16.3 of newtonX documentation for order of output values.
    # parse actual data
    for line in f:
        stripline = line.strip()
        # TODO: FIXME: unfortunately, Newton-X reports current step number, time and
        #       active state at the end of a timestep.
        #       Currently we use this to set the step number and time for the
        #       following time step. This is confusing and possibly vulnerable to
        #       strangely-formatted data -- is there a guarantee that time steps are
        #       always in order? There's almost certainly a better way to do this.
        #       For example, instead of assuming times follow expected order,
        #       lookahead to next FINISHING
        # THIS MIGHT BE SOLVED WITH ABOVE TMP STORAGE AND ONLY WRITING UPON READING THE FINISHING LINE
        if stripline.startswith("FINISHING"):
            # Set active state for _current_ time step
            tmp_astate[ts] = int(stripline.split()[8])

            # Set step number and time for _the following_ time step
            time = float(stripline.split()[4])
            # TODO: Why can't we take the time step denoted in the log?
            ts = int(round(time / delta_t))
            # logging.debug(f"{ts}=round({time}/{delta_t})")

            actual_max_ts = max(actual_max_ts, ts)
            tmp_full_forces[ts] = tmp_forces
            tmp_full_energy[ts] = tmp_energy
            tmp_full_nacs[ts] = tmp_nacs
            logging.debug(f"finished ts {ts - 1}")

        elif stripline.startswith("Gradient vectors"):
            for iatom in range(natoms):
                tmp_forces[iatom] = [float(n) for n in next(f).strip().split()]

        elif stripline.startswith("Nonadiabatic coupling vectors"):
            # TODO: FIXME: Are we sure that all NACS are in identical order in all formats?
            for icomb in range(math.comb(nstates, 2)):
                # Order is: V(from, to),iatom, dir
                # Increase steps from rightmost dimension to leftmost.
                # I.e. each line has x,y,z value.
                # First natoms lines have the values for the different atoms in first state combination
                # Each block of natoms represents successive state combinations
                for iatom in range(natoms):
                    tmp_nacs[icomb, iatom] = [float(n) for n in next(f).strip().split()]

        elif stripline.startswith("Energy ="):
            for istate in range(nstates):
                tmp_energy[istate] = float(next(f).strip())

    dataset.assign_coords(astate=tmp_astate)
    dataset.forces.values = tmp_full_forces
    dataset.energy.values = tmp_full_energy
    dataset.nacs.values = tmp_full_nacs

    return actual_max_ts + 1

    # TODO: Are these units even correct?
    # return #xr.Dataset(
    #     #{
    #         "energy": (
    #             ["ts", "state"],
    #             energy,
    #             {"units": "hartree", "unitdim": "Energy"},
    #         ),
    #         # 'dip_all': (['ts', 'state', 'state2', 'direction'], dip_all),
    #         # 'dip_perm': (['ts', 'state', 'direction'], dip_perm),
    #         # 'dip_trans': (['ts', 'statecomb', 'direction'], dip_trans),
    #         # 'sdiag': (['ts'], sdiag),
    #         "astate": (["ts"], astate, {"long_name": "active state"}),
    #         "forces": (
    #             ["ts", "atom", "direction"],
    #             forces,
    #             {"units": "hartree/bohr", "unitdim": "Force"},
    #         ),
    #         # 'has_forces': (['ts'], has_forces),
    #         # 'phases': (['ts', 'state'], phases),
    #         "nacs": (
    #             ["ts", "statecomb", "atom", "direction"],
    #             nacs,
    #             {"long_name": "nonadiabatic couplings", "units": "au"},
    #         ),
    #     },
    #     coords=coords,
    #     attrs={
    #         "max_ts": max_ts,
    #         "real_tmax": real_tmax,
    #         "delta_t": delta_t,
    #         "completed": completed,
    #     },
    # )
