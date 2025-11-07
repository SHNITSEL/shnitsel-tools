from io import TextIOWrapper
import pathlib
from typing import NamedTuple, Tuple
import numpy as np
from shnitsel.io.shared.trajectory_setup import OptionalTrajectorySettings, RequiredTrajectorySettings, assign_optional_settings, assign_required_settings, create_initial_dataset
from shnitsel.io.shared.variable_flagging import is_variable_assigned, mark_variable_assigned
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

    # TODO: FIXME: Use other available files to read input if nx.log is not available, e.g. RESULTS/dyn.out, RESULTS/dyn.xyz, RESULTS/en.dat From those we can get most information anyway.
    # TODO: FIXME: Read basis from JOB_NAD files
    with open(path_obj / "RESULTS" / "nx.log") as f:
        settings = parse_settings_from_nx_log(f)
        default_attributes = get_default_input_attributes(
            "newtonx", loading_parameters)
        # Add time dimension

        trajectory = create_initial_dataset(
            settings.num_steps, settings.num_states, settings.num_atoms, "newtonx", loading_parameters
        )

        trajectory = trajectory.assign_coords(
            {
                "time": ("time",
                         np.arange(0, settings.num_steps) * settings.delta_t,
                         default_attributes[str("time")]
                         ),
            }
        )
        mark_variable_assigned(trajectory["time"])

        # Read several datasets into trajectory and get first indicator of actual performed steps
        actual_steps, trajectory = parse_nx_log_data(f, trajectory)

        if actual_steps < settings.num_steps:
            # Filter only assigned timesteps
            trajectory.attrs["completed"] = False
            trajectory = trajectory.isel(time=slice(0, actual_steps))
        elif actual_steps > settings.num_steps:
            raise ValueError(
                f"Trajectory data at {path_obj} contained data for {actual_steps} frames, which is more than the initially denoted {settings.num_steps} frames. Cannot allocate space after the fact."
            )

    with open(os.path.join(traj_path, "RESULTS", "dyn.xyz")) as f:
        atNames, atNums, atXYZ = parse_xyz(f)

    trajectory.atNames.values = atNames
    mark_variable_assigned(trajectory["atNames"])
    trajectory.atNums.values = atNums
    mark_variable_assigned(trajectory["atNums"])
    trajectory.atXYZ.values = atXYZ
    mark_variable_assigned(trajectory["atXYZ"])

    # Set all settings we require to be present on the trajectory
    # TODO: FIXME: Check if we can actually derive the number of singlets, doublets, or triplets from newtonx output.
    required_settings = RequiredTrajectorySettings(
        settings.t_max,
        settings.delta_t,
        min(settings.num_steps, actual_steps),
        settings.completed,
        "newtonx",
        "dynamic",
        settings.newtonx_version,
        -1,
        -1,
        -1)
    assign_required_settings(trajectory, required_settings)

    optional_settings = OptionalTrajectorySettings(
        has_forces=is_variable_assigned(trajectory["forces"]))
    assign_optional_settings(trajectory, optional_settings)

    return trajectory


class NewtonXSettingsResult(NamedTuple):
    t_max: float
    delta_t: float
    num_steps: int
    num_atoms: int
    num_states: int
    completed: bool
    newtonx_version: str


def parse_settings_from_nx_log(f) -> NewtonXSettingsResult:
    completed = True
    newtonx_version = "unknown"

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
        line_stripped = line.strip()
        # Newton-X version 2.2 (build 5, 2018-04-11)
        if line_stripped.startswith("Newton-X version"):
            parts = [x.strip() for x in line_stripped.split()]
            newtonx_version = parts[2]
        if line_stripped.startswith("version"):
            if newtonx_version == "unknown":
                parts = [x.strip() for x in line_stripped.split()]
                newtonx_version = parts[1].strip(",")
        if line_stripped.startswith("Initial parameters:"):
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

    return NewtonXSettingsResult(real_tmax, delta_t, nsteps, natoms, nstates, completed, newtonx_version)


def _create_initial_dataset(
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
        # NOTE: No dipoles in NewtonX output per default.
        # "dip_all": ["state", "state2", "direction"],
        # "dip_perm": ["state", "direction"],
        # "dip_trans": ["statecomb", "direction"],
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
        # "dip_perm": np.nan,
        # "dip_trans": np.nan,
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
        # "state2": nstates,
        "atom": natoms,
        "direction": 3,
        "statecomb": math.comb(nstates, 2),
    }

    coords: dict | xr.Dataset = {
        "state": (states := np.arange(1, nstates + 1)),
        # "state2": states,
        "atom": np.arange(natoms),
        "direction": ["x", "y", "z"],
    }

    coords = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(
            states, 2), names=["from", "to"]),
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
    default_attributes = get_default_input_attributes(
        "newtonx", loading_parameters)

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
            (default_attributes[varname]
             if varname in default_attributes else {}),
        )
        for varname, dims in template.items()
    }

    res_dataset = xr.Dataset(datavars, coords)

    # Try and set some default attributes on the coordinates for the dataset
    for coord_name in res_dataset.coords:
        if coord_name in default_attributes:
            res_dataset[coord_name].attrs.update(
                default_attributes[str(coord_name)])

    res_dataset.attrs["input_format"] = "newtonx"
    res_dataset.attrs["input_type"] = "dynamic"
    res_dataset.attrs["has_forces"] = "active_only"

    # TODO: FIXME: No dipoles read?
    logging.error(
        "Attention: NewtonX parsing currently does not read dipole information"
    )

    return res_dataset


def parse_nx_log_data(
    f: TextIOWrapper, dataset: xr.Dataset
) -> Tuple[int, xr.Dataset]:  # Tuple[int, xr.Dataset]:
    """Function to parse the nx.log data into a dataset from the input stream f.

    Will return the total number of actual timesteps read and the resulting dataset.
    Usual read data includes: forces, active state ("astate")

    Args:
        f (TextIOWrapper): Input filestream of a nx.log file
        dataset (xr.Dataset): The dataset to parse the data into

    Returns:
        Tuple[int, xr.Dataset]: The total number of actual timesteps read and the resulting dataset after applying modifications.
    """
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
    tmp_times = np.zeros((ntimesteps,))
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
            t_astate = int(stripline.split()[8])

            # Set step number and time for _the following_ time step
            t_time = float(stripline.split()[4])

            # Figure out current time step
            ts = int(round(t_time / delta_t))
            # Assign all values in this time step
            tmp_astate[ts] = t_astate
            tmp_times[ts] = tmp_times
            tmp_full_forces[ts] = tmp_forces
            tmp_full_energy[ts] = tmp_energy
            tmp_full_nacs[ts] = tmp_nacs

            actual_max_ts = max(actual_max_ts, ts)
            logging.debug(f"finished ts {ts}")

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
                    tmp_nacs[icomb, iatom] = [
                        float(n) for n in next(f).strip().split()]

        elif stripline.startswith("Energy ="):
            for istate in range(nstates):
                tmp_energy[istate] = float(next(f).strip())

    dataset = dataset.assign_coords(
        {
            "astate": ("time",
                       tmp_astate,
                       dataset["astate"].attrs
                       ),
        }
    )
    mark_variable_assigned(dataset["astate"])

    dataset.forces.values = tmp_full_forces
    mark_variable_assigned(dataset["forces"])
    dataset.energy.values = tmp_full_energy
    mark_variable_assigned(dataset["energy"])
    dataset.nacs.values = tmp_full_nacs
    mark_variable_assigned(dataset["nacs"])
    # dataset = dataset.assign_coords(time=tmp_times)

    return actual_max_ts + 1, dataset  # , dataset

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
