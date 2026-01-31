from io import TextIOWrapper
import pathlib
from typing import Any, NamedTuple, Tuple
import numpy as np
from shnitsel.io.shared.trajectory_setup import (
    OptionalTrajectorySettings,
    RequiredTrajectorySettings,
    assign_optional_settings,
    assign_required_settings,
    create_initial_dataset,
)
from shnitsel.io.shared.variable_flagging import (
    is_variable_assigned,
    mark_variable_assigned,
)
import xarray as xr
import logging
import os
import re
import math

from shnitsel.units.definitions import length

from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
from shnitsel.data.atom_helpers import get_atom_number_from_symbol

from ..xyz import parse_xyz


# TODO: FIXME: use loading_parameters to configure units and state names
def parse_newtonx(
    traj_path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read a NewtonX trajectory directory into a Dataset with standard shnitsel annotations and units

    Parameters
    ----------
    pathlist : PathOptionsType
        Path to the NewtonX trajectory output
    loading_parameters : LoadingParameters | None, optional
        Parameter settings for e.g. standard units or state names.

    Returns
    -------
    xr.Dataset
        The Dataset object containing all of the loaded data in default shnitsel units
    """

    path_obj: pathlib.Path = make_uniform_path(traj_path)  # type: ignore

    misc_settings, newtonx_res = parse_nx_misc_input_settings(path_obj)
    actual_steps = -1

    # TODO: FIXME: Use other available files to read input if nx.log is not available, e.g. RESULTS/dyn.out, RESULTS/dyn.xyz, RESULTS/en.dat From those we can get most information anyway.
    # TODO: FIXME: Read basis from JOB_NAD files
    # Add time dimension
    if (
        newtonx_res.num_states <= 0
        or newtonx_res.num_atoms <= 0
        or newtonx_res.num_steps <= 0
    ):
        logging.error(
            "Could not extract trajectory settings like number of atoms, number of states or number of steps."
        )
        raise FileNotFoundError(
            f"Could not find key files like `nx.log` in NewtonX directory:{path_obj}"
        )

    trajectory, default_format_attributes = create_initial_dataset(
        newtonx_res.num_steps,
        newtonx_res.num_states,
        newtonx_res.num_atoms,
        "newtonx",
        loading_parameters,
    )

    nxlog_path = path_obj / "RESULTS" / "nx.log"
    moldynlog_path = path_obj / "moldyn.log"
    if nxlog_path.is_file() or moldynlog_path.is_file():
        # Switch to one existing file
        if nxlog_path.is_file():
            used_nx_path = nxlog_path
        else:
            used_nx_path = moldynlog_path
        with open(used_nx_path) as f:
            # Read several datasets into trajectory and get first indicator of actual performed steps
            actual_steps, trajectory = parse_nx_log_data(
                f, trajectory, newtonx_res, default_format_attributes
            )

            if actual_steps < newtonx_res.num_steps:
                # Filter only assigned timesteps
                trajectory.attrs["completed"] = False
                trajectory = trajectory.isel(time=slice(0, actual_steps))
            elif actual_steps > newtonx_res.num_steps:
                raise ValueError(
                    f"Trajectory data at {path_obj} contained data for {actual_steps} frames, which is more than the initially denoted {newtonx_res.num_steps} frames. Cannot allocate space after the fact."
                )
    else:
        logging.warning(
            "Missing %(path)s for input of comprehensive settings, newtonx version and observables like forces and NACs.",
            {'path': nxlog_path},
        )

    dyn_xyz_path = path_obj / "RESULTS" / "dyn.xyz"
    if dyn_xyz_path.is_file():
        with open(dyn_xyz_path) as f:
            atNames, atNums, atXYZ = parse_xyz(f)
            trajectory.atNames.values = atNames
            mark_variable_assigned(trajectory["atNames"])
            trajectory.atNums.values = atNums
            mark_variable_assigned(trajectory["atNums"])
            trajectory.atXYZ.values = atXYZ
            trajectory.atXYZ.attrs['units'] = length.Bohr  # noqa: F821
            mark_variable_assigned(trajectory["atXYZ"])
    else:
        logging.info(
            "Did not find %(path)s for position, atom names and types input.",
            {'path': dyn_xyz_path},
        )

    en_dat_path = path_obj / "RESULTS" / "en.dat"

    if en_dat_path.is_file():
        try:
            trajectory = parse_en_data(
                en_dat_path, trajectory, default_format_attributes
            )

        except FileNotFoundError:
            logging.warning(
                "Could not read %(path)s for E_kin, times and E_pot<state> input.",
                {'path': en_dat_path},
            )
            pass
    else:
        logging.info(
            "Did not find %(path)s for E_kin, times and E_pot<state> input.",
            {'path': en_dat_path},
        )

    dyn_out_path = path_obj / "RESULTS" / "dyn.out"
    if dyn_out_path.is_file():
        with open(dyn_out_path) as f:
            trajectory = parse_dyn_out(f, trajectory)
    else:
        logging.info(
            "Did not find %(path)s for E_kin and velocities input.",
            {'path': dyn_out_path},
        )
        pass

    if not is_variable_assigned(trajectory["time"]):
        # Time is not initilized with a variable, hence we need to apply default attributes here.
        trajectory = trajectory.assign_coords(
            {
                "time": (
                    "time",
                    np.arange(0, newtonx_res.num_steps) * newtonx_res.delta_t,
                    default_format_attributes[str("time")],
                ),
            }
        )
        mark_variable_assigned(trajectory["time"])

    # Set all settings we require to be present on the trajectory
    # TODO: FIXME: Check if we can actually derive the number of singlets, doublets, or triplets from newtonx output.
    required_settings = RequiredTrajectorySettings(
        newtonx_res.t_max,
        newtonx_res.delta_t,
        min(newtonx_res.num_steps, actual_steps)
        if actual_steps > 0
        else newtonx_res.num_steps,
        newtonx_res.completed,
        "newtonx",
        "dynamic",
        newtonx_res.newtonx_version,
        newtonx_res.num_states,
        0,
        0,
    )
    trajectory = assign_required_settings(trajectory, required_settings)

    # TODO: FIXME: Check if newtonx always prints only the active state forces or sometimes may include other forces.
    optional_settings = OptionalTrajectorySettings(
        has_forces="active_only"
        if is_variable_assigned(trajectory["forces"])
        else None,
        misc_input_settings=misc_settings,
    )
    trajectory = assign_optional_settings(trajectory, optional_settings)

    return trajectory


class NewtonXSettingsResult(NamedTuple):
    """Class to keep track of key settings from the NewtonX log files"""

    t_max: float
    delta_t: float
    num_steps: int
    num_atoms: int
    num_states: int
    completed: bool
    newtonx_version: str


def parse_settings_from_nx_log(
    f: TextIOWrapper,
) -> Tuple[NewtonXSettingsResult, dict[str, Any]]:
    """Function to parse key settings from the NewtonX RESULTS/nx.log file

    Parameters
    ----------
    f : TextIOWrapper
        The input stream of lines from nx.log

    Returns
    -------
    NewtonXSettingsResult
        Key setting parameters from nx.log
    dict[str, Any]
        A Dict of all relevant settings loaded from nx.log
    """
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

    logging.debug("found real_tmax: %(tmax)f", {'tmax': real_tmax})
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

    return NewtonXSettingsResult(
        real_tmax, delta_t, nsteps, natoms, nstates, completed, newtonx_version
    ), settings


def parse_en_data(
    endata_path: pathlib.Path, dataset: xr.Dataset, default_attributes: dict
) -> xr.Dataset:
    """Function to read energy data from en.dat in case the reading from nx.log has not succeeded:

    Parameters
    ----------
    endata_path : pathlib.Path
        Path to a RESULTS/en.dat file
    dataset : xr.Dataset
        The dataset to update
    default_attributes : dict
        Default attributes for the newtonx format

    Returns
    -------
    xr.Dataset
        the updated dataset
    """

    endata = np.loadtxt(endata_path, ndmin=2)

    times = endata[:, 0]
    en_tot = endata[:, -1]
    en_active = endata[:, -2]
    en_states = endata[:, 1:-2]

    e_kin = en_tot - en_active

    # print(default_attributes["time"])

    if not is_variable_assigned(dataset["time"]):
        logging.debug("Assigning time from en.dat")
        dataset = dataset.assign_coords(
            {
                "time": (
                    "time",
                    times,
                    default_attributes["time"],
                ),
            }
        )
        mark_variable_assigned(dataset["time"])
    if not is_variable_assigned(dataset["energy"]):
        logging.debug("Assigning energy from en.dat")
        dataset["energy"].values = en_states
        mark_variable_assigned(dataset["energy"])
    if not is_variable_assigned(dataset["e_kin"]):
        logging.debug("Assigning e_kin from en.dat")
        dataset["e_kin"].values = e_kin
        mark_variable_assigned(dataset["e_kin"])
    return dataset


def parse_dyn_out(f: TextIOWrapper, dataset: xr.Dataset) -> xr.Dataset:
    """Function to gather dynamics data (e_kin and velocities) from the RESULTS/dyn.out file, if it exists.

    Parameters
    ----------
    f :TextIOWrapper
        The stream of lines from the dyn.out file
    dataset : xr.Dataset
        The dataset to write the data to

    Returns
    -------
    xr.Dataset
        The updated dataset after the read
    """

    natoms = dataset.sizes["atom"]
    nstates = dataset.sizes["state"]
    ntimesteps = dataset.sizes["time"]

    ts: int = 0

    tmp_pos_in_bohr = np.zeros((ntimesteps, natoms, 3))
    tmp_e_kin = np.zeros(
        (
            ntimesteps,
            nstates,
        )
    )
    tmp_velocities_in_au = np.zeros((ntimesteps, natoms, 3))

    tmp_e_pot = np.zeros(
        (
            ntimesteps,
            nstates,
        )
    )

    tmp_atNames = []
    tmp_atMasses = []
    has_atNames = False

    tmp_astate = np.full((ntimesteps,), -1)

    tmp_forces_in_au = np.zeros((ntimesteps, nstates, natoms, 3))

    has_positions = False
    has_velocities = False
    has_forces = False
    has_ekin = False
    has_epot = False
    has_astate = False

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
        if stripline.startswith("STEP"):
            parts = stripline.split()
            # Figure out current time step
            ts = int(parts[1].strip())

            state_prefix = "on state "
            prefix_offset = stripline.find(state_prefix)
            if prefix_offset > 0:
                state_parts = stripline[prefix_offset + len(state_prefix) :].split()
                astate = int(state_parts[0])
                has_astate = True
            else:
                astate = -1

            tmp_astate[ts] = astate

            # logging.debug(f"Starting ts {ts}")

        elif stripline.find("geometry:") > 0:
            has_positions = True
            for iatom in range(natoms):
                line_parts = next(f).strip().split()
                atName = line_parts[0]
                tmp_pos_in_bohr[ts][iatom] = [float(n) for n in line_parts[2:-1]]
                tmp_mass = float(line_parts[-1])

                if not has_atNames:
                    tmp_atNames.append(atName)
                    tmp_atMasses.append(tmp_mass)

            has_atNames = True
        elif stripline.find("velocity:") > 0:
            has_velocities = True
            for iatom in range(natoms):
                tmp_velocities_in_au[ts][iatom] = [
                    float(n) for n in next(f).strip().split()
                ]
        elif stripline.find("acceleration:") > 0:
            if not has_atNames:
                logging.error(
                    "Found acceleration data before finding atom mass information. Input may be corrupted."
                )
            if not tmp_astate[ts] > 0:
                logging.error(
                    "Found acceleration data before finding active state information. Input may be corrupted."
                )
            else:
                has_forces = True
                for iatom in range(natoms):
                    tmp_forces_in_au[ts][tmp_astate[ts] - 1][iatom] = [
                        float(n) for n in next(f).strip().split()
                    ]
                    tmp_forces_in_au[ts][tmp_astate[ts] - 1][iatom] *= tmp_atMasses[
                        iatom
                    ]

        elif stripline.startswith("Time") and stripline.find("Ekin:") > 0:
            #     Time    Etot         Ekin       Epot E0,      E1, ...
            # %       0.50   -232.354045      0.055764   -232.578726   -232.409809   -232.393353
            next_line = next(f).strip()
            next_line_parts = next_line.split()

            tmp_e_kin[ts] = float(next_line_parts[2])
            tmp_e_pot[ts] = [float(next_line_parts[3 + i]) for i in range(nstates)]

            has_ekin = True
            has_epot = True

    if has_astate and not is_variable_assigned(dataset["astate"]):
        # logging.info("Currently, the E_pot from dyn.out is unused")
        logging.debug("Assigning astate from dyn.out")
        dataset.astate.values = tmp_astate
        mark_variable_assigned(dataset["astate"])

    if has_epot and not is_variable_assigned(dataset["energy"]):
        # logging.info("Currently, the E_pot from dyn.out is unused")
        logging.debug("Assigning energy from dyn.out")
        dataset.energy.values = tmp_e_pot
        mark_variable_assigned(dataset["energy"])

    if has_ekin and not is_variable_assigned(dataset["e_kin"]):
        logging.debug("Assigning e_kin from dyn.out")
        dataset.e_kin.values = tmp_e_kin
        mark_variable_assigned(dataset["e_kin"])

    if has_velocities and not is_variable_assigned(dataset["velocities"]):
        logging.info("Assigning velocities from dyn.out")
        dataset.velocities.values = tmp_velocities_in_au
        mark_variable_assigned(dataset["velocities"])

    if has_forces and not is_variable_assigned(dataset["forces"]):
        logging.info("Assigning forces from dyn.out")
        dataset.forces.values = tmp_forces_in_au
        mark_variable_assigned(dataset["forces"])

    if has_positions and not is_variable_assigned(dataset["atXYZ"]):
        # If no other option to set positions has been used, use this source
        from shnitsel.units.definitions import length

        logging.debug("Assigning positions from dyn.out")
        dataset["atXYZ"].values = tmp_pos_in_bohr
        dataset["atXYZ"].attrs["unit"] = length.Bohr
        mark_variable_assigned(dataset["atXYZ"])

    if has_atNames and not is_variable_assigned(dataset["atNames"]):
        tmp_atNums = [get_atom_number_from_symbol(sym) for sym in tmp_atNames]
        dataset["atNames"].values = tmp_atNames
        dataset["atNums"].values = tmp_atNums
        logging.debug("Assigning atom names and types from dyn.out")
        mark_variable_assigned(dataset["atNames"])
        mark_variable_assigned(dataset["atNums"])

    return dataset


def parse_nx_log_data(
    f: TextIOWrapper,
    dataset: xr.Dataset,
    settings: NewtonXSettingsResult,
    default_attributes: dict,
) -> Tuple[int, xr.Dataset]:  # Tuple[int, xr.Dataset]:
    """Function to parse the nx.log data into a dataset from the input stream f.

    Will return the total number of actual timesteps read and the resulting dataset.
    Usual read data includes: forces, active state ("astate")

    Parameters
    ----------
    f : TextIOWrapper
        Input filestream of a nx.log file
    dataset : xr.Dataset
        The dataset to parse the data into
    settings: NewtonXSettingsResult
        The previously parsed settings
    default_attributes :  dict
        Default attributes for the newtonx format

    Returns
    -------
    Tuple[int, xr.Dataset]
        The total number of actual timesteps read and the resulting dataset after applying modifications.
    """
    # f should be after the setting

    natoms = dataset.sizes["atom"]
    nstates = dataset.sizes["state"]
    nstatecomb = dataset.sizes["statecomb"]
    ntimesteps = dataset.sizes["time"]

    actual_max_ts: int = -1

    delta_t = settings.delta_t

    ts: int = 0
    time: float = 0

    tmp_astate = np.zeros((ntimesteps,), dtype=int)
    tmp_times = np.zeros((ntimesteps,))
    tmp_forces = np.zeros((natoms, 3))
    tmp_energy = np.zeros((nstates,))
    tmp_nacs = np.zeros((nstatecomb, natoms, 3))
    tmp_full_forces = np.zeros((ntimesteps, nstates, natoms, 3))
    tmp_full_energy = np.zeros(
        (
            ntimesteps,
            nstates,
        )
    )
    tmp_full_nacs = np.zeros((ntimesteps, nstatecomb, natoms, 3))

    full_has_energy = False
    full_has_forces = False
    full_has_nacs = False

    step_has_energy = False
    step_has_forces = False
    step_has_nacs = False
    step_has_nacs_norm = False

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
            tmp_times[ts] = t_time
            # Assign this to only the active state
            if step_has_forces:
                tmp_full_forces[ts][t_astate - 1] = tmp_forces
                full_has_forces = True
            if step_has_energy:
                tmp_full_energy[ts] = tmp_energy
                full_has_energy = True
            if step_has_nacs:
                tmp_full_nacs[ts] = tmp_nacs
                full_has_nacs = True

            step_has_energy = False
            step_has_forces = False
            step_has_nacs = False

            actual_max_ts = max(actual_max_ts, ts)
            # logging.debug(f"finished ts {ts}")

        elif stripline.lower().startswith("gradient vectors"):
            step_has_forces = True
            for iatom in range(natoms):
                tmp_forces[iatom] = [float(n) for n in next(f).strip().split()]

        elif stripline.lower().startswith("nonadiabatic coupling vectors"):
            # NOTE: The label for full entries is nonadiabatic coupling vectors
            # There also exist entries for nonadiabatic coupling terms that are the normed full vectors across all atoms.

            step_has_nacs = True
            for icomb in range(math.comb(nstates, 2)):
                # Order is: V(from, to),iatom, dir
                # Increase steps from rightmost dimension to leftmost.
                # I.e. each line has x,y,z value.
                # First natoms lines have the values for the different atoms in first state combination
                # Each block of natoms represents successive state combinations
                for iatom in range(natoms):
                    tmp_nacs[icomb, iatom] = [float(n) for n in next(f).strip().split()]
        elif stripline.lower().startswith("nonadiabatic coupling terms"):
            # TODO: FIXME: consider allowing to load nacs_norm directly from input instead.
            if not step_has_nacs_norm:
                step_has_nacs_norm = True
                logging.warning(
                    "We currently do not support reading pre-normalized/scalarized NACs terms (i.e. non-vectors) in our NewtonX input."
                )
        elif stripline.lower().startswith("energy ="):
            for istate in range(nstates):
                tmp_energy[istate] = float(next(f).strip())
            step_has_energy = True

    if not is_variable_assigned(dataset.astate):
        dataset["astate"].values = tmp_astate
        mark_variable_assigned(dataset["astate"])

    if full_has_forces and not is_variable_assigned(dataset.forces):
        logging.debug("Assigning forces from nx.log")
        dataset.forces.values = tmp_full_forces
        mark_variable_assigned(dataset["forces"])
    if full_has_energy and not is_variable_assigned(dataset.energy):
        logging.debug("Assigning energy from nx.log")
        dataset.energy.values = tmp_full_energy
        mark_variable_assigned(dataset["energy"])
    if full_has_nacs and not is_variable_assigned(dataset.nacs):
        logging.debug("Assigning nacs from nx.log")
        dataset.nacs.values = tmp_full_nacs
        mark_variable_assigned(dataset["nacs"])

    if not is_variable_assigned(dataset["time"]):
        dataset = dataset.assign_coords(
            {"time": ("time", tmp_times, default_attributes["time"])}
        )
        mark_variable_assigned(dataset["time"])

    logging.debug("Finish reading NewtonX log.")

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


def parse_nx_misc_input_settings(
    path: pathlib.Path,
) -> tuple[dict, NewtonXSettingsResult]:
    """Function to parse various input settings from the newtonx trajectory directory.

    Parameters
    ----------
    path : pathlib.Path
        The path of the base trajectory folder. Should Contain `RESULTS`, `control.dyn` and either `JOB_AD` or `JOB_NAD`

    Returns
    -------
    dict
        The collected miscallenous settings. For specific settings, the keys "nx.log", "control.dyn"
        allow for searching through specific settings if the files were present.
    NewtonXSettingsResult
        The combined key settings if they could be extracted.
    """

    res = {}

    extracted_dt = -1.0
    extracted_tmax = -1.0
    extracted_nsteps = -1
    extracted_natoms = -1
    extracted_nstates = -1
    extracted_completed = False
    extracted_newtonx_version = "unknown"

    control_dyn_path = path / "control.dyn"
    if control_dyn_path.is_file():
        control_dyn_settings = {}
        with open(control_dyn_path) as dyn:
            lines = [l.strip() for l in dyn.readlines()]

            for l in lines:
                if l.find("=") >= 0:
                    key, val = re.split(" *= *", l.strip(), maxsplit=1)
                    key = key.strip()
                    val = val.strip()
                    if key != "prog":
                        if "." in val:
                            val = float(val)
                        else:
                            val = int(val)
                    else:
                        prog_val = float(val)
                        if prog_val >= 13.0:
                            logging.warning(
                                "We currently do not support reading Triplet states from NewtonX trajectories. \n"
                                "If you require reading of triplet state input from NewtonX format, please open an issue in the Shnitsel-Tools github repo"
                                " and provide an example configuration/setup.\n"
                                "For now, all state multiplicities will be assumed to be singlets."
                            )
                    control_dyn_settings[key] = val

        res["control.dyn"] = control_dyn_settings

    endata_path = path / "RESULTS" / "en.dat"

    if endata_path.is_file():
        endata_res = np.loadtxt(endata_path, ndmin=2)

        if endata_res is not None and len(endata_res) > 0:
            times = endata_res[:, 0]
            extracted_nsteps = max(extracted_nsteps, endata_res.shape[0])
            extracted_tmax = max(extracted_tmax, np.max(times))
            extracted_dt = times[1] - times[0]
            extracted_nstates = (
                endata_res.shape[1] - 3
            )  # No times, active energy and total energy

    dyn_xyz_path = path / "RESULTS" / "dyn.xyz"
    if dyn_xyz_path.is_file():
        with open(dyn_xyz_path) as f:
            while len((line := f.readline()).strip()) == 0:
                pass

            # Number of atoms on first line
            extracted_natoms = int(line)

    veloc_path = path / "veloc"
    if veloc_path.is_file():
        veloc_data = np.loadtxt(veloc_path, ndmin=2)
        # Number of atoms on first line
        extracted_natoms = veloc_data.shape[0]

    # TODO: FIXME: Add input parsing from QM directories?
    jobnad_path = path / "JOB_NAD"
    jobad_path = path / "JOB_AD"

    results_path_log = path / "RESULTS" / "nx.log"
    if results_path_log.is_file():
        with open(results_path_log) as f:
            newtonx_res, misc_settings = parse_settings_from_nx_log(f)
            res["nx.log"] = misc_settings

            extracted_dt = newtonx_res.delta_t
            extracted_tmax = max(extracted_nsteps, newtonx_res.t_max)
            extracted_nsteps = max(extracted_nsteps, newtonx_res.num_steps)
            extracted_natoms = newtonx_res.num_atoms
            extracted_nstates = newtonx_res.num_states
            extracted_completed = newtonx_res.completed
            extracted_newtonx_version = newtonx_res.newtonx_version
    else:
        logging.warning(
            "The NewtonX log file at %(path)s is missing. Key settings, gradients, NACs, etc. will not be extracted.",
            {'path': results_path_log},
        )

    return res, NewtonXSettingsResult(
        extracted_tmax,
        extracted_dt,
        extracted_nsteps,
        extracted_natoms,
        extracted_nstates,
        extracted_completed,
        extracted_newtonx_version,
    )
