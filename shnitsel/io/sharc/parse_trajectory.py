from io import TextIOWrapper
import pathlib
from typing import Any, Dict, List, Tuple
import numpy as np
import xarray as xr
from itertools import combinations, permutations
import pandas as pd
import logging
import os
import re

from shnitsel.io.shared.helpers import (
    PathOptionsType,
    make_uniform_path,
)
from shnitsel.data.atom_helpers import get_atom_number_from_symbol
from shnitsel.io.sharc.qm_helpers import (
    INTERFACE_READERS,
)
from shnitsel.data.state_helpers import set_sharc_state_type_and_name_defaults

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
from shnitsel.units.definitions import (
    distance,
    _distance_unit_scales,
)


from shnitsel.io.shared.helpers import LoadingParameters


def read_traj(
    traj_path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read a single SHARC trajectory directory

    Parameters
    ----------
    traj_path : PathOptionsType
        The path to load the trajectory form
    loading_parameters : LoadingParameters | None, optional
        Parameter settings for e.g. standard units or state names.

    Returns
    -------
    xr.Dataset
        The parsed SHARC directory as a Dataset
    """
    # TODO: FIXME: use loading_parameters to configure units and state names

    path_obj: pathlib.Path = make_uniform_path(traj_path)

    # Read some settings from input
    # In particular, if a trajectory is extended by increasing
    # tmax and resuming, the header of output.dat will give
    # only the original nsteps, leading to an ndarray IndexError
    input_path = path_obj / "input"
    output_path = path_obj / "output.dat"
    output_listing_path = path_obj / "output.lis"
    output_log_path = path_obj / "output.log"
    geom_file = path_obj / "geom"
    veloc_file = path_obj / "veloc"

    sharc_version = "unknown"
    nsinglets: int = 0
    ndoublets: int = 0
    ntriplets: int = 0

    state_multiplicities = None
    state_charges = None

    delta_t: float | None = None
    t_max: float | None = None
    nsteps: int = 0
    natoms: int | None = None
    nstates: int | None = None
    energy_offset: float | None = None
    variables_listings = None
    completed = True

    misc_settings = {}

    if input_path.is_file():
        with open(input_path) as f:
            settings = parse_input_settings(f.readlines())
        delta_t = float(settings["stepsize"])
        t_max = float(settings["tmax"])
        nsteps = max(nsteps, int(round(t_max / delta_t)) + 1)
        energy_offset = float(settings["ezero"])

        if "nstates" in settings:
            state_mult_array = [int(x.strip()) for x in settings["nstates"].split()]
            settings["nstates"] = state_mult_array
            max_mult = len(state_mult_array)

            nsinglets = state_mult_array[0] if max_mult >= 1 else 0
            ndoublets = state_mult_array[1] if max_mult >= 2 else 0
            ntriplets = state_mult_array[2] if max_mult >= 3 else 0
            state_multiplicities = state_mult_array
        if "charge" in settings:
            charge_settings = settings["charge"]
            logging.debug(
                "Found charge info in output data: %(charge_settings)s",
                {'charge_settings': charge_settings},
            )
            state_charges = [int(x.strip()) for x in charge_settings.split()]

        misc_settings["input"] = settings

    if output_path.is_file():
        with open(output_path) as f:
            settings = parse_output_settings(f)
        # Unit of dtstep is completely unclear.
        # delta_t = float(settings["dtstep"]) *time_

        nsteps = max(nsteps, int(settings["nsteps"]) + 1)
        natoms = int(settings["natom"])

        energy_offset = settings["ezero"]
        sharc_version = settings["SHARC_version"]

        if "nstates" in settings:
            state_mult_array = [int(x.strip()) for x in settings["nstates"].split()]
            max_mult = len(state_mult_array)

            nsinglets = state_mult_array[0] if max_mult >= 1 else 0
            ndoublets = state_mult_array[1] if max_mult >= 2 else 0
            ntriplets = state_mult_array[2] if max_mult >= 3 else 0
            state_multiplicities = state_mult_array
        if "charge" in settings:
            charge_settings = settings["charge"]
            logging.debug(
                "Found charge info in output data: %(charge_settings)s",
                {'charge_settings': charge_settings},
            )
            state_charges = [int(x.strip()) for x in charge_settings.split()]
        misc_settings["output.dat"] = settings

    if output_listing_path.is_file():
        settings, variables_listings = parse_output_listings(output_listing_path)
        delta_t = float(settings["delta_t"])
        t_max = float(settings["t_max"])
        nsteps = max(nsteps, int(settings["nsteps"]) + 1)
        misc_settings["output.lis"] = settings

    if output_log_path.is_file():
        with open(output_log_path) as f:
            settings = parse_output_log(f)

            if "version" in settings:
                sharc_version = settings["version"]
            if "SHARC_version" in settings:
                sharc_version = settings["SHARC_version"]

            delta_t = float(settings["stepsize"])
            t_max = float(settings["tmax"])
            nsteps = int(t_max / delta_t) + 1
            energy_offset = settings["ezero"]

            if "nstates" in settings:
                state_mult_array = [int(x.strip()) for x in settings["nstates"].split()]
                max_mult = len(state_mult_array)

                nsinglets = state_mult_array[0] if max_mult >= 1 else 0
                ndoublets = state_mult_array[1] if max_mult >= 2 else 0
                ntriplets = state_mult_array[2] if max_mult >= 3 else 0
        misc_settings["output.log"] = settings

    if energy_offset is None:
        raise FileNotFoundError(
            f"Could not detect ezero offset for SHARC data from path {traj_path}. Make sure, output.dat, input or output.log are present."
        )

    # TODO: FIXME: Check if the factors 1/2/3 are correct or if we should just sum up the states?
    # In principle, there are 2 states per doublet and three per triplet. However, not all will be used.
    # See SHARC documentation for state order
    nstates = nsinglets + ndoublets * 2 + ntriplets * 3

    # Try other sources for the number of atoms
    if natoms is None:
        if geom_file.is_file():
            with open(geom_file) as f:
                natoms = int(
                    round(
                        np.sum([1 if len(x.strip()) > 0 else 0 for x in f.readlines()])
                    )
                )
    if natoms is None:
        if veloc_file.is_file():
            with open(veloc_file) as f:
                natoms = int(
                    round(
                        np.sum([1 if len(x.strip()) > 0 else 0 for x in f.readlines()])
                    )
                )

    if nsteps is None or natoms is None:
        raise FileNotFoundError(
            "Could not find enough information to deduce the number of atoms or steps."
        )

    trajectory, default_format_attributes = create_initial_dataset(
        nsteps, nstates, natoms, "sharc", loading_parameters
    )

    if output_path.is_file():
        # Try and parse full output from output.dat
        with open(os.path.join(traj_path, "output.dat")) as f:
            completed_dat, max_ts_dat, trajectory = parse_trajout_dat(
                f, trajectory_in=trajectory, loading_parameters=loading_parameters
            )

        completed = completed and completed_dat

        # Filter out only the assigned part if the trajectory did not complete
        if not completed:
            # Filter by index not by ts
            # res = res.sel(ts=res.ts <= res.attrs["max_ts"])
            trajectory = trajectory.isel(time=slice(0, max_ts_dat))

        nsteps = min(trajectory.sizes["time"], max_ts_dat)

    # Starting with SHARC 4.1, there will be a "start.time" file with a time offset
    start_time_file = path_obj / "start.time"
    time_offset = 0.0
    if start_time_file.exists():
        with open(start_time_file) as f:
            lines = f.readlines()
            time_offset = float(lines[0])

    if variables_listings is not None:
        if nsteps > len(variables_listings["time"]):
            completed = False

        if not is_variable_assigned(trajectory.time):
            if "time" in variables_listings:
                # logging.debug(f"Time attributes: {trajectory.time.attrs}")
                trajectory.coords["time"] = (
                    "time",
                    variables_listings["time"] + time_offset,
                    default_format_attributes["time"],
                )
                mark_variable_assigned(trajectory["time"])
        if not is_variable_assigned(trajectory.astate):
            if "astate" in variables_listings:
                trajectory.coords["astate"] = (
                    "time",
                    variables_listings["astate"],
                    default_format_attributes["astate"],
                )
                mark_variable_assigned(trajectory["astate"])

    # TODO: Note that for consistency, we renamed the ts dimension to time to agree with other format
    if (
        not is_variable_assigned(trajectory.atNames)
        or not is_variable_assigned(trajectory.atNums)
        or not is_variable_assigned(trajectory.atXYZ)
    ):
        with open(os.path.join(traj_path, "output.xyz")) as f:
            atNames, atNums, atXYZ = parse_trajout_xyz(nsteps, f)

            trajectory.coords["atNames"] = ("atom", atNames, trajectory.atNames.attrs)
            mark_variable_assigned(trajectory.atNames)

            trajectory.coords["atNums"] = ("atom", atNums, trajectory.atNums.attrs)
            mark_variable_assigned(trajectory.atNums)

            # These are in Angstrom for some weird reason... why not consistent? Who knows.
            trajectory["atXYZ"][...] = (
                atXYZ
                * _distance_unit_scales[distance.Angstrom]
                / _distance_unit_scales[distance.Bohr]
            )
            mark_variable_assigned(trajectory.atXYZ)

    # Set all settings we require to be present on the trajectory
    required_settings = RequiredTrajectorySettings(
        t_max,
        delta_t,
        nsteps,
        completed,
        "sharc",
        "dynamic",
        sharc_version,
        nsinglets,
        ndoublets
        * 2,  # To agree with numbering conventions of the other systems, we need to multiply by 2.
        ntriplets
        * 3,  # To agree with numbering conventions of the other systems, we need to multiply by 3.
    )

    trajectory = assign_required_settings(trajectory, required_settings)

    optional_settings = OptionalTrajectorySettings(
        has_forces=is_variable_assigned(trajectory["forces"]),
        misc_input_settings=misc_settings if len(misc_settings) > 0 else None,
    )

    if state_multiplicities is not None:
        charges = state_charges
        if charges is None:
            logging.debug("Missing charge information.")
            main_version = int(sharc_version.split(".")[0])
            if main_version < 4:
                logging.info(
                    "For sharc before version 4.0, we will attempt to extract charge data from QM interface settings."
                )

                qm_path = path_obj / "QM"
                for int_name, int_reader in INTERFACE_READERS.items():
                    # logging.debug(f"Trying format: {int_name}")
                    res_dict = int_reader(trajectory, qm_path)
                    # logging.debug(f"Res qm data: {res_dict}")

                    if "theory_basis" in res_dict:
                        optional_settings.theory_basis_set = res_dict["theory_basis"]
                    if "est_level" in res_dict:
                        optional_settings.est_level = res_dict["est_level"]
                    if "charge" in res_dict:
                        charges = res_dict["charge"]

                    if charges is not None:
                        logging.info(
                            "Found charge data from the %(int_name)s interface",
                            {'int_name': int_name},
                        )
                        break
            else:
                # Assume we are uncharged if no charge data found.
                logging.info(
                    "We assume there is no charge because no charge information was found"
                )
                charges = 0
        else:
            logging.debug("Charge information present.")

        trajectory = set_sharc_state_type_and_name_defaults(
            trajectory, state_multiplicities, charges
        )
    trajectory =assign_optional_settings(trajectory, optional_settings)

    return trajectory


def parse_output_settings(f: TextIOWrapper) -> dict[str, Any]:
    """Function to parse settings from the `output.dat` file as far as they are available.

    Read settings and other info from the file.

    Parameters
    ----------
    f : TextIOWrapper
        File wrapper providing the `output.dat` file contents

    Returns
    -------
    dict[str, Any]
        A key-value dictionary, where the keys are the names of the settings
    """
    settings = {}
    for line in f:
        stripped = line.strip()

        if stripped.startswith("****"):
            break
        parsed = stripped.split()

        if len(parsed) == 2:
            settings[parsed[0]] = parsed[1]
        elif len(parsed) > 2:
            settings[parsed[0]] = " ".join(parsed[1:])
        elif len(parsed) == 1:
            settings[parsed[0]] = True

    # We have reached the end of the settings block:

    return settings


def parse_output_listings(path: pathlib.Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Function to parse settings from the `output.lis` file as far as they are available.

    This is used to read the delta_t variable if it hasn't been set otherwise.

    Parameters
    ----------
    path : pathlib.Path
        Path to the `output.lis` file.

    Returns
    -------
    tuple[dict[str, Any], dict[str, Any]]
        First, a key-value dictionary, where the keys are the names of the settings like delta_t, nsteps and t_max. 
        Then a key_value dictionary with names of variables and their values extracted from the file.
    """
    settings = {}
    variables = {}

    lis_data = np.loadtxt(path, ndmin=2)

    steps = lis_data[:, 0]
    nsteps = int(round(np.max(steps))) + 1
    settings["nsteps"] = nsteps

    times = lis_data[:, 1]
    settings["delta_t"] = times[1] - times[0]
    settings["t_max"] = np.max(times)

    active_state = np.array([int(round(x)) for x in lis_data[:, 2]])

    epot_relative_active = lis_data[:, 5]

    variables["astate"] = active_state
    variables["time"] = times
    variables["active_state_E"] = epot_relative_active

    return settings, variables


def parse_output_log(f: TextIOWrapper) -> Dict[str, Any]:
    """Function to parse settings from the `output.log` file as far as they are available.

    This is used to read the version of sharc and input settings if unavailable elsewhere.

    Parameters
    ----------
    f : TextIOWrapper
        The file input stream from `output.log`

    Returns
    -------
    dict[str, Any]
        A key-value dictionary, where the keys are the names of the settings like delta_t, t_max, or version.
    """
    settings = {}

    sharc_version = None

    for line in f:
        stripped = line.strip()

        if stripped.startswith("Version:"):
            sharc_version = stripped.split()[1]

        if stripped.startswith("Input File"):
            break

    # Skip underline
    line = next(f)

    input_lines = []
    for line in f:
        line_stripped = line.strip()
        if line_stripped.startswith("====="):
            break
        if len(line_stripped) > 0:
            input_lines.append(line_stripped)

    settings = parse_input_settings(input_lines)

    if sharc_version is not None:
        settings["version"] = sharc_version

    return settings


# TODO: FIXME: The format differs between versions of SHARC. 2.0, 3.0 and 4.0 have different meta-information about the referenced state in their output.dat.
# TODO: FIXME: They also have different multipliers and naming convention than the other formats
def parse_trajout_dat(
    f: TextIOWrapper,
    trajectory_in: xr.Dataset,
    loading_parameters: LoadingParameters | None = None,
) -> Tuple[bool, int, xr.Dataset]:
    """Function to parse the contents of an 'output.dat' in a sharc trajectory output directory into a Dataset.

    Parameters
    ----------
    f : TextIOWrapper
        A file wrapper providing the contents of 'output.dat'.
    nsteps : int | None, optional
        The number of maximum steps expected. Defaults to None.

    Raises
    ------
    ValueError
        Raised if not enough steps are found in the output.dat file

    Returns
    -------
    xr.Dataset
        A flag to indicate if the full trajectory has been read, the number of steps that acutally were read and the full dataset with unit attributes and further helpful attributes applied.
    """
    settings = parse_output_settings(f)
    nsteps = trajectory_in.sizes["time"]

    nsteps_output_dat = int(settings["nsteps"]) + 1  # let's not forget ts=0
    # logging.debug(f"(From input file) nsteps = {nsteps}")
    nsteps = max(nsteps_output_dat, nsteps)

    natoms = int(settings["natom"])  # yes, really 'natom', not 'natoms'!
    # logging.debug(f"natoms = {natoms}")
    energy_offset_zero = float(settings["ezero"])
    # logging.debug(f"energy_offset_zero = {energy_offset_zero}")
    nstates = trajectory_in.sizes["state"]
    # logging.debug(f"nstates = {nstates}")
    nstates = trajectory_in.sizes["state"]

    expect_socs = False
    if "spinorbit" in settings or "nospinorbit" not in settings:
        logging.info("Expecting SOCs in SHARC")
        expect_socs = True

    # Read atomic numbers and names from file
    # ! Atomic numbers
    # 0.6000000000000E+001
    # 0.7000000000000E+001
    # 0.1000000000000E+001
    # 0.1000000000000E+001
    # 0.1000000000000E+001
    # 0.1000000000000E+001
    # ! Elements
    # C
    # N
    # H
    # H
    # H
    # H
    # ! Atomic masses
    # 0.2187466181995E+005
    # 0.2552603505759E+005
    # 0.1837143472948E+004
    # 0.1837143472948E+004
    # 0.1837143472948E+004
    # 0.1837143472948E+004
    atNames = np.full((natoms,), "")
    atNums = np.full((natoms,), "")
    for line in f:
        stripped = line.strip()

        if stripped.startswith("****"):
            break

        if stripped.startswith("! Atomic numbers"):
            for i in range(natoms):
                atNums[i] = int(round(float(next(f).strip())))
        if stripped.startswith("! Elements"):
            for i in range(natoms):
                atNames[i] = next(f).strip()

    trajectory_in.coords["atNames"] = (
        "atom",
        atNames,
        trajectory_in.coords["atNames"].attrs,
    )
    mark_variable_assigned(trajectory_in["atNames"])
    trajectory_in.coords["atNums"] = (
        "atom",
        atNums,
        trajectory_in.coords["atNums"].attrs,
    )
    mark_variable_assigned(trajectory_in["atNums"])

    idx_table_nacs = {
        (si, sj): idx
        for idx, (si, sj) in enumerate(combinations(range(1, nstates + 1), 2))
    }

    idx_table_socs = {
        (si, sj): idx
        for idx, (si, sj) in enumerate(permutations(range(1, nstates + 1), 2))
    }

    # now we know the number of steps, we can initialize the data arrays:
    ts = -1
    max_ts = -1

    energy_assigned = False
    force_assigned = False
    dipole_assigned = False
    phases_assigned = False
    e_kin_assigned = False
    sdiag_assigned = False
    astate_assigned = False
    nacs_assigned = False
    socs_assigned = False

    sharc_version_parts = [int(x) for x in settings["SHARC_version"].split(".")]
    _sharc_main_version = sharc_version_parts[0]

    tmp_dip_all = np.full((nsteps, nstates, nstates, 3), np.nan)
    tmp_energy = np.full_like(trajectory_in.energy.values, np.nan)
    tmp_forces = np.full_like(trajectory_in.forces.values, np.nan)
    tmp_phases = np.full_like(trajectory_in.phases.values, np.nan)
    tmp_e_kin = np.full_like(trajectory_in.e_kin.values, np.nan)
    tmp_sdiag = np.full_like(trajectory_in.sdiag.values, 0, dtype=np.int32)
    tmp_astate = np.full_like(trajectory_in.astate.values, 0, dtype=np.int32)
    tmp_nacs = np.full_like(trajectory_in.nacs.values, np.nan)
    tmp_socs = np.full_like(trajectory_in.socs, 0 + 0j)

    # skip through until initial step:
    for line in f:
        if line.startswith("! 0 Step"):
            ts = int(next(f).strip())
            if ts != 0:
                logging.warning(
                    "Initial timestep's index is not 0 but %(ts)d", {'ts': ts}
                )
            max_ts = max(max_ts, ts)
            break

    for index, line in enumerate(f):
        if line[0] != "!":
            continue

        if line.startswith("! 0 Step"):
            # update `ts` to current timestep #
            new_ts = int(next(f).strip())
            if new_ts != (ts or 0) + 1:
                logging.warning(
                    "Non-consecutive timesteps: %(ts)d -> %(next_ts)d",
                    {'ts': ts, 'new_ts': new_ts},
                )
            ts = new_ts
            max_ts = max(max_ts, ts)
            # logging.debug(f"timestep = {ts}")

        if line.startswith("! 1 Hamiltonian"):
            energy_assigned = True

            for istate in range(nstates):
                stripline = next(f).strip()
                float_entries = [float(x) for x in stripline.split()]

                # Energy needs to be offset by energy_offset_zero
                tmp_energy[ts, istate] = float_entries[istate * 2] + energy_offset_zero
                if expect_socs:
                    socs_assigned = True

                    for jstate in range(nstates):
                        if istate == jstate:
                            continue

                        # State id and state index off by one
                        full_comb = idx_table_socs[(istate + 1, jstate + 1)]

                        tmp_socs[ts, full_comb] = np.complex128(
                            float_entries[jstate * 2], float_entries[jstate * 2 + 1]
                        )

        if line.startswith("! 3 Dipole moments"):
            dipole_assigned = True
            direction = {"X": 0, "Y": 1, "Z": 2}[line.strip().split()[4]]
            for istate in range(nstates):
                linecont = next(f).strip().split()
                # delete every second element in list (imaginary values, all zero)
                tmp_dip_all[ts, istate, :, direction] = [
                    float(i) for i in linecont[::2]
                ]

        if line.startswith("! 4 Overlap matrix"):
            found_overlap = False
            phasevector = np.ones((nstates))

            wvoverlap = np.zeros((nstates, nstates))
            for j in range(nstates):
                linecont = next(f).strip().split()
                # delete every second element in list (imaginary values, all zero)
                wvoverlap[j] = [float(n) for n in linecont[::2]]

            for istate in range(nstates):
                if np.abs(wvoverlap[istate, istate]) >= 0.5:
                    found_overlap = True
                    if wvoverlap[istate, istate] >= 0.5:
                        phasevector[istate] = +1
                    else:
                        phasevector[istate] = -1

            if found_overlap:
                phases_assigned = True
                tmp_phases[ts] = phasevector

        if line.startswith("! 7 Ekin"):
            e_kin_assigned = True
            tmp_e_kin[ts] = float(next(f).strip())

        if line.startswith("! 8 states (diag, MCH)"):
            pair = next(f).strip().split()
            sdiag_assigned = True
            astate_assigned = True
            tmp_sdiag[ts] = int(pair[0])
            tmp_astate[ts] = int(pair[1])

        if line.startswith("! 15 Gradients (MCH)"):
            state = int(line.strip().split()[-1]) - 1
            force_assigned = True

            for atom in range(natoms):
                tmp_forces[ts, state, atom] = [
                    float(n) for n in next(f).strip().split()
                ]

        if line.startswith("! 16 NACdr matrix element"):
            linecont = line.strip().split()
            si, sj = int(linecont[-2]), int(linecont[-1])

            nacs_assigned = True

            if si < sj:  # elements (si, si) are all zero; elements (sj, si) = -(si, sj)
                sc = idx_table_nacs[(si, sj)]  # statecomb index
                for atom in range(natoms):
                    tmp_nacs[ts, sc, atom, :] = [
                        float(n) for n in next(f).strip().split()
                    ]
            else:  # we can skip the block
                for _ in range(natoms):
                    next(f)

    # has_forces = forces.any(axis=(1, 2, 3))

    if nacs_assigned:
        trajectory_in["nacs"].values = tmp_nacs
        mark_variable_assigned(trajectory_in["nacs"])
    if socs_assigned:
        trajectory_in["socs"].values = tmp_socs
        mark_variable_assigned(trajectory_in["socs"])
        mark_variable_assigned(trajectory_in["full_statecomb"])
        mark_variable_assigned(trajectory_in["full_statecomb_from"])
        mark_variable_assigned(trajectory_in["full_statecomb_to"])
    if force_assigned:
        trajectory_in["forces"].values = tmp_forces
        mark_variable_assigned(trajectory_in["forces"])
    if sdiag_assigned:
        trajectory_in["sdiag"].values = tmp_sdiag
        mark_variable_assigned(trajectory_in["sdiag"])
    if astate_assigned:
        trajectory_in["astate"].values = tmp_astate
        mark_variable_assigned(trajectory_in["astate"])
    if phases_assigned:
        trajectory_in["phases"].values = tmp_phases
        mark_variable_assigned(trajectory_in["phases"])
    if dipole_assigned:
        # post-processing
        # np.diagonal swaps state and direction, so we transpose them back
        trajectory_in.dip_perm.values[:] = np.diagonal(
            tmp_dip_all, axis1=1, axis2=2
        ).transpose(0, 2, 1)
        idxs_dip_trans = (slice(None), *np.triu_indices(nstates, k=1), slice(None))
        trajectory_in.dip_trans.values[:] = tmp_dip_all[idxs_dip_trans]
        mark_variable_assigned(trajectory_in["dip_perm"])
        mark_variable_assigned(trajectory_in["dip_trans"])
    if energy_assigned:
        trajectory_in["energy"].values = tmp_energy
        mark_variable_assigned(trajectory_in["energy"])
    if e_kin_assigned:
        # For now, we do not include e_kin or velocities.
        # mark_variable_assigned(trajectory_in["e_kin"])
        pass

    if not (max_ts + 1 <= nsteps):
        raise ValueError(
            f"Metadata declared {nsteps=} timesteps, but the "
            f"greatest timestep index was {max_ts + 1=}"
        )
    completed = max_ts + 1 == nsteps

    return completed, max_ts + 1, trajectory_in


def parse_trajout_xyz(
    nsteps: int, f: TextIOWrapper
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read atom names, atom numbers and positions for each time step up until a maximum of `nsteps` from an `output.xyz` file and returm them.

    Parameters
    ----------
    nsteps : int
        The maximum number of steps to be read
    f : TextIOWrapper
        The input file wrapper providing the contents of an `output.xyz` file.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Tuple of (atom_names, atom_numbers, atom_positions) as numpy arrays.
        Only atom_positions has the first index indicate the time step, the second the atom and the third the direction.
        Other entries are 1d arrays.
    """
    # TODO: FIXME: inputs in xyz appear to be in angstrom contrary to Bohr in .log and .dat
    first = next(f)
    assert first.startswith(" " * 6)
    natoms = int(first.strip())

    atNames = np.full((natoms), "")
    atNums = np.full((natoms), -1)
    atXYZ = np.full((nsteps, natoms, 3), np.nan)

    ts = 0

    for index, line in enumerate(f):
        if "t=" in line:
            assert ts < nsteps, (
                f"Excess time step at ts={ts}, for a maximum nsteps={nsteps}"
            )
            for atom in range(natoms):
                linecont = re.split(" +", next(f).strip())
                if ts == 0:
                    atNames[atom] = linecont[0]
                atNums[atom] = get_atom_number_from_symbol(linecont[0])
                atXYZ[ts, atom] = [float(n) for n in linecont[1:]]
            ts += 1

    return (atNames, atNums, atXYZ)


def parse_input_settings(input_lines: list[str]) -> dict[str, Any]:
    """Function to parse settings from the `input` file.

    Can be provided the contents of the file as found in the `output.log` file.

    Parameters
    ----------
    input_lines : list[str]
        The lines of the input file.

    Returns
    -------
    dict[str, Any]
        A key-value dictionary, where the keys are the first words in each line
    """
    settings = {}
    for line in input_lines:
        parsed = line.strip().split()
        if len(parsed) == 2:
            settings[parsed[0]] = parsed[1]
        elif len(parsed) > 2:
            settings[parsed[0]] = " ".join(parsed[1:])
        elif len(parsed) == 1:
            settings[parsed[0]] = True
    return settings
