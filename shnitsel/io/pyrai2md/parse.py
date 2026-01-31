from io import TextIOWrapper
import logging
import os
import pathlib
import re
from typing import Any, Dict, Tuple
import xarray as xr
import pandas as pd
import numpy as np
from pyparsing import nestedExpr
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

from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
from shnitsel.data.atom_helpers import get_atom_number_from_symbol

_int_pattern = re.compile(r"^[-+]?[0-9]+$")
_float_pattern = re.compile(r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$")
_double_whitespace_pattern = re.compile(r"\s{2,}")
_colon_whitespace_pattern = re.compile(r":\s{1,}")
_setting_name_pattern = re.compile(r"(\w+([\s_]{1}\w+)*)")


_re_nac_header_line = re.compile(
    r"&nonadiabatic coupling vectors (?P<state_1>\d+) *- *(?P<state_2>\d+) in Hartree/Bohr +M *= *(?P<mult_1>\d+) */ *(?P<mult_2>\d+)"
)
_re_soc_line = re.compile(
    r"((?P<compute_missing>Not computed)|<H>=(?P<coupling>(\d\.)+)) +(?P<state_1>\d+) *- *(?P<state_2>\d+) *in cm-1 M1 *= *(?P<mult_1>\d+) *M2 *= *(?P<mult_2>\d+)"
)


def parse_pyrai2md(
    traj_path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read a trajector of the PyrAI2md format.

    Parameters
    ----------
    pathlist : PathOptionsType
        Path to the directory containing a PyrAI2md output file list
    loading_parameters : LoadingParameters | None, optional
        Parameter settings for e.g. standard units or state names.

    Returns
    -------
    xr.Dataset
        The Dataset object containing all of the loaded data in default shnitsel units
    """

    # TODO: FIXME: Pyrai2md has a different understanding of the number of triplets: It states per state whether it is singlet, doublet, triplet, no further multiplicities are considered.

    path_obj: pathlib.Path = make_uniform_path(traj_path)
    # TODO: FIXME: Check if there are other files of pyrai2md trajectories to read information from.
    md_energies_paths = list(path_obj.glob("*.md.energies"))
    if (n := len(md_energies_paths)) != 1:
        raise FileNotFoundError(
            "Expected to find a single file ending with '.md.energies' "
            f"but found {n} files: {md_energies_paths}"
        )
    log_paths = list(path_obj.glob("*.log"))
    if (n := len(log_paths)) != 1:
        raise FileNotFoundError(
            "Expected to find a single file ending with '.log' "
            f"but found {n} files: {log_paths}"
        )

    with open(log_paths[0]) as f:
        settings = read_pyrai2md_settings_from_log(f)
        # pprint(settings)

    state_ids = np.array(settings["global"]["State order"])
    state_types = np.array(settings["global"]["Multiplicity"])

    nstates = len(state_ids)
    # TODO: FIXME: This sounds somewhat wrong in relation to how we set the state names. There are not enough states for ntriplets*3 T states!
    nsinglets = np.sum(state_types == 1)
    ndoublets = np.sum(state_types == 2)
    ntriplets = np.sum(state_types == 3)

    nsteps = settings["md"]["Step"]
    natoms = settings["global"]["Active atoms"]
    delta_t = settings["md"]["Dt (au)"]

    trajectory, default_format_attributes = create_initial_dataset(
        nsteps, nstates, natoms, "pyrai2md", loading_parameters
    )

    # TODO: FIXME: apply state names

    with xr.set_options(keep_attrs=True):
        trajectory, max_ts1, times = parse_md_energies(md_energies_paths[0], trajectory)

        with open(os.path.join(log_paths[0])) as f:
            trajectory, max_ts2 = parse_observables_from_log(f, trajectory)

    completed = (max_ts1 == nsteps) and (max_ts2 == nsteps)
    real_max_ts = min(max_ts1, max_ts2)

    # Cut to actual size
    trajectory = trajectory.isel(time=slice(0, real_max_ts))

    trajectory = trajectory.assign_coords(
        {
            "time": ("time", times, default_format_attributes["time"]),
            "state_types": (
                "state",
                state_types,
                default_format_attributes["state_types"],
            ),
        }
    )
    mark_variable_assigned(trajectory["time"])
    mark_variable_assigned(trajectory["state_types"])
    # print(repr(trajectory))
    # print(repr(trajectory.energy))
    # print(repr(trajectory.forces))
    # print(repr(trajectory.atXYZ))

    # Set all settings we require to be present on the trajectory
    required_settings = RequiredTrajectorySettings(
        nsteps * delta_t,
        delta_t,
        real_max_ts,
        completed,
        "pyrai2md",
        "dynamic",
        settings["version"],
        nsinglets,
        ndoublets,
        ntriplets,
    )
    trajectory = assign_required_settings(trajectory, required_settings)

    optional_settings = OptionalTrajectorySettings(
        has_forces=is_variable_assigned(trajectory["forces"]),
        misc_input_settings=settings,
    )
    trajectory= assign_optional_settings(trajectory, optional_settings)

    return trajectory


def parse_md_energies(
    path: pathlib.Path, trajectory_in: xr.Dataset
) -> Tuple[xr.Dataset, int, np.ndarray]:
    """Function to parse energy and time information into the provided trajectory.

    Returns the number of discovered time steps and the resulting trajectory

    Parameters
    ----------
    path : pathlib.Path
        The path to a "*.md.energies" file.
    trajectory_in : xr.Dataset
        The trajectory to store the data in

    Returns
    -------
    Tuple[xr.Dataset, int, np.ndarray]
        The resulting state of the trajectory after reading and the number of time steps actually found in the data. Finally, the actual values of absolute time read from the energy file for each step
    """
    df = pd.read_csv(path, sep=r"\s+", header=None, skiprows=1).set_index(0)
    df.index.name = "time"
    # We convert later
    # df.index *= 0.5 / 20.67  # convert a.u. to fs
    energy = df.loc[:, 4:].values
    # nstates = len(energy.columns)
    times = df.index.values  # df.loc[:, 0].values
    # print("E:", energy)
    # print("t:", times)

    num_ts = df.shape[0]
    e_kin = df.loc[:, 2].values

    trajectory_in["energy"].values[:num_ts] = energy
    mark_variable_assigned(trajectory_in["energy"])
    trajectory_in["e_kin"].values[:num_ts] = e_kin
    mark_variable_assigned(trajectory_in["e_kin"])
    return trajectory_in, num_ts, times


def read_pyrai2md_settings_from_log(f: TextIOWrapper) -> dict[str, Any]:
    """Function to read the settings from a pyrai2md log file.

    Parameters
    ----------
    f : TextIOWrapper
        The input file stream.

    Returns
    -------
    dict[str, Any]
        The resulting dictionary of settings
    """
    settings: Dict[str, Any] = {}
    settings["global"] = {}

    def decode_setting_string(value_string: str) -> Any:
        # Catch integers and floats first
        if _int_pattern.match(value_string):
            return int(value_string)
        elif _float_pattern.match(value_string):
            return float(value_string)
        elif value_string.startswith("'"):
            # We have a delimited string
            res_string = value_string[
                value_string.find("'") + 1 : value_string.rfind("'")
            ]
            return str(res_string)
        elif value_string.find("[") == -1:
            # We have a string if there isn't an array happening
            return str(value_string)
        else:
            # We have an array to decode:
            try:
                relevant_part = value_string[
                    value_string.find("[") : value_string.rfind("]") + 1
                ]
                decoded_array = nestedExpr("[", "]").parseString(relevant_part).asList()

                # Cut off surrounding array
                if len(decoded_array) > 0:
                    decoded_array = decoded_array[0]

                def recurse_parse(data):
                    if isinstance(data, list):
                        return [
                            converted
                            for x in data
                            if not x == ","
                            and (converted := recurse_parse(x)) is not None
                        ]
                    else:
                        return decode_setting_string(data.strip(","))

                # logging.debug(f"Settings array after decomposition: {decoded_array}")
                parsed_array = recurse_parse(decoded_array)
                # logging.debug(f"Settings array after parsing: {parsed_array}")
                return parsed_array

            except Exception as e:
                logging.debug(
                    "Encountered error while parsing pyrai2md settings array: %(e)s",
                    {'e': e},
                )
                return None

    has_had_section = False
    current_section = None
    current_section_settings = None

    while curr_line := next(f):
        curr_stripped = curr_line.strip()
        if curr_stripped.startswith("Iter:"):
            # We reached the beginning of actual frames
            break
        if len(curr_stripped) == 0:
            # Skip empty lines
            continue

        if curr_stripped.startswith("version:"):
            # Read version string
            version_string = curr_stripped[len("version:") :].strip()
            settings["version"] = version_string
        elif curr_stripped.startswith("&"):
            # Beginning of a section Header
            section_name = curr_stripped[1:].strip()
            current_section = section_name
            has_had_section = True
        elif curr_stripped.startswith("---"):
            # Beginning or end of a section body
            if current_section is not None:
                if current_section_settings is None:
                    # Start section parameters
                    current_section_settings = {}
                else:
                    # End section parameters
                    if current_section is not None:
                        settings[current_section] = current_section_settings
                        # Reset section
                        current_section = None
                        current_section_settings = None
        else:
            # We are in a section block
            if len(curr_stripped) > 0:
                kv_split = _double_whitespace_pattern.split(curr_stripped)
                if len(kv_split) < 2:
                    kv_split = _colon_whitespace_pattern.split(curr_stripped)
                parts = [x.strip().strip(":") for x in kv_split]

                if len(parts) == 1:
                    logging.debug(
                        "No value set for setting %(section)s.%(setting)s",
                        {'section': current_section, 'setting': parts[0]},
                    )
                elif len(parts) >= 2:
                    key = parts[0]
                    if not _setting_name_pattern.match(key):
                        logging.debug(
                            "Skipping key %(key)s because it did not conform to usual naming conventions",
                            {'key': key},
                        )
                        continue

                    value = "  ".join(parts[1:])
                    if current_section_settings is not None:
                        current_section_settings[key] = decode_setting_string(value)
                    elif has_had_section:
                        # Global settings have different arrays...
                        if len(parts) > 2:
                            settings["global"][key] = [
                                decode_setting_string(x) for x in parts[1:]
                            ]
                        else:
                            settings["global"][key] = decode_setting_string(parts[1])

    logging.debug("Parsed pyrai2md settings: %(settings)s", {'settings': settings})
    return settings


def parse_observables_from_log(
    f: TextIOWrapper, trajectory_in: xr.Dataset
) -> Tuple[xr.Dataset, int]:
    """Function to read multiple observables from a PyrAI2md log file.

    Returns the trajectory with added observables data.
    To calculate the NACs accurately, the trajectory requires state-energy data to be set so that we can
    renormalize the NACs data with energy deltas.

    Parameters
    ----------
    f : TextIOWrapper
        The file input stream from which the log is read
    trajectory_in : xr.Dataset
        The initial trajectory state

    Raises
    ------
    ValueError
        If no state information could be read from an iteration's output in a frame
    ValueError
        Number of steps could not be read from the trajectory input
    ValueError
        Multiple end messages found in log.

    Returns
    -------
    Tuple[xr.Dataset,int]
        The updated trajectory and the true final time step seen in the log
    """
    # Read MD settings:
    #  &md
    # -------------------------------------------------------
    #  Initial state:              2
    #  Initialize random velocity  0
    #  Temperature (K):            300
    #  Step:                       4000
    #  Dt (au):                    20.67

    expected_nsteps: int | None = trajectory_in.sizes["time"]
    nstates: int = trajectory_in.sizes["state"]
    natoms: int = trajectory_in.sizes["atom"]
    nstatecomb: int = trajectory_in.sizes["statecomb"]
    nfull_statecomb: int = trajectory_in.sizes["full_statecomb"]

    if expected_nsteps is None:
        raise ValueError("Could not read `nsteps` from trajectory")

    # Read final settings:
    # *---------------------------------------------------*
    # |                                                   |
    # |          Nonadiabatic Molecular Dynamics          |
    # |                                                   |
    # *---------------------------------------------------*
    #
    #
    # State order:         1   2   3
    # Multiplicity:        1   1   1
    #
    # QMMM key:         None
    # QMMM xyz          Input
    # Active atoms:     45
    # Inactive atoms:   0
    # Link atoms:       0
    # Highlevel atoms:  45
    # Midlevel atoms:   0
    # Lowlevel atoms:   0

    while not next(f).startswith(" *---"):
        pass

    # Set up numpy arrays
    explicit_ts = np.full((expected_nsteps,), -1, dtype=int)
    astate = np.full((expected_nsteps), -1, dtype=int)
    forces = np.full((expected_nsteps, nstates, natoms, 3), np.nan)
    atXYZ = np.full((expected_nsteps, natoms, 3), np.nan)
    atNames = np.full((natoms), "", dtype=str)
    got_atNames = False
    # TODO: Use velocities?
    veloc = np.full((expected_nsteps, natoms, 3), np.nan)
    dcmat = np.full((expected_nsteps, nstates, nstates), np.nan)
    nacs = np.full((expected_nsteps, nstates, nstates), np.nan)
    nacs = np.full((expected_nsteps, nstatecomb, natoms, 3), np.nan)
    socs = np.full((expected_nsteps, nfull_statecomb), np.nan + 0j, dtype=np.complex128)

    has_forces = False
    has_veloc = False
    has_dcms = False
    has_nacs = False
    has_socs = False
    has_positions = False

    has_energy_deltas = is_variable_assigned(trajectory_in.energy)
    energies = trajectory_in.energy.values

    energy_deltas = {
        (i + 1, j + 1): energies[:, i] - energies[:, j]
        for i in range(nstates)
        for j in range(nstates)
    }

    state_comb_order = trajectory_in.statecomb.values
    state_comb_dict: dict[tuple[int, int], int] = {
        pair: i for i, pair in enumerate(state_comb_order)
    }

    idx_table_socs = {
        (si, sj): idx
        for idx, (si, sj) in enumerate(trajectory_in.full_statecomb.values)
    }
    # logging.debug(state_comb_order)
    # logging.debug(state_comb_dict)
    # logging.debug(idx_table_socs)

    # TODO: FIXME: Read variable units from the file and compare to expected values or override.
    ts_idx = -1
    end_msg_count = 0
    for line in f:
        line = line.strip()
        # The start of a timestep
        # Iter:        1  Ekin =           0.1291084223229551 au T =   300.00 K dt =         20 CI:   3
        # Root chosen for geometry opt   2
        if line.startswith("Iter:"):
            ts_idx += 1
            explicit_ts[ts_idx] = int(line.strip().split()[1])

            for _ in range(10):
                # Get active state
                # A surface hopping is not allowed
                # **
                # At state:   2
                line = next(f).strip()
                if line.startswith("At state"):
                    astate[ts_idx] = int(line.strip().split()[2])
                    break
                # A surface hopping event happened
                # **
                # From state:   2 to state:   3 *
                elif line.startswith("From state"):
                    astate[ts_idx] = int(line.strip().split()[5])
                    break
            else:
                raise ValueError(f"No state info found for Iter: {ts_idx + 1}")

        # Positions:
        #   &coordinates in Angstrom
        # -------------------------------------------------------------------------------
        # C          0.5765950000000000     -0.8169010000000000     -0.0775610000000000
        # C          1.7325100000000000     -0.1032670000000000      0.1707480000000000
        # -------------------------------------------------------------------------------
        if line.startswith("&coordinates"):
            hline = next(f)
            has_positions = True

            assert hline.startswith("---")

            for iatom in range(natoms):
                content = next(f).strip().split()
                atXYZ[ts_idx, iatom] = np.asarray(content[1:], dtype=float)
                if not got_atNames:
                    atNames[iatom] = str(content[0])

            got_atNames = True

            hline = next(f)
            assert hline.startswith("---")

        # Velocities:
        #   &velocities in Bohr/au
        # -------------------------------------------------------------------------------
        # C          0.0003442000000000      0.0001534200000000     -0.0000597200000000
        # C         -0.0005580000000000      0.0003118300000000     -0.0000154900000000
        # -------------------------------------------------------------------------------
        if line.startswith("&velocities"):
            hline = next(f)
            has_veloc = True
            assert hline.startswith("---")
            for iatom in range(natoms):
                veloc[ts_idx, iatom] = np.asarray(
                    next(f).strip().split()[1:], dtype=float
                )
            hline = next(f)
            assert hline.startswith("---")

        # Forces:
        #   &gradient state               1 in Eh/Bohr
        # -------------------------------------------------------------------------------
        # C         -0.0330978534152795      0.0073099255379017      0.0082666356536386
        # C          0.0313629524413876      0.0196036465968827      0.0060952442704520
        # -------------------------------------------------------------------------------
        if line.startswith("&gradient"):
            istate = int(line.strip().split()[2]) - 1
            hline = next(f)
            has_forces = True
            assert hline.startswith("---")
            for iatom in range(natoms):
                next_line = next(f).strip()
                if next_line.lower().find("not computed") > 0:
                    # Skip if no computation has happened
                    break
                forces[ts_idx, istate, iatom] = np.asarray(
                    next_line.split()[1:], dtype=float
                )
            hline = next(f)
            assert hline.startswith("---")

        # Derivative coupling matrix:
        #  &derivative coupling matrix
        # -------------------------------------------------------------------------------
        #       0.0000000000000000       0.0000000000000004      -0.0000000000000001
        #      -0.0000000000000004       0.0000000000000000       0.0000000000000003
        #       0.0000000000000001      -0.0000000000000003       0.0000000000000000
        # -------------------------------------------------------------------------------

        # '&derivative coupling matrix
        # -------------------------------------------------------------------------------
        # %s-------------------------------------------------------------------------------
        # '% (print_matrix(self.traj.nac))
        if line.startswith("&derivative coupling matrix"):
            hline = next(f)
            has_dcms = True
            assert hline.startswith("---")
            for istate1 in range(nstates):
                dcmat[ts_idx, istate1] = np.asarray(
                    next(f).strip().split(), dtype=float
                )
            hline = next(f)
            assert hline.startswith("---")

        # ' &nonadiabatic coupling vectors %3d - %3d in Hartree/Bohr M = %1d / %1d
        # -------------------------------------------------------------------------------
        # %s-------------------------------------------------------------------------------
        # ' % (s1 + 1, s2 + 1, m1, m2, print_coord(np.concatenate((self.traj.atoms, coupling), axis=1)))
        #
        # or
        # ' &nonadiabatic coupling vectors %3d - %3d in Hartree/Bohr M = %1d / %1d
        # -------------------------------------------------------------------------------
        # Not computed
        # -------------------------------------------------------------------------------
        # '
        if line.startswith("&nonadiabatic coupling vectors"):
            has_nacs = True
            match = _re_nac_header_line.match(line)
            if match:
                from_state = int(match.group("state_1"))
                to_state = int(match.group("state_2"))
                from_state_mult = int(match.group("mult_1"))
                to_state_mult = int(match.group("mult_2"))
                comb_index = state_comb_dict[(from_state, to_state)]

                if not has_energy_deltas:
                    raise ValueError(
                        "To normalize PyRAI2md NACs, we need to have energy delta values. Please make sure the energy delta values can be read before the NAC entries in the `.log` file."
                    )

                # TODO: FIXME: Check if this needs the opposite sign
                delta_e = energy_deltas[(from_state, to_state)]
                nextline = next(f).strip()
                assert nextline.startswith("---")
                for iatom in range(natoms):
                    nextline = next(f).strip()
                    # We normalize with the \Delta E_ij to get 1/Bohr units
                    vec = [float(x) / delta_e for x in nextline.split()]
                    nacs[ts_idx, comb_index, iatom, :] = vec

                nextline = next(f).strip()
                assert nextline.startswith("---")

            else:
                logging.warning(
                    "Malformed NAC line in PyRAI2md trajectory: %(line)s",
                    {'line': line},
                )

        # soc_info = ''
        #         for n, pair in enumerate(self.traj.soc_coupling):
        #             s1, s2 = pair
        #             m1 = self.traj.statemult[s1]
        #             m2 = self.traj.statemult[s2]
        #             try:
        #                 coupling = self.traj.soc[n]
        #                 soc_info += '  <H>=%10.4f            %3d - %3d in cm-1 M1 = %1d M2 = %1d\n' % (
        #                     coupling, s1 + 1, s2 + 1, m1, m2)

        #             except IndexError:
        #                 soc_info += '  Not computed              %3d - %3d in cm-1 M1 = %1d M2 = %1d\n' % (
        #                     s1 + 1, s2 + 1, m1, m2)

        #         if len(self.traj.soc_coupling) > 0:
        # log_info += '
        # &spin-orbit coupling
        # -------------------------------------------------------------------------------
        # %s-------------------------------------------------------------------------------
        # ' % soc_info
        if line.startswith("&spin-orbit coupling"):
            has_socs = True
            nextline = next(f).strip()
            assert nextline.startswith("---")
            nextline = next(f).strip()
            while not nextline.startswith("---"):
                match = _re_soc_line.match(nextline)

                if match:
                    from_state = int(match.group("state_1"))
                    to_state = int(match.group("state_2"))
                    from_state_mult = int(match.group("mult_1"))
                    to_state_mult = int(match.group("mult_2"))

                    if "compute_missing" in match.groupdict():
                        logging.info(
                            "Soc missing for %(from_state)d -> %(to_state)d (mults: %(from_state_mult)d -> %(to_state_mult)d)",
                            {
                                'from_state': from_state,
                                'to_state': to_state,
                                'from_state_mult': from_state_mult,
                                'to_state_mult': to_state_mult,
                            },
                        )
                    else:
                        socs[ts_idx, idx_table_socs[(from_state, to_state)]] = (
                            float(match.group("coupling")) + 0j
                        )

                nextline = next(f).strip()

        # Surface hopping information at the end of each timestep:
        #  &surface hopping information
        # -------------------------------------------------------
        #
        #     Random number:             0.15725129
        #     Accumulated probability:   0.00000000
        #     state mult  level   probability
        #     1     1     1       0.00000000
        #     2     1     2       0.00000000
        #     3     1     3       0.00000000
        #
        #
        # -------------------------------------------------------
        if line.startswith("&surface hopping information"):
            hline = next(f)
            assert hline.startswith("---")
            # We don't currently parse this
            while not next(f).startswith("---"):
                pass

        # Completion indicator:
        # Nonadiabatic Molecular Dynamics End:  2025-04-13 01:12:26 Total:     0 days    15 hours    59 minutes    20 seconds
        if line.startswith("Nonadiabatic Molecular Dynamics End:"):
            end_msg_count += 1

    if end_msg_count > 1:
        raise ValueError(
            'Completion message "Nonadiabatic Molecular Dynamics End:" appeared '
            f"{end_msg_count} times"
        )

    real_max_ts = explicit_ts.max()
    trajectory_in["astate"].values = astate
    mark_variable_assigned(trajectory_in["astate"])
    if has_veloc:
        trajectory_in["velocities"].values = veloc
        mark_variable_assigned(trajectory_in["velocities"])

    if has_dcms:
        logging.warning(
            "DCM currently not processed on PyRAI2md trajectories due to shape mismatch"
        )
        # TODO: FIXME: Deal with dcm shape and adding dcm to default template
        # trajectory_in["dcm"].values = dcmat
        # mark_variable_assigned(trajectory_in["dcm"])

    if has_nacs:
        trajectory_in["nacs"].values = nacs
        mark_variable_assigned(trajectory_in["nacs"])
    else:
        # NOTE: Use /PyrAI2md/Dynamics/aimd.py to complete the reading of NACS/dcm
        logging.info("No NACS available for PyrAI2md file")

    if has_socs:
        logging.warning(
            "SOCs from PyrAI2md files have not been tested. There may be a mismatch in dimensionality."
        )
        # TODO: FIXME: Deal with soc shape
        trajectory_in["socs"].values = socs
        mark_variable_assigned(trajectory_in["socs"])
        mark_variable_assigned(trajectory_in["full_statecomb"])
        mark_variable_assigned(trajectory_in["full_statecomb_from"])
        mark_variable_assigned(trajectory_in["full_statecomb_to"])

    if has_forces:
        trajectory_in["forces"].values = forces
        mark_variable_assigned(trajectory_in["forces"])
    if has_positions:
        trajectory_in["atXYZ"].values = atXYZ
        mark_variable_assigned(trajectory_in["atXYZ"])
    trajectory_in.attrs["completed"] = (
        end_msg_count == 1 or real_max_ts >= expected_nsteps
    )

    # TODO: FIXME: Do we need Phases to be included?
    trajectory_in = trajectory_in.assign_coords(
        {
            "atNames": ("atom", atNames, trajectory_in.atNames.attrs),
            "atNums": (
                "atom",
                [get_atom_number_from_symbol(x) for x in atNames],
                trajectory_in.atNums.attrs,
            ),
        }
    )
    mark_variable_assigned(trajectory_in["atNames"])
    mark_variable_assigned(trajectory_in["atNums"])
    return (
        trajectory_in,
        real_max_ts,
    )

    """return xr.Dataset(
        {
            # 'dip_all': (['ts', 'state', 'state2', 'direction'], dip_all),
            # 'dip_perm': (['ts', 'state', 'direction'], dip_perm),
            # 'dip_trans': (['ts', 'statecomb', 'direction'], dip_trans),
            # 'sdiag': (['ts'], sdiag),
            "astate": (["ts"], astate, {"long_name": "active state"}),
            "forces": (
                ["ts", "state", "atom", "direction"],
                forces,
                {"units": "hartree/bohr", "unitdim": "Force"},
            ),
            # 'has_forces': (['ts'], has_forces),
            # 'phases': (['ts', 'state'], phases),
            # 'nacs': (
            #     ['ts', 'statecomb', 'atom', 'direction'],
            #     nacs,
            #     {'long_name': "nonadiabatic couplings", 'units': "au"},
            # ),
            "atXYZ": (["ts", "atom", "direction"], atXYZ),
            "dcmat": (["ts", "state", "state2"], dcmat),
        },
        coords=coords,
        attrs={
            "max_ts": ,
            "completed": end_msg_count == 1,
        },
    )"""
