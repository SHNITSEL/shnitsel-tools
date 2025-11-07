from io import TextIOWrapper
import logging
import math
import os
from itertools import combinations
from glob import glob
import pathlib
import re
from typing import Any, Dict, Tuple
import xarray as xr
import pandas as pd
import numpy as np
from pyparsing import nestedExpr
from pprint import pprint

from shnitsel.io.helpers import (
    LoadingParameters,
    PathOptionsType,
    get_atom_number_from_symbol,
    make_uniform_path,
)
from shnitsel.units.definitions import get_default_input_attributes


def parse_pyrai2md(
    traj_path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read a trajector of the PyrAI2md format.

    Args:
        pathlist (PathOptionsType): Path to the directory containing a PyrAI2md output file list
        loading_parameters (LoadingParameters | None, optional): Parameter settings for e.g. standard units or state names.

    Returns:
        xr.Dataset: The Dataset object containing all of the loaded data in default shnitsel units
    """
    # TODO: FIXME: 
    logging.warning("No NACS available for PyrAI2md")

    path_obj: pathlib.Path = make_uniform_path(traj_path)
    # TODO: FIXME: use loading_parameters to configure units and state names
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
        pprint(settings)

    state_ids = np.array(settings["global"]["State order"])
    state_types = np.array(settings["global"]["Multiplicity"])

    nstates = len(state_ids)
    nsinglets = np.sum(state_types == 1)
    ndoublets = np.sum(state_types == 2)
    ntriplets = np.sum(state_types == 3)

    nsteps = settings["md"]["Step"]
    natoms = settings["global"]["Active atoms"]
    delta_t = settings["md"]["Dt (au)"]

    default_attributes = get_default_input_attributes("pyrai2md", loading_parameters)

    trajectory = create_initial_dataset(nsteps, nstates, natoms, loading_parameters)

    trajectory.attrs["delta_t"] = delta_t
    # TODO: FIXME: Make sure this is actually a common naming scheme
    trajectory.attrs["t_max"] = np.nan
    trajectory.attrs["max_ts"] = nsteps
    trajectory.attrs["completed"] = True
    trajectory.attrs["input_format_version"] = settings["version"]

    trajectory.attrs["num_singlets"] = nsinglets
    trajectory.attrs["num_doublets"] = ndoublets
    trajectory.attrs["num_triplets"] = ntriplets

    # TODO: FIXME: apply state names

    with xr.set_options(keep_attrs=True):
        with open(os.path.join(log_paths[0])) as f:
            trajectory, max_ts2 = parse_observables_from_log(f, trajectory)
        trajectory, max_ts1, times = parse_md_energies(md_energies_paths[0], trajectory)

    real_max_ts = min(max_ts1, max_ts2)

    # Cut to actual size
    trajectory = trajectory.isel(time=slice(0, real_max_ts))
    trajectory.assign_coords(time=times)
    trajectory.coords["time"].attrs.update(default_attributes["time"])

    # TODO: FIXME: conflicting dimension sizes "time". We need to deal with trajectory not finishing its full run.
    # One test trajectory did not finish its full number of steps as denoted in the log, so we need to check if the trajectory has finished
    # before sizing the output.

    return trajectory


def create_initial_dataset(
    nsteps: int,
    nstates: int,
    natoms: int,
    loading_parameters: LoadingParameters | None,
    **kwargs,
) -> xr.Dataset:
    """Function to initialize an `xr.Dataset` with appropriate variables and coordinates to acommodate loaded data.

    Args:
        nsteps (int): The number of expected steps in this trajectory
        nstates (int): Number of states within the datasets
        natoms (int): The number of atoms within the datasets

    Returns:
        xr.Dataset: An xarray Dataset with appropriately sized DataArrays and coordinates also including default attributes for all variables.
    """
    template = {
        "energy": ["time", "state"],
        # "dip_all": ["time", "state", "state2", "direction"],
        # "dip_perm": ["time", "state", "direction"],
        # "dip_trans": ["time", "statecomb", "direction"],
        "forces": ["time", "state", "atom", "direction"],
        # 'has_forces': ['placeholder'],
        # 'has_forces': [],
        "phases": ["time", "state"],
        "nacs": ["time", "statecomb", "atom", "direction"],
        "astate": ["time"],
        "atXYZ": ["time", "atom", "direction"],
        "atNames": ["atom"],
        "atNums": ["atom"],
        "state_names": ["state"],
        "state_types": ["state"],
    }

    template_default_values = {
        "energy": np.nan,
        "dip_all": np.nan,
        "dip_perm": np.nan,
        "dip_trans": np.nan,
        "forces": np.nan,
        "phases": np.nan,
        "nacs": np.nan,
        "astate": 0,
        "atXYZ": np.nan,
        "atNames": "",
        "atNums": -1,
        "state_names": "",
        "state_types": 0,
    }

    if natoms == 0:
        # This probably means that check_dims() couldn't find natoms,
        # so we don't expect properties with an atom dimension.
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
        "time": nsteps,
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
    default_sharc_attributes = get_default_input_attributes(
        "pyrai2md", loading_parameters
    )

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
                default_sharc_attributes[varname]
                if varname in default_sharc_attributes
                else {}
            ),
        )
        for varname, dims in template.items()
    }

    res_dataset = xr.Dataset(datavars, coords)

    # Try and set some default attributes on the coordinates for the dataset
    for coord_name in res_dataset.coords:
        if coord_name in default_sharc_attributes:
            res_dataset[coord_name].attrs.update(
                default_sharc_attributes[str(coord_name)]
            )

    res_dataset.attrs["input_format"] = "pyrai2md"
    res_dataset.attrs["input_type"] = "dynamic"

    return res_dataset


def parse_md_energies(
    path: pathlib.Path, trajectory_in: xr.Dataset
) -> Tuple[xr.Dataset, int, np.ndarray]:
    """Function to parse energy and time information into the provided trajectory.

    Returns the number of discovered time steps and the resulting trajectory

    Args:
        path (pathlib.Path): The path to a "*.md.energies" file.
        trajectory_in (xr.Dataset): The trajectory to store the data in

    Returns:
        Tuple[xr.Dataset, int, np.ndarray]: The resulting state of the trajectory after reading and the number of time steps actually found in the data. Finally, the actual values of absolute time read from the energy file for each step
    """
    df = pd.read_csv(path, sep=r"\s+", header=None, skiprows=1).set_index(0)
    df.index.name = "time"
    # We convert later
    # df.index *= 0.5 / 20.67  # convert a.u. to fs
    energy = df.loc[:, 4:].values
    # nstates = len(energy.columns)
    times = df.index.values  # df.loc[:, 0].values

    num_ts = df.shape[0]

    trajectory_in["energy"].values[:num_ts] = energy
    return trajectory_in, num_ts, times
    return (
        xr.Dataset.from_dataframe(energy)
        .to_array("state")
        .assign_coords(state=np.arange(1, nstates + 1))
    )


_int_pattern = re.compile(r"^[-+]?[0-9]+$")
_float_pattern = re.compile(r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$")
_double_whitespace_pattern = re.compile(r"\s{2,}")
_colon_whitespace_pattern = re.compile(r":\s{1,}")
_setting_name_pattern = re.compile(r"(\w+([\s_]{1}\w+)*)")


def read_pyrai2md_settings_from_log(f: TextIOWrapper) -> Dict[str, Any]:
    """Function to read the settings from a pyrai2md log file.

    Args:
        f (TextIOWrapper): The input file stream.

    Returns:
        Dict[str, Any]: The resulting dictionary of settings
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
                    f"Encountered error while parsing pyrai2md settings array: {e}"
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
                        f"No value set for setting {current_section}.{parts[0]}"
                    )
                elif len(parts) >= 2:
                    key = parts[0]
                    if not _setting_name_pattern.match(key):
                        logging.debug(
                            f"Skipping key {key} because it did not conform to usual naming conventions"
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

    logging.info(f"Parsed pyrai2md settings: {settings}")
    return settings


def parse_observables_from_log(
    f: TextIOWrapper, trajectory_in: xr.Dataset
) -> Tuple[xr.Dataset, int]:
    """Function to read multiple observables from a PyrAI2md log file.

    Returns the trajectory with added observables data

    Args:
        f (TextIOWrapper): The file input stream from which the log is read
        trajectory_in (xr.Dataset): The initial trajectory state

    Raises:
        ValueError: If no state information could be read from an iteration's output in a frame
        ValueError: Number of steps could not be read from the trajectory input
        ValueError: Multiple end messages found in log.

    Returns:
        Tuple[xr.Dataset,int]: The updated trajectory and the true final time step seen in the log
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
    veloc = np.full((expected_nsteps, natoms, 3), np.nan)
    dcmat = np.full((expected_nsteps, nstates, nstates), np.nan)

    ts_idx = -1
    end_msg_count = 0
    for line in f:
        # The start of a timestep
        # Iter:        1  Ekin =           0.1291084223229551 au T =   300.00 K dt =         20 CI:   3
        # Root chosen for geometry opt   2
        if line.startswith("  Iter:"):
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
                raise ValueError(f"No state info found for Iter: {ts_idx+1}")

        # Positions:
        #   &coordinates in Angstrom
        # -------------------------------------------------------------------------------
        # C          0.5765950000000000     -0.8169010000000000     -0.0775610000000000
        # C          1.7325100000000000     -0.1032670000000000      0.1707480000000000
        # -------------------------------------------------------------------------------
        if line.startswith("  &coordinates"):
            hline = next(f)
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
        if line.startswith("  &velocities"):
            hline = next(f)
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
        if line.startswith("  &gradient"):
            istate = int(line.strip().split()[2]) - 1
            hline = next(f)
            assert hline.startswith("---")
            for iatom in range(natoms):
                forces[ts_idx, istate, iatom] = np.asarray(
                    next(f).strip().split()[1:], dtype=float
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
        if line.startswith("  &derivative coupling matrix"):
            hline = next(f)
            assert hline.startswith("---")
            for istate1 in range(nstates):
                dcmat[ts_idx, istate1] = np.asarray(
                    next(f).strip().split(), dtype=float
                )
            hline = next(f)
            assert hline.startswith("---")

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
        if line.startswith("  &surface hopping information"):
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

    trajectory_in["forces"].values = forces
    trajectory_in["atXYZ"].values = atXYZ
    trajectory_in.attrs["completed"] = (
        end_msg_count == 1 or real_max_ts >= expected_nsteps
    )

    # TODO: FIXME: Do we need Phases to be included?

    return (
        trajectory_in.assign_coords(
            {
                "atNames": atNames,
                "atNums": [get_atom_number_from_symbol(x) for x in atNames],
            }
        ),
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
