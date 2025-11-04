from io import TextIOWrapper
from typing import Any, Dict, Tuple
import numpy as np
import xarray as xr
from itertools import combinations
import pandas as pd
import logging
import os
import re
import math

from shnitsel.io.helpers import PathOptionsType, get_atom_number_from_symbol
from shnitsel.units.definitions import get_default_input_attributes
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults


from shnitsel.io.helpers import LoadingParameters


def read_traj(
    traj_path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read a single SHARC trajectory directory

    Args:
        traj_path (PathOptionsType): The path to load the trajectory form
        loading_parameters (LoadingParameters | None, optional): Parameter settings for e.g. standard units or state names.

    Returns:
        xr.Dataset: The parsed SHARC directory as a Dataset
    """
    # TODO: FIXME: use loading_parameters to configure units and state names

    # Read some settings from input
    # In particular, if a trajectory is extended by increasing
    # tmax and resuming, the header of output.dat will give
    # only the original nsteps, leading to an ndarray IndexError
    input_path = os.path.join(traj_path, "input")
    if os.path.isfile(input_path):
        with open(input_path) as f:
            settings = parse_input_settings(f)
        delta_t = float(settings["stepsize"])
        tmax = float(settings["tmax"])
        nsteps = int(tmax / delta_t) + 1
    else:
        delta_t = None
        tmax = None
        nsteps = None

    with open(os.path.join(traj_path, "output.dat")) as f:
        single_traj = parse_trajout_dat(
            f, nsteps=nsteps, loading_parameters=loading_parameters
        )

    # TODO: Note that for consistency, we renamed the ts dimension to time to agree with other format

    nsteps = single_traj.sizes["time"]

    with open(os.path.join(traj_path, "output.xyz")) as f:
        atNames, atNums, atXYZ = parse_trajout_xyz(nsteps, f)

    single_traj.coords["atNames"] = "atom", atNames
    single_traj.coords["atNums"] = "atom", atNums

    single_traj["atXYZ"][...] = atXYZ

    if delta_t is not None:
        single_traj.attrs["delta_t"] = delta_t

    if tmax is not None:
        single_traj.attrs["t_max"] = tmax

    single_traj.attrs["input_format"] = "sharc"
    single_traj.attrs["input_type"] = "dynamic"

    return convert_all_units_to_shnitsel_defaults(single_traj)


def parse_trajout_dat(
    f: TextIOWrapper,
    nsteps: int | None = None,
    loading_parameters: LoadingParameters = None,
) -> xr.Dataset:
    """Function to parse the contents of an 'output.dat' in a sharc trajectory output directory into a Dataset.

    Args:
        f (TextIOWrapper): A file wrapper providing the contents of 'output.dat'.
        nsteps (int | None, optional): The number of maximum steps expected. Defaults to None.

    Raises:
        ValueError: Raised if not enough steps are found in the output.dat file

    Returns:
        xr.Dataset: The full dataset with unit attributes and further helpful attributes applied.
    """
    settings = {}
    for line in f:
        if line.startswith("*"):
            break

        parsed = line.strip().split()
        if len(parsed) == 2:
            settings[parsed[0]] = parsed[1]
        elif len(parsed) > 2:
            settings[parsed[0]] = parsed[1:]
        else:
            logging.warning("Key without value in settings of output.dat")

    nsteps_output_dat = int(settings["nsteps"]) + 1  # let's not forget ts=0
    if nsteps is None or nsteps < nsteps_output_dat:
        nsteps = nsteps_output_dat
        logging.debug(f"nsteps = {nsteps}")
    else:
        logging.debug(f"(From input file) nsteps = {nsteps}")
    natoms = int(settings["natom"])  # yes, really 'natom', not 'natoms'!
    logging.debug(f"natoms = {natoms}")
    ezero = float(settings["ezero"])
    logging.debug(f"ezero = {ezero}")
    state_settings = [int(s) for s in settings["nstates_m"]]
    state_settings += [0] * (3 - len(state_settings))
    nsinglets, ndoublets, ntriplets = state_settings
    nstates = nsinglets + 2 * ndoublets + 3 * ntriplets
    logging.debug(f"nstates = {nstates}")

    idx_table_nacs = {
        (si, sj): idx
        for idx, (si, sj) in enumerate(combinations(range(1, nstates + 1), 2))
    }

    template = {
        "energy": ["time", "state"],
        "e_kin": ["time"],
        "dip_all": ["time", "state", "state2", "direction"],
        "dip_perm": ["time", "state", "direction"],
        "dip_trans": ["time", "statecomb", "direction"],
        "forces": ["time", "state", "atom", "direction"],
        # 'has_forces': ['placeholder'],
        # 'has_forces': [],
        "phases": ["time", "state"],
        "nacs": ["time", "statecomb", "atom", "direction"],
        "atXYZ": ["time", "atom", "direction"],
        "atNames": ["atom"],
        "atNums": ["atom"],
        "state_names": ["state"],
        "state_types": ["state"],
        "astate": ["time"],
        "sdiag": ["time"],
    }
    dim_lengths = {
        "time": nsteps,
        "state": nstates,
        "state2": nstates,
        "atom": natoms,
        "direction": 3,
        "statecomb": math.comb(nstates, 2),
    }

    template_default_values = {
        "energy": np.nan,
        "e_kin": np.nan,
        "dip_all": np.nan,
        "dip_perm": np.nan,
        "dip_trans": np.nan,
        "sdiag": -1,
        "astate": -1,
        "forces": np.nan,
        "phases": np.nan,
        "nacs": np.nan,
        "atXYZ": np.nan,
        "state_names": "",
        "atNames": "",
        "atNums": -1,
        "state_types": 0,
    }

    def default_fill_prop(name):
        return default_fill(template[name], template_default_values[name])

    def default_fill(dims, default_value):
        return np.full([dim_lengths[d] for d in dims], default_value)

    # now we know the number of steps, we can initialize the data arrays:
    positions = default_fill_prop("atXYZ")
    energy = default_fill_prop("energy")
    e_kin = default_fill_prop("e_kin")
    dip_all = default_fill_prop("dip_all")
    phases = default_fill_prop("phases")
    astate = default_fill_prop("astate")
    sdiag = default_fill_prop("sdiag")
    forces = default_fill_prop("forces")
    nacs = default_fill_prop("nacs")

    atNums = default_fill_prop("atNums")
    atNames = default_fill_prop("atNames")
    state_names = default_fill_prop("state_names")

    state_types = default_fill_prop("state_types")

    state_types[:nsinglets] = 1
    state_names[:nsinglets] = [f"S{i}" for i in range(nsinglets)]
    state_types[nsinglets : nsinglets + 2 * ndoublets] = 2
    state_names[nsinglets : nsinglets + 2 * ndoublets] = [
        f"D{i}" for i in range(2 * ndoublets)
    ]
    state_types[nsinglets + 2 * ndoublets :] = 3
    state_names[nsinglets + 2 * ndoublets :] = [f"T{i}" for i in range(3 * ntriplets)]

    max_ts = -1

    # skip through until initial step:
    for line in f:
        if line.startswith("! 0 Step"):
            ts = int(next(f).strip())
            if ts != 0:
                logging.warning("Initial timestep's index is not 0")
            max_ts = max(max_ts, ts)
            break

    for index, line in enumerate(f):
        if line[0] != "!":
            continue

        if line.startswith("! 0 Step"):
            # update `ts` to current timestep #
            new_ts = int(next(f).strip())
            if new_ts != (ts or 0) + 1:
                logging.warning(f"Non-consecutive timesteps: {ts} -> {new_ts}")
            ts = new_ts
            max_ts = max(max_ts, ts)
            logging.debug(f"timestep = {ts}")

        if line.startswith("! 1 Hamiltonian"):
            for istate in range(nstates):
                energy[ts, istate] = float(next(f).strip().split()[istate * 2]) + ezero

        if line.startswith("! 3 Dipole moments"):
            direction = {"X": 0, "Y": 1, "Z": 2}[line.strip().split()[4]]
            for istate in range(nstates):
                linecont = next(f).strip().split()
                # delete every second element in list (imaginary values, all zero)
                dip_all[ts, istate, :, direction] = [float(i) for i in linecont[::2]]

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
                phases[ts] = phasevector

        if line.startswith("! 7 Ekin"):
            e_kin[ts] = float(next(f).strip())

        if line.startswith("! 8 states (diag, MCH)"):
            pair = next(f).strip().split()
            sdiag[ts] = int(pair[0])
            astate[ts] = int(pair[1])

        if line.startswith("! 15 Gradients (MCH)"):
            state = int(line.strip().split()[-1]) - 1

            for atom in range(natoms):
                forces[ts, state, atom] = [float(n) for n in next(f).strip().split()]

        if line.startswith("! 16 NACdr matrix element"):
            linecont = line.strip().split()
            si, sj = int(linecont[-2]), int(linecont[-1])

            if si < sj:  # elements (si, si) are all zero; elements (sj, si) = -(si, sj)
                sc = idx_table_nacs[(si, sj)]  # statecomb index
                for atom in range(natoms):
                    nacs[ts, sc, atom, :] = [float(n) for n in next(f).strip().split()]
            else:  # we can skip the block
                for _ in range(natoms):
                    next(f)

    # post-processing
    # np.diagonal swaps state and direction, so we transpose them back
    dip_perm = np.diagonal(dip_all, axis1=1, axis2=2).transpose(0, 2, 1)
    idxs_dip_trans = (slice(None), *np.triu_indices(nstates, k=1), slice(None))
    dip_trans = dip_all[idxs_dip_trans]
    # has_forces = forces.any(axis=(1, 2, 3))

    if not (max_ts + 1 <= nsteps):
        raise ValueError(
            f"Metadata declared {nsteps=} timesteps, but the "
            f"greatest timestep index was {max_ts + 1=}"
        )
    completed = max_ts + 1 == nsteps

    # Currently 1-based numbering corresponding to internal SHARC usage.
    # Ultimately aiming to replace numbers with labels ('S0', 'S1', ...),
    # but that has disadvantages in postprocessing.
    coords: dict | xr.Dataset = {
        "time": np.arange(nsteps),
        "state": (states := np.arange(1, nstates + 1)),
        "state2": states,
        "atom": np.arange(natoms),
        "direction": ["x", "y", "z"],
    }

    statecomb = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(states, 2), names=["from", "to"]),
        dim="statecomb",
    )

    coords = statecomb.merge(coords)

    default_sharc_attributes = get_default_input_attributes("sharc", loading_parameters)

    # TODO: FIXME: Use input names and input units and only apply attributes if there are attributes

    res = xr.Dataset(
        {
            "atXYZ": (
                template["atXYZ"],
                positions,
                default_sharc_attributes["atXYZ"],
            ),
            "energy": (
                template["energy"],
                energy,
                default_sharc_attributes["energy"],
            ),
            "e_kin": (
                template["e_kin"],
                e_kin,
                default_sharc_attributes["e_kin"],
            ),
            "dip_perm": (
                template["dip_perm"],
                dip_perm,
                default_sharc_attributes["dip_perm"],
            ),
            "dip_trans": (
                template["dip_trans"],
                dip_trans,
                default_sharc_attributes["dip_trans"],
            ),
            "sdiag": (
                template["sdiag"],
                sdiag,
                default_sharc_attributes["sdiag"],
            ),
            "astate": (template["astate"], astate, default_sharc_attributes["astate"]),
            "forces": (template["forces"], forces, default_sharc_attributes["forces"]),
            # 'has_forces': (['ts'], has_forces),
            "phases": (template["phases"], phases, default_sharc_attributes["phases"]),
            "nacs": (template["nacs"], nacs, default_sharc_attributes["nacs"]),
            "atNums": (template["atNums"], atNums, default_sharc_attributes["atNums"]),
            "atNames": (
                template["atNames"],
                atNames,
                default_sharc_attributes["atNames"],
            ),
            "state_names": (
                template["state_names"],
                state_names,
                default_sharc_attributes["state_names"],
            ),
            "state_types": (
                template["state_types"],
                state_types,
                (
                    default_sharc_attributes["state_types"]
                    if "state_types" in default_sharc_attributes
                    else {}
                ),
            ),
        },
        coords=coords,
        attrs={"max_ts": max_ts, "completed": completed},
    )

    for coord_name in res.coords:
        if coord_name in default_sharc_attributes:
            res[coord_name].attrs.update(default_sharc_attributes[str(coord_name)])

    res.attrs["input_format"] = "sharc"
    res.attrs["input_type"] = "dynamic"

    res.attrs["num_singlets"] = nsinglets
    res.attrs["num_doublets"] = ndoublets
    res.attrs["num_triplets"] = ntriplets

    if not completed:
        # Filter by index not by ts
        # res = res.sel(ts=res.ts <= res.attrs["max_ts"])
        res = res.isel(time=slice(0, max_ts + 1))

    return res


def parse_trajout_xyz(
    nsteps: int, f: TextIOWrapper
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read atom names, atom numbers and positions for each time step up until a maximum of `nsteps` from an `output.xyz` file and returm them.

    Args:
        nsteps (int): The maximum number of steps to be read
        f (TextIOWrapper): The input file wrapper providing the contents of an `output.xyz` file.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
            Tuple of (atom_names, atom_numbers, atom_positions) as numpy arrays.
            Only atom_positions has the first index indicate the time step, the second the atom and the third the direction.
            Other entries are 1d arrays.
    """
    first = next(f)
    assert first.startswith(" " * 6)
    natoms = int(first.strip())

    atNames = np.full((natoms), "")
    atNums = np.full((natoms), -1)
    atXYZ = np.full((nsteps, natoms, 3), np.nan)

    ts = 0

    for index, line in enumerate(f):
        if "t=" in line:
            assert (
                ts < nsteps
            ), f"Excess time step at ts={ts}, for a maximum nsteps={nsteps}"
            for atom in range(natoms):
                linecont = re.split(" +", next(f).strip())
                if ts == 0:
                    atNames[atom] = linecont[0]
                atNums[atom] = get_atom_number_from_symbol(linecont[0])
                atXYZ[ts, atom] = [float(n) for n in linecont[1:]]
            ts += 1

    return (atNames, atNums, atXYZ)


def parse_input_settings(f: TextIOWrapper) -> Dict[str, Any]:
    """Function to parse settings from the `input` file

    Args:
        f (TextIOWrapper): File wrapper providing the input file contents

    Returns:
        Dict[str, Any]: A key-value dictionary, where the keys are the first words in each line
    """
    settings = {}
    for line in f:
        parsed = line.strip().split()
        if len(parsed) == 2:
            settings[parsed[0]] = parsed[1]
        elif len(parsed) > 2:
            settings[parsed[0]] = parsed[1:]
        elif len(parsed) == 1:
            settings[parsed[0]] = True
    return settings
