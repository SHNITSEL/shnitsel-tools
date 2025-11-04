import pathlib
from shnitsel.io.helpers import LoadingParameters, make_uniform_path
from io import TextIOWrapper
import numpy as np
import pandas as pd
import xarray as xr
import logging
import os
import re
import math
from itertools import product, combinations
from glob import glob
from typing import Dict, List, NamedTuple, Any, Sequence, Set, Tuple
from tqdm.auto import tqdm
from shnitsel.io.helpers import (
    PathOptionsType,
    dip_sep,
    __atnum2symbol__,
    get_triangular,
    ConsistentValue,
    get_atom_number_from_symbol,
)
from shnitsel.io.xyz import get_dipoles_per_xyz
from shnitsel._contracts import needs
from shnitsel.units.definitions import get_default_input_attributes
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults

_re_grads = re.compile("[(](?P<nstates>[0-9]+)x(?P<natoms>[0-9]+)x3")
_re_nacs = re.compile("[(](?P<nstates>[0-9]+)x[0-9]+x(?P<natoms>[0-9]+)x3")


class IcondPath(NamedTuple):
    idx: int
    path: pathlib.Path
    # prefix: str | None


def nans(*dims):
    return np.full(dims, np.nan)


def list_iconds(
    iconds_path: str | os.PathLike = "./iconds/", glob_expr: str = "**/ICOND_*"
) -> List[IcondPath]:
    """Retrieve a list of all potential initial condition directories to be parsed given the input path and matching patter.

    Args:
        iconds_path (str | os.PathLike, optional): The path where to look for initial condition directories. Defaults to './iconds/'.
        glob_expr (str, optional): The pattern for finding initial conditions. Defaults to '**/ICOND_*'.

    Raises:
        FileNotFoundError: If no directories match the pattern.

    Returns:
        List[IcondPath]: The list of Tuples of the parsed ID and the full path of the initial conditions
    """
    path_obj: pathlib.Path = make_uniform_path(iconds_path)

    dirs = list(
        path_obj.glob(
            glob_expr,
            # recursive=True
        )
    )
    if len(dirs) == 0:
        raise FileNotFoundError(
            f"The search '{glob_expr}' didn't match any directories "
            f"under {iconds_path=} "
            f"relative to working directory '{os.getcwd()}'"
        )
    icond_paths = sorted(
        [
            candidate_directory
            for candidate_directory in dirs
            if (qm_out := candidate_directory / "QM.out").exists() and qm_out.is_file()
        ],
        key=lambda x: x.name,
    )

    # FIXME: This fails if the pattern of the path changes
    return [IcondPath(int(ipath.name[6:]), ipath) for ipath in icond_paths]


def dims_from_QM_out(f: TextIOWrapper) -> Tuple[int | None, int | None]:
    """Function to also read the relevant dimensions (number of atoms, number of states) from WM.out

    Args:
        f (TextIOWrapper): The QM.out file to read from

    Returns:
        Tuple[int,int]: First the number of states, second the number of atoms or respectively None if not found.
    """
    # faced with redundancy, use it to ensure consistency
    nstates = ConsistentValue("nstates", weak=True)
    natoms = ConsistentValue("natoms", weak=True)

    for index, line in enumerate(f):
        if line.startswith("! 1 Hamiltonian Matrix"):
            nstates.v = int(next(f).split(" ")[0])
        elif line.startswith("! 2 Dipole Moment Matrices"):
            dim = re.split(" +", next(f).strip())
            nstates.v = int(dim[0])
        elif line.startswith("! 3 Gradient Vectors"):
            info = _re_grads.search(line)
            assert info is not None
            nstates.v, natoms.v = map(int, info.group("nstates", "natoms"))
        elif line.startswith("! 5 Non-adiabatic couplings"):
            info = _re_nacs.search(line)
            assert info is not None
            nstates.v, natoms.v = map(int, info.group("nstates", "natoms"))

    return nstates.v, natoms.v


def dims_from_QM_log(log: TextIOWrapper) -> Tuple[int, int, int, int, int]:
    """Function to retrieve the listed number of states and the number of atoms from the Qm.log file of initial conditions

    Args:
        log (TextIOWrapper): Input file handle to read the log file contents from

    Raises:
        ValueError: If the log file contains an inconsistently structured Line about states but does neither specify Singlet nor Triplet states.

    Returns:
        Tuple[int,int,int,int,int]: First the number of states, then the number of atoms. This is followed by the number of singlet, doublet and triplet states (if available). If a value is not available, it will default to 0.
    """
    nstates = ConsistentValue("nstates", weak=True)
    nstates_singlet = ConsistentValue("nstates_singlet", weak=True)
    nstates_doublet = ConsistentValue("nstates_doublet", weak=True)
    nstates_triplet = ConsistentValue("nstates_triplet", weak=True)

    natoms = ConsistentValue("natoms", weak=True)
    for line in log:
        if line.startswith("States:"):
            linecont = line.strip().split()
            if "Singlet" in linecont and "Triplet" not in linecont:
                nsinglets = int(linecont[1])
                ntriplets = 0
            elif "Singlet" in linecont and "Triplet" in linecont:
                nsinglets = int(linecont[1])
                ntriplets = int(linecont[3])
            elif "Triplet" in linecont and "Singlet" not in linecont:
                ntriplets = int(linecont[1])
                nsinglets = 0
            else:
                raise ValueError(f"Invalid State line in QM.log: {line}")

            # calculate total number of states
            nstates.v = nsinglets + (3 * ntriplets)
            nstates_singlet.v = nsinglets
            nstates_triplet.v = ntriplets

        elif line.startswith("Found Geo!"):
            linecont = re.split(" ", line.strip())
            natoms.v = int(linecont[-1][0:-1])

    num_states = nstates.v if nstates.v is not None else 0
    num_atoms = natoms.v if natoms.v is not None else 0
    num_singlets = nstates_singlet.v if nstates_singlet.v is not None else 0
    num_doublets = nstates_doublet.v if nstates_doublet.v is not None else 0
    num_triplets = nstates_triplet.v if nstates_triplet.v is not None else 0

    return num_states, num_atoms, num_singlets, num_doublets, num_triplets


def check_dims(pathlist: Sequence[pathlib.Path]) -> Tuple[int, int, int, int, int]:
    """Function to obtain the number of atoms and states across all input paths.

    Will only return the tuple of (number_states, number_atoms) if these numbers are consistent across all paths.
    Otherwise, an error will be raised.

    Args:
        pathlist (Sequence[pathlib.Path]): The list of paths belonging to the same system to check for consistent dimensions

    Raises:
        FileNotFoundError: If the number of valid input paths in pathlist is zero, a FileNotFoundError is raised
        ValueError: If the number of states and the number of atoms does not agree across all systems a ValueError is raised

    Returns:
        Tuple[int, int, int, int, int]: The number of states and the number of atoms, then the number of singlets, doublets and triplets in this order
    """

    nstates = ConsistentValue("nstates", ignore_none=True)
    nstates_singlet = ConsistentValue("nstates_singlet", weak=True)
    nstates_doublet = ConsistentValue("nstates_doublet", weak=True)
    nstates_triplet = ConsistentValue("nstates_triplet", weak=True)

    natoms = ConsistentValue("natoms", ignore_none=True)
    for path in pathlist:
        try:
            with open(path / "QM.out") as f:
                nstates.v, natoms.v = dims_from_QM_out(f)
        except FileNotFoundError:
            pass
        try:
            with open(path / "QM.log") as f:
                (
                    nstates.v,
                    natoms.v,
                    nstates_singlet.v,
                    nstates_doublet.v,
                    nstates_triplet.v,
                ) = dims_from_QM_log(f)
        except FileNotFoundError:
            pass

    if not nstates.defined or not natoms.defined:
        raise FileNotFoundError(
            "Pathlist empty or no valid path found within pathlist for initial condition input"
        )

    num_singlets = nstates_singlet.v if nstates_singlet.v is not None else 0
    num_doublets = nstates_doublet.v if nstates_doublet.v is not None else 0
    num_triplets = nstates_triplet.v if nstates_triplet.v is not None else 0

    return nstates.v, natoms.v, num_singlets, num_doublets, num_triplets


def dir_of_iconds(
    path: PathOptionsType,
    *,
    levels: int = 1,
    id_subset: Set[int] | None = None,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset:
    """Function to retrieve all initial conditions from subdirectories of the provided `path`.


    Args:
        path (str | os.PathLike, optional): The path to a directory containing the directories with initial conditions. Defaults to './iconds/'.
        levels (int, optional): Currently unused. Intended to denote how many levels down initial conditions will be looked for. Defaults to 1.
        id_subset (Set[int]  |  None, optional): An optional set of ids to restrict the anaylsis to. No initial conditions with an id outside of this set will be loaded. Defaults to None.

    Returns:
        xr.Dataset: The resulting dataset containing all of the initial conditions as separate trajectories
    """
    initial_condition_paths: List[IcondPath] = list_iconds(path)
    if id_subset is not None:
        initial_condition_paths = [
            icond for icond in initial_condition_paths if icond.idx in id_subset
        ]

    return read_iconds_multi_directory(initial_condition_paths, loading_parameters)


def finalize_icond_dataset(
    dataset: xr.Dataset, loading_parameters: LoadingParameters
) -> xr.Dataset:
    """Function to expand the initial conditions dataset with a time dimension.

    Also sets the default unit on the time dimension

    Args:
        dataset (xr.Dataset): The initial conditions dataset. Should not have a "time" dimension yet.
        loading_parameters (LoadingParameters): Loading parameters to override units

    Returns:
        xr.Dataset: The modified dataset
    """
    isolated_keys = [
        "atNames",
        "atNums",
        "state",
        "statecomb",
        "state_names",
        "state_types",
    ]

    if "time" not in dataset.coords:
        dataset_res = dataset.set_coords(isolated_keys)
        dataset_res = dataset_res.expand_dims("time")
        dataset_res = dataset_res.assign_coords(time=("time", [0.0]))

        default_sharc_attributes = get_default_input_attributes(
            "sharc", loading_parameters
        )
        dataset_res["time"].attrs.update(default_sharc_attributes["time"])

    # Set completed flag
    dataset_res.atts["completed"] = True
    return dataset_res


def create_icond_dataset(
    indices: List[int] | None,
    nstates: int,
    natoms: int,
    loading_parameters: LoadingParameters | None,
    **kwargs,
) -> xr.Dataset:
    """Function to initialize an `xr.Dataset` with appropriate variables and coordinates to acommodate loaded data.

    Args:
        indices (List[int]): List of indices for the different initial conditions
        nstates (int): Number of states within the datasets
        natoms (int): The number of atoms within the datasets

    Returns:
        xr.Dataset: An xarray Dataset with appropriately sized DataArrays and coordinates also including default attributes for all variables.
    """
    template = {
        "energy": ["state"],
        "dip_all": ["state", "state2", "direction"],
        "dip_perm": ["state", "direction"],
        "dip_trans": ["statecomb", "direction"],
        "forces": ["state", "atom", "direction"],
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
        "dip_all": np.nan,
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

    if indices is not None and len(indices) > 0:
        niconds = len(indices)
        for varname, dims in template.items():
            dims.insert(0, "icond")
    else:
        niconds = 0

    dim_lengths = {
        "icond": niconds,
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
        "has_forces": (
            "icond",
            x if (x := kwargs.get("has_forces")) is not None else nans(niconds),
        ),
    }

    if indices is not None:
        coords["icond"] = indices

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
    default_sharc_attributes = get_default_input_attributes("sharc", loading_parameters)

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

    res_dataset.attrs["input_format"] = "sharc"
    res_dataset.attrs["input_type"] = "static"

    return res_dataset


# TODO: FIXME: Make function read a single initial condition set and use the generic aggregation functions to combine them.
def read_iconds_individual(
    path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read initial a single initial condition directory into a Dataset with standard shnitsel annotations and units

    Args:
        path (PathOptionsType): The path to a initial conditions directory
        loading_parameters (LoadingParameters | None, optional): Parameter settings for e.g. standard units or state names.

    Returns:
        xr.Dataset: The Dataset object containing all of the loaded data from the initial condition in default shnitsel units
    """
    # TODO: FIXME: use loading_parameters to configure units and state names
    logging.info("Ensuring consistency of ICONDs dimensions")
    nstates, natoms, nsinglets, ndoublets, ntriplets = check_dims([path])
    iconds = create_icond_dataset(
        None, nstates, natoms, loading_parameters=loading_parameters
    )

    # Set information on the singlet, doublet and triplet states, if available
    iconds["state_types"][:nsinglets] = 1
    iconds["state_types"][nsinglets : nsinglets + 2 * ndoublets] = 2
    iconds["state_types"][nsinglets + 2 * ndoublets :] = 3

    iconds.attrs["num_singlets"] = nsinglets
    iconds.attrs["num_doublets"] = ndoublets
    iconds.attrs["num_triplets"] = ntriplets

    logging.info("Reading ICONDs data into Dataset...")

    with open(path / "QM.out") as f:
        parse_QM_out(f, out=iconds, loading_parameters=loading_parameters)

    try:
        with open(path / "QM.log") as f:
            parse_QM_log_geom(f, out=iconds)
    except FileNotFoundError:
        # This should be an error. We probably cannot recover from this and action needs to be taken
        logging.error(
            f"""no QM.log file found in {path}.
            This is currently used to determine geometry.
            Eventually, user-inputs will be accepted as an alternative.
            See https://github.com/SHNITSEL/db-workflow/issues/3"""
        )
        logging.warning(
            f"No positional information found in {path}, the loaded trajectory does not contain positional data 'atXYZ'."
        )
        return None

    return convert_all_units_to_shnitsel_defaults(
        finalize_icond_dataset(iconds, loading_parameters=loading_parameters)
    )


def read_iconds_multi_directory(
    pathlist: List[IcondPath],
    indices: List[int] | None = None,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset:
    """Function to read initial condition directories into a Dataset with standard shnitsel annotations and units

    Args:
        pathlist (List[IcondPath]): Lists of Initial condition paths from which we can parse the results into the dataset
        indices (List[int] | None, optional): Optional. The list of indices of initial conditions we want to limit the parsing to. Defaults to None.
        loading_parameters (LoadingParameters | None, optional): Parameter settings for e.g. standard units or state names.

    Returns:
        xr.Dataset: The Dataset object containing all of the loaded data in default shnitsel units
    """
    # TODO: FIXME: use loading_parameters to configure units and state names
    logging.info("Ensuring consistency of ICONDs dimensions")
    nstates, natoms, nsinglets, ndoublets, ntriplets = check_dims(
        [path for _, path in pathlist]
    )
    logging.info("Allocating Dataset for ICONDs")
    if indices is None:
        indices = [p.idx for p in pathlist]
    iconds = create_icond_dataset(
        indices, nstates, natoms, loading_parameters=loading_parameters
    )

    # Set information on the singlet, doublet and triplet states, if available
    iconds["state_types"][:nsinglets] = 1
    iconds["state_types"][nsinglets : nsinglets + 2 * ndoublets] = 2
    iconds["state_types"][nsinglets + 2 * ndoublets :] = 3

    iconds.attrs["num_singlets"] = nsinglets
    iconds.attrs["num_doublets"] = ndoublets
    iconds.attrs["num_triplets"] = ntriplets

    logging.info("Reading ICONDs data into Dataset...")

    for icond_index, path in tqdm(pathlist):
        with open(path / "QM.out") as f:
            parse_QM_out(
                f,
                out=iconds.sel(icond=icond_index),
                loading_parameters=loading_parameters,
            )

    for icond_index, path in tqdm(pathlist):
        try:
            with open(path / "QM.log") as f:
                parse_QM_log_geom(f, out=iconds.sel(icond=icond_index))
        except FileNotFoundError:
            # This should be an error. We probably cannot recover from this and action needs to be taken
            logging.error(
                f"""no QM.log file found in {path}.
                This is currently used to determine geometry.
                Eventually, user-inputs will be accepted as an alternative.
                See https://github.com/SHNITSEL/db-workflow/issues/3"""
            )
            logging.warning(
                f"No positional information found in {path}, the loaded trajectory does not contain full positional data 'atXYZ'."
            )
            return None

    return convert_all_units_to_shnitsel_defaults(
        finalize_icond_dataset(iconds, loading_parameters=loading_parameters)
    )


def parse_QM_log(log: TextIOWrapper) -> Dict[str, Any]:
    """Function to parse main information from the QM.log file

    Args:
        log (TextIOWrapper): Input file wrapper to read the information from

    Raises:
        ValueError: If there are neither singlet nor triplet states listed in the QM.log file

    Returns:
        Dict[str, Any]: Dictionary with key information about the system
    """
    info: Dict[str, Any] = {}
    for line in log:
        if line.startswith("States:"):
            linecont = re.split(" +|\t", line.strip())
            if "Singlet" in linecont and "Triplet" not in linecont:
                nsinglets = int(linecont[2])
                ntriplets = 0
            elif "Singlet" in linecont and "Triplet" in linecont:
                nsinglets = int(linecont[2])
                ntriplets = int(linecont[5])
            elif "Triplet" in linecont and "Singlet" not in linecont:
                ntriplets = int(linecont[2])
                nsinglets = 0
            else:
                raise ValueError(
                    "QM.log file is malformed. States have neither singlet nor triplet states listed"
                )

            # calculate total number of states
            nstates = nsinglets + (3 * ntriplets)

            info["nStates"] = nstates
            info["nSinglets"] = nsinglets
            info["nTriplets"] = ntriplets
            nnacs = int(nsinglets * (nsinglets - 1) / 2) + int(
                ntriplets * (ntriplets - 1) / 2
            )
            info["nNACS"] = nnacs
            info["nDipoles"] = int(nsinglets + ntriplets + nnacs)

        elif line.startswith("Method:"):
            linecont = re.split(" +|\t", line.strip())
            method = linecont[2]

            info["method"] = method

        elif line.startswith("Found Geo!"):
            linecont = re.split(" ", line.strip())
            natom = int(linecont[-1][0:-1])

            info["nAtoms"] = natom

        elif line.startswith("Geometry in Bohrs:"):
            # NB. Geometry is indeed in bohrs!
            atnames = []
            atxyz = np.zeros((natom, 3))
            for i in range(natom):
                geometry_line = re.split(" +", next(log).strip())
                atnames.append(geometry_line[0])
                atxyz[i] = [float(geometry_line[j]) for j in range(1, 4)]

            info["atNames"] = atnames
            info["atNums"] = [get_atom_number_from_symbol(n) for n in atnames]
            info["atXYZ"] = atxyz

    return info


def parse_QM_log_geom(f: TextIOWrapper, out: xr.Dataset):
    """Read geometry into an xr.Dataset object from the provided file input stream `f`.

    f must be the contents of a `QM.log` file.

    Args:
        f (TextIOWrapper): File wrapper for a `QM.log` file's contents
        out (xr.Dataset): The dataset to write the resulting geometry to
    """

    # NB. Geometry is indeed in bohrs!
    while not next(f).startswith("Geometry in Bohrs:"):
        pass

    for i in range(out.sizes["atom"]):
        geometry_line = next(f).strip().split()
        atom_symbol = geometry_line[0].strip()
        out["atNames"][i] = atom_symbol
        out["atNums"][i] = get_atom_number_from_symbol(atom_symbol)
        out["atXYZ"][i] = map(float, geometry_line[1:4])


def parse_QM_out(
    f: TextIOWrapper,
    out: xr.Dataset | None = None,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset | None:
    """Function to read all information about forces, energies, dipoles, nacs, etc. from initial condition QM.out files.

    if ``out=None`` is provided, a new Dataset is constructed and returned by the function

    Args:
        f (TextIOWrapper): File input of the QM.out file to parse the data from
        out (xr.Dataset  |  None, optional): Optional target Dataset to write the loaded data into. Defaults to None.
        loading_parameters (LoadingParameters,optional): Optional loading parameters to override variable mappings and units.

    Returns:
        xr.Dataset | None: If a new Dataset was constructed instead of being written to `out`, it will be returned.
    """
    res: xr.Dataset | dict[str, np.ndarray]
    if out is not None:
        # write data directly into dataset
        res = out
    else:
        # write data as ndarrays into dict, then make dataset after parsing
        res = {}

    res["has_forces"] = np.array([0])
    nstates = ConsistentValue("nstates")
    natoms = ConsistentValue("natoms")

    for index, line in enumerate(f):
        if line.startswith("! 1 Hamiltonian Matrix"):
            # get number of states from dimensions of Hamiltonian
            nstates.v = int(next(f).split(" ")[0])
            if out is None:
                res["energy"] = nans(nstates.v)

            for istate in range(nstates.v):
                energyline = re.split(" +", next(f).strip())
                res["energy"][istate] = float(energyline[2 * istate])

        elif line.startswith("! 2 Dipole Moment Matrices"):
            dim = re.split(" +", next(f).strip())
            n = int(dim[0])
            m = int(dim[1])

            if out is None:
                res["dip_all"] = nans(n, m, 3)
                res["dip_perm"] = nans(n, 3)
                res["dip_trans"] = nans(math.comb(n, 2), 3)

            res["dip_all"][:, :, 0] = get_dipoles_per_xyz(f, n, m)
            next(f)
            res["dip_all"][:, :, 1] = get_dipoles_per_xyz(f, n, m)
            next(f)
            res["dip_all"][:, :, 2] = get_dipoles_per_xyz(f, n, m)

            res["dip_perm"][:], res["dip_trans"][:] = dip_sep(np.array(res["dip_all"]))

        elif line.startswith("! 3 Gradient Vectors"):
            res["has_forces"] = np.array([1])

            search_res = _re_grads.search(line)
            assert search_res is not None
            get_dim = search_res.group
            nstates.v = int(get_dim("nstates"))
            natoms.v = int(get_dim("natoms"))

            if out is None:
                res["forces"] = nans(nstates.v, natoms.v, 3)

            for istate in range(nstates.v):
                next(f)
                for atom in range(natoms.v):
                    res["forces"][istate][atom] = [
                        float(entry) for entry in next(f).strip().split()
                    ]

        elif line.startswith("! 5 Non-adiabatic couplings"):
            search_res = _re_nacs.search(line)
            assert search_res is not None
            get_dim = search_res.group
            nstates.v = int(get_dim("nstates"))
            natoms.v = int(get_dim("natoms"))

            if out is None:
                res["nacs"] = nans(math.comb(nstates.v, 2), natoms.v, 3)

            nacs_all = nans(nstates.v, nstates.v, natoms.v, 3)

            for bra, ket in product(range(nstates.v), range(nstates.v)):
                # TODO info currently unused, but keep the `next(f)` no matter what!
                nac_multi = int(re.split(" +", next(f).strip())[-1])  # noqa: F841

                for atom in range(natoms.v):
                    nacs_line = re.split(" +", next(f).strip())
                    nacs_all[bra, ket, atom] = [float(n) for n in nacs_line]

            # all nacs, i.e., nacs of all singlet and triplet states
            # all diagonal elements are zero (self-coupling, e.g. S1 and S1)
            # off-diagonal elements contain couplings of different states (e.g. S0 and S1)
            # in principle one has here the full matrix for the nacs between all singlet and triplet states
            # in the following we extract only the upper triangular elements of the matrix

            res["nacs"][:] = get_triangular(nacs_all)

        elif line.startswith("! 6 Overlap matrix"):
            nlines = int(re.split(" +", next(f).strip())[0])
            assert nlines == nstates.v

            found_overlap = False
            phasevector = np.ones((nlines))

            wvoverlap = np.zeros((nlines, nlines))
            for j in range(nlines):
                linecont = [float(n) for n in re.split(" +", next(f).strip())]
                vec = [n for n in linecont[::2]]
                assert len(vec) == nlines
                wvoverlap[j] = vec

            for istate in range(nlines):
                if np.abs(wvoverlap[istate, istate]) >= 0.5:
                    found_overlap = True
                    if wvoverlap[istate, istate] >= 0.5:
                        res["phases"][istate] = +1
                    else:
                        res["phases"][istate] = -1

            if found_overlap:
                res["phases"][:] = phasevector
                pass

        elif line.startswith("! 8 Runtime"):
            next(f)

    if out is None:
        if not res["has_forces"]:
            res["forces"] = nans(natoms.v, 3)

        assert isinstance(res, dict)
        return create_icond_dataset(
            indices=None,
            nstates=nstates.v,
            natoms=natoms.v,
            loading_parameters=loading_parameters,
            **res,
        )

        # xr.Dataset(
        #     {
        #         'energy': (['state'], energy),
        #         'dip_all': (['state', 'state2', 'direction'], dip_all),
        #         'dip_perm': (['state', 'direction'], dip_perm),
        #         'dip_trans': (['statecomb', 'direction'], dip_trans),
        #         'forces': (['atom', 'direction'], forces),
        #         'has_forces': ([], has_forces),
        #         'phases': (['state'], phases),
        #         'nacs': (['statecomb', 'atom', 'direction'], nacs)
        #     },
        #     coords={
        #         'state': np.arange(1, nstates.v+1),
        #         'state2': np.arange(1, nstates.v+1),
        #         'atom': np.arange(natoms.v),
        #         'statecomb': np.arange(math.comb(nstates.v, 2)),
        #         'direction': ['x', 'y', 'z']
        #     }
        # )
    else:
        # all the data has already been written to `out`
        # no need to return anything
        return None


@needs(dims={"icond"}, coords={"icond"}, not_dims={"time"})
def iconds_to_frames(iconds: xr.Dataset) -> xr.Dataset:
    """Function to convert the `icond` coordinate into a `trajid` coordinate and also build a combined `frame`+`time` multiindex as `frame`

    Will exempt atNames and atNums from the rearrangement.

    Args:
        iconds (xr.Dataset): input dataset to replace the dimension in

    Raises:
        ValueError: Raised if at least one array has size 0 in one coordinate

    Returns:
        xr.Dataset: The transformed dataset
    """
    for name, var in iconds.data_vars.items():
        shape = var.data.shape
        if 0 in shape:
            raise ValueError(
                f"Variable '{name}' has shape {shape} which contains 0. "
                "Please remove this variable before converting to frames. "
                "Note: An empty variable could indicate a problem with parsing."
            )

    # if 'atNames' in iconds.data_vars and 'atNames' not in iconds.coords:
    #    iconds = iconds.assign_coords(atNames=iconds.atNames)

    isolated_keys = ["atNames", "atNums", "state_names", "state_types"]

    res = iconds.rename_dims(icond="trajid").rename_vars(icond="trajid")

    for var in res.data_vars:
        if var not in isolated_keys:
            res[var] = res[var].expand_dims("time")

    return res.assign_coords(time=("time", [0.0])).stack(frame=["trajid", "time"])
