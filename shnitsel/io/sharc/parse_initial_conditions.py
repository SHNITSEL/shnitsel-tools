import pathlib
from shnitsel.io.shared.helpers import LoadingParameters, make_uniform_path
from io import TextIOWrapper
import numpy as np
import xarray as xr
import logging
import os
import re
from itertools import product
from typing import Dict, List, NamedTuple, Any, Tuple
from shnitsel.io.shared.helpers import (
    PathOptionsType,
    dip_sep,
    get_triangular,
    ConsistentValue,
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
from shnitsel.io.xyz import get_dipoles_per_xyz
from shnitsel._contracts import needs

_re_grads = re.compile("[(](?P<nstates>[0-9]+)x(?P<natoms>[0-9]+)x3")
_re_nacs = re.compile("[(](?P<nstates>[0-9]+)x[0-9]+x(?P<natoms>[0-9]+)x3")
_re_state_map_entry = re.compile(
    r"(?P<state_id>\d+): +\[(?P<multiplicity>\d+), +(?P<n_mult>\d+), *(?P<unknown>[^\[]*)\]"
)


class IcondPath(NamedTuple):
    idx: int
    path: pathlib.Path
    # prefix: str | None


def nans(*dims):
    return np.full(dims, np.nan)


def list_iconds(
    iconds_path: str | os.PathLike = "./iconds/", glob_expr: str = "**/ICOND_*"
) -> list[IcondPath]:
    """Retrieve a list of all potential initial condition directories to be parsed given the input path and matching patter.

    Parameters
    ----------
    iconds_path : str | os.PathLike, optional
        The path where to look for initial condition directories. Defaults to './iconds/'.
    glob_expr : str, optional
        The pattern for finding initial conditions. Defaults to '**/ICOND_*'.

    Raises:
        FileNotFoundError: If no directories match the pattern.

    Returns
    -------
    list[IcondPath]
        The list of Tuples of the parsed ID and the full path of the initial conditions
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


def dims_from_QM_out(f: TextIOWrapper) -> tuple[int | None, int | None]:
    """Function to also read the relevant dimensions (number of atoms, number of states) from WM.out

    Parameters
    ----------
    f : TextIOWrapper
        The QM.out file to read from

    Returns
    -------
    tuple[int,int]
        First the number of states, second the number of atoms or respectively None if not found.
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


def dims_from_QM_log(log: TextIOWrapper) -> tuple[int, int, int, int, int]:
    """Function to retrieve the listed number of states and the number of atoms from the Qm.log file of initial conditions

    Parameters
    ----------
    log : TextIOWrapper
        Input file handle to read the log file contents from

    Raises
    ------
    ValueError
        If the log file contains an inconsistently structured Line about states but does neither specify Singlet nor Triplet states.

    Returns
    -------
    tuple[int,int,int,int,int]
        First the number of states, then the number of atoms. This is followed by the number of singlet, doublet and triplet states (if available). If a value is not available, it will default to 0.
    """
    nstates = ConsistentValue("nstates", weak=True)
    nstates_singlet = ConsistentValue("nstates_singlet", weak=True)
    nstates_doublet = ConsistentValue("nstates_doublet", weak=True)
    nstates_triplet = ConsistentValue("nstates_triplet", weak=True)

    natoms = ConsistentValue("natoms", weak=True)
    for line in log:
        line = line.strip()
        if line.startswith("States:"):
            linecont = line[len("States:") :].strip().split()

            num_parts = len(linecont)
            if num_parts % 2 != 0:
                logging.warning(
                    "Found `States:` line with odd number of entries. Expected pairs `<number_of_states> <state_multiplicity>`."
                )
            nsinglets = 0
            ndoublets = 0
            ntriplets = 0
            for i in range(0, num_parts, 2):
                type_str = linecont[i + 1].strip().lower()
                if type_str == "singlet":
                    nsinglets = int(linecont[i])
                elif type_str == "doublet":
                    ndoublets = int(linecont[i])
                elif type_str == "doublet":
                    ntriplets = int(linecont[i])
                else:
                    raise ValueError(
                        f"Invalid State line in QM.log: {line}. Found unknown multiplicity: `{type_str}`."
                    )

            # calculate total number of states
            nstates.v = nsinglets + (2 * ndoublets) + (3 * ntriplets)
            nstates_singlet.v = nsinglets
            nstates_doublet.v = ndoublets
            nstates_triplet.v = ntriplets

        elif line.startswith("Found Geo!"):
            linecont = re.split(" ", line.strip())
            natoms.v = int(linecont[-1][0:-1])
        elif line.startswith("statemap:"):
            # state_id: [multiplicity, index_mult, ?]
            # {1: [1, 1, 0.0], 2: [1, 2, 0.0]}
            state_map = line[len("statemap:") :].strip()
            # print(state_map)

            matches = _re_state_map_entry.findall(state_map)

            if matches:
                states_count = [0 for _ in range(0, 3)]
                for match in matches:
                    state_id, multiplicity, n_mult, unknown = match

                    state_id = int(state_id)
                    multiplicity = int(multiplicity)
                    n_mult = int(n_mult)
                    unknown = float(unknown)

                    if multiplicity >= 4:
                        logging.warning(
                            "Found state of above triplet multiplicity: %(state_id)d (full state descriptor: {(%(multiplicity)d, %(n_mult)d, %(unknown)f)}). \n State will be ignored.",
                            {
                                'state_id': state_id,
                                'multiplicity': multiplicity,
                                'n_mult': n_mult,
                                'unknown': unknown,
                            },
                        )
                    else:
                        states_count[multiplicity - 1] += 1

                nsinglets, ndoublets, ntriplets = states_count

                nstates.v = nsinglets + (2 * ndoublets) + (3 * ntriplets)
                nstates_singlet.v = nsinglets
                nstates_doublet.v = ndoublets
                nstates_triplet.v = ntriplets

    num_states = nstates.v if nstates.v is not None else 0
    num_atoms = natoms.v if natoms.v is not None else 0
    num_singlets = nstates_singlet.v if nstates_singlet.v is not None else 0
    num_doublets = nstates_doublet.v if nstates_doublet.v is not None else 0
    num_triplets = nstates_triplet.v if nstates_triplet.v is not None else 0

    return num_states, num_atoms, num_singlets, num_doublets, num_triplets


def check_dims(pathlist: List[pathlib.Path]) -> tuple[int, int, int, int, int]:
    """Function to obtain the number of atoms and states across all input paths.

    Will only return the tuple of (number_states, number_atoms) if these numbers are consistent across all paths.
    Otherwise, an error will be raised.

    Parameters
    ----------
    pathlist : Sequence[pathlib.Path]
        The list of paths belonging to the same system to check for consistent dimensions

    Raises
    ------
    FileNotFoundError
        If the number of valid input paths in pathlist is zero, a FileNotFoundError is raised
    ValueError
        If the number of states and the number of atoms does not agree across all systems a ValueError is raised

    Returns
    -------
    tuple[int, int, int, int, int]
        The number of states and the number of atoms, then the number of singlets, doublets and triplets in this order
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
            with open(path / "QM.in") as f:
                info = parse_QM_in(f)

                if "num_atoms" in info:
                    natoms.v = int(info["num_atoms"])
                if "num_states" in info:
                    nstates.v = info["num_states"]
                if "num_singlets" in info:
                    nstates_singlet.v = info["num_singlets"]
                if "num_doublets" in info:
                    nstates_doublet.v = info["num_doublets"]
                if "num_triplets" in info:
                    nstates_triplet.v = info["num_triplets"]
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


def finalize_icond_dataset(
    dataset: xr.Dataset,
    loading_parameters: LoadingParameters,
    default_format_attributes: dict[str, dict[str, Any]],
) -> xr.Dataset:
    """Function to expand the initial conditions dataset with a time dimension.

    Also sets the default unit on the time dimension based on `default_format_attributes`.

    Parameters
    ----------
    dataset : xr.Dataset)
        The initial conditions dataset. Should not have a "time" dimension yet.
    loading_parameters : LoadingParameters)
        Loading parameters to override units
    default_format_attributes : dict[str, dict[str, Any]]
        Default attributes to set on variables, mostly used to set the time dimension attributes

    Returns
    -------
    xr.Dataset
        The modified dataset
    """

    if "time" not in dataset.coords:
        dataset_res = dataset.expand_dims("time")
        dataset_res = dataset_res.assign_coords(time=("time", [0.0]))

        dataset_res["time"].attrs.update(default_format_attributes["time"])
        mark_variable_assigned(dataset_res["time"])
    else:
        dataset_res = dataset

    # Set completed flag
    dataset_res.attrs["completed"] = True
    return dataset_res


def read_iconds_individual(
    path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Function to read initial a single initial condition directory into a Dataset with standard shnitsel annotations and units

    Parameters
    ----------
    path : PathOptionsType
        The path to a initial conditions directory
    loading_parameters : LoadingParameters | None, optional
        Parameter settings for e.g. standard units or state names.

    Returns
    -------
    xr.Dataset
        The Dataset object containing all of the loaded data from the initial condition in default shnitsel units
    """
    from ...units.definitions import length

    path_obj: pathlib.Path = make_uniform_path(path)  # type: ignore
    # Read settings and initial setup from QM.in
    qm_in_path = path_obj / "QM.in"
    info = {}
    if qm_in_path.is_file():
        with open(qm_in_path) as f:
            info = parse_QM_in(f)

    # logging.info("Ensuring consistency of ICONDs dimensions")
    nstates, natoms, nsinglets, ndoublets, ntriplets = check_dims([path_obj])
    # logging.debug(
    #     f"Found {nstates} States, among which: S/D/T = {nsinglets}/{ndoublets}/{ntriplets}"
    # )

    # NOTE: Currently no way to determine the version of SHARC that wrote the iconds from QM.in and QM.out. only set from QM.log
    sharc_version = "unknown"

    # Create dataset
    iconds, default_format_attributes = create_initial_dataset(
        0, nstates, natoms, "sharc", loading_parameters
    )

    # logging.info("Reading ICONDs data into Dataset...")

    with open(path_obj / "QM.out") as f:
        parse_QM_out(f, out=iconds, loading_parameters=loading_parameters)

        # if we have found the version, use it.
    try:
        with open(path_obj / "QM.log") as f:
            parse_QM_log_geom(f, out=iconds)

        if "input_format_version" in iconds:
            sharc_version = iconds.attrs["input_format_version"]
    except FileNotFoundError:
        # This should be an error. We probably cannot recover from this and action needs to be taken
        logging.info(
            """no `QM.log` file found in %(path)s. 
            This is mainly used to determine geometry.\n
            Attempting to read from `QM.in` instead """,
            {'path': path},
            # Eventually, user-inputs will be accepted as an alternative.
            # See https://github.com/SHNITSEL/db-workflow/issues/3"""
        )

        try:
            # TODO: FIXME: Figure out unit of positions in QM.in
            if "atNames" in info:
                iconds["atNames"][:] = (atnames := info["atNames"])
                mark_variable_assigned(iconds.atNames)
                iconds["atNums"][:] = [get_atom_number_from_symbol(n) for n in atnames]
                mark_variable_assigned(iconds.atNums)
            if "atXYZ" in info:
                iconds["atXYZ"][:, :] = info["atXYZ"]
                if "unit" in info:
                    # We should set the unit accordingly if a unit is specified in QM.in
                    unit_name = info["unit"].lower()
                    if unit_name == "angstrom":
                        iconds.atXYZ.attrs["units"] = length.Angstrom
                    elif unit_name == "bohr":
                        iconds.atXYZ.attrs["units"] = length.Bohr
                    else:
                        logging.warning(
                            "Unsupported input length/distance unit in QM.in: %(unit)s. Unit on the position is assumed to be of unit %(assumed_unit)s",
                            {
                                'unit': unit_name,
                                'assumed_unit': default_format_attributes['atXYZ'][
                                    'units'
                                ],
                            },
                        )

                mark_variable_assigned(iconds.atXYZ)
        except FileNotFoundError:
            logging.warning(
                "No positional information found in %(path)s/QM.in nor in %(path)s/QM.log, the loaded trajectory does not contain positional data 'atXYZ'.",
                {'path': path},
            )

    # iconds.attrs["delta_t"] = 0.0
    # iconds.attrs["t_max"] = 0.0
    # iconds.attrs["max_ts"] = 1

    # Set all settings we require to be present on the trajectory
    required_settings = RequiredTrajectorySettings(
        0.0,
        0.0,
        1,
        True,
        "sharc",
        "static",
        sharc_version,
        nsinglets,
        ndoublets * 2,  # Factor 2 to make consistent with other tools
        ntriplets * 3,  # Factor 2 to make consistent with other tools
    )

    iconds = assign_required_settings(iconds, required_settings)

    optional_settings = OptionalTrajectorySettings(
        has_forces=is_variable_assigned(iconds["forces"]),
        misc_input_settings={"QM.in": info} if len(info) > 0 else None,
        est_level=info["method"] if "method" in info else None,
        theory_basis_set=info["basis"] if "basis" in info else None,
    )

    if "states" in info:
        multiplicities = info["states"]

        charges = None
        if "charge" in info:
            charges = info["charge"]
        else:
            main_version = (
                0 if sharc_version == "unknown" else int(sharc_version.split(".")[0])
            )
            if main_version < 4:
                logging.info(
                    "For sharc before version 4.0, we will attempt to extract charge data from QM interface settings."
                )

                qm_path = path_obj / "QM"
                for int_name, int_reader in INTERFACE_READERS.items():
                    res_dict = int_reader(iconds, qm_path)

                    if "theory_basis" in res_dict:
                        optional_settings.theory_basis_set = res_dict["theory_basis"]
                    if "est_level" in res_dict:
                        optional_settings.est_level = res_dict["est_level"]
                    if "charge" in res_dict:
                        charges = res_dict["charge"]

                    if charges is not None:
                        logging.info(f"Found charge data from the {int_name} interface")
                        break
            else:
                # Assume we are uncharged if no charge data found.
                logging.info(
                    "We assume there is no charge because no charge information was found"
                )
                charges = 0

        iconds = set_sharc_state_type_and_name_defaults(iconds, multiplicities, charges)

    iconds = assign_optional_settings(iconds, optional_settings)

    return finalize_icond_dataset(
        iconds,
        loading_parameters=loading_parameters,
        default_format_attributes=default_format_attributes,
    )


def parse_QM_in(qm_in: TextIOWrapper) -> dict[str, Any]:
    """Function to read settings of initial conditions from QM.in file.

    Will attempt to read key settings and initial positions that would usually be read form QM.log.

    Parameters
    ----------
    qm_in : TextIOWrapper
        File stream of the found `QM.in` file containing some key settings and positional information.

    Raises
    ------
    FileNotFoundError
        If parts of the file are malformed or missing.

    Returns
    -------
    dict[str, Any]
        the resulting settings in a key-value pair. Contains `num_atoms`, `num_states`, `atNames` and `atXYZ`.
    """

    info: dict[str, Any] = {}

    # Example format:
    # 12
    # Initial condition ICOND_00000/
    # C 0.00000000 0.00000000 0.00000000
    # C 2.83836864 0.00000000 0.00000000
    # C 4.30592995 2.06320298 0.00000000
    # C 3.39262531 4.73886620 0.00037795
    # H -0.72924531 -1.91996174 -0.00018897
    # H -0.75872504 0.96111471 -1.65766776
    # H -0.75872504 0.96073676 1.65785673
    # H 3.75072841 -1.84437270 0.00000000
    # H 6.35023567 1.81621578 0.00000000
    # H 1.34208349 4.84865929 0.00018897
    # H 4.09352473 5.74968071 -1.65502214
    # H 4.09333576 5.74911379 1.65596700
    # unit bohr
    # states 3
    # init
    # savedir ./SAVE/

    # H
    # DM
    # NACDR
    # GRAD

    # Get all non-empty lines
    lines = [l.strip() for l in qm_in.readlines() if len(l.strip()) > 0]
    if len(lines) < 1:
        raise FileNotFoundError("QM.in did not contain all necessary information")

    info["num_atoms"] = (num_atoms := int(lines[0]))

    if len(lines) < num_atoms + 2:
        raise FileNotFoundError("QM.in did not contain all necessary information")

    atXYZ = np.full((num_atoms, 3), np.nan)
    atNames = np.full((num_atoms), "")

    for i in range(num_atoms):
        line_parts = [x.strip() for x in lines[i + 2].split()]
        atNames[i] = line_parts[0]
        atXYZ[i] = [float(x) for x in line_parts[1:4]]

    info["atXYZ"] = atXYZ
    info["atNames"] = atNames

    for i in range(num_atoms + 2, len(lines)):
        if len(lines[i]) > 0:
            line_parts = [x.strip() for x in lines[i].split()]
            key = line_parts[0].lower()
            if len(line_parts) > 1:
                value = " ".join(line_parts[1:])
            else:
                value = True
            info[key] = value

    if "states" in info:
        state_string = info["states"]
        state_parts = state_string.split()
        info["states"] = [int(s) for s in state_parts]
        max_mult = len(state_parts)
        nsinglets = 0
        ndoublets = 0
        ntriplets = 0

        if max_mult >= 1:
            nsinglets = int(state_parts[0])
        if max_mult >= 2:
            ndoublets = int(state_parts[1])
        if max_mult >= 3:
            ntriplets = int(state_parts[2])

        # Need to be multiplied, but some may be unused? See sharc state order in documentation
        info["num_states"] = nsinglets * 1 + ndoublets * 2 + ntriplets * 3
        info["num_singlets"] = nsinglets
        info["num_doublets"] = ndoublets
        info["num_triplets"] = ntriplets
        info["max_multiplicity"] = max_mult
        if max_mult > 3:
            info["nums_higher_multiplicities"] = [int(x) for x in state_parts[3:]]
    if "charge" in info:
        info["charge"] = [int(c) for c in info["charge"].split()]

    # print(info)
    return info


def parse_QM_log(log: TextIOWrapper) -> dict[str, Any]:
    """Function to parse main information from the QM.log file

    Parameters
    ----------
    log : TextIOWrapper
        Input file wrapper to read the information from

    Raises
    ------
    ValueError
        If there are neither singlet nor triplet states listed in the QM.log file

    Returns
    -------
    dict[str, Any]
        Dictionary with key information about the system
    """
    info: dict[str, Any] = {}
    for line in log:
        if line.startswith("SHARC_version") or line.startswith("Version"):
            version_num = line.split()[1]
            info["input_format_version"] = version_num

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
            method_indicator = "Method:"
            method_parts = line[len(method_indicator)].strip().split("/")

            info["method"] = method_parts[0]
            if len(method_parts) > 1:
                info["basis"] = method_parts[1]

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

    Parameters
    ----------
    f : TextIOWrapper
        File wrapper for a `QM.log` file's contents
    out : xr.Dataset
        The dataset to write the resulting geometry to
    """

    # NB. Geometry is indeed in bohrs!
    while not (line := next(f).strip()).startswith("Geometry in Bohrs:"):
        if line.startswith("SHARC_version") or line.startswith("Version"):
            version_num = line.split()[1]
            out.attrs["input_format_version"] = version_num
        pass

    natoms = out.sizes["atom"]
    nstates = out.sizes["state"]
    tmp_names = np.full((natoms,), "")
    tmp_nums = np.full((natoms,), -1)
    tmp_positions = np.full((natoms, 3), 0.0)

    for i in range(out.sizes["atom"]):
        geometry_line = next(f).strip().split()
        atom_symbol = geometry_line[0].strip()
        tmp_names[i] = atom_symbol
        tmp_nums[i] = get_atom_number_from_symbol(atom_symbol)
        tmp_positions[i] = list(map(float, geometry_line[1:4]))

    out["atNames"][:] = tmp_names
    out["atNums"][:] = tmp_nums
    out["atXYZ"][:] = tmp_positions

    mark_variable_assigned(out.atNames)
    mark_variable_assigned(out.atNums)
    mark_variable_assigned(out.atXYZ)

    while (line := next(f).strip()).find("Final Results") < 0:
        if line.startswith("statemap:"):
            state_types = np.full((nstates,), 0)
            # state_id: [multiplicity, index_mult, ?]
            # {1: [1, 1, 0.0], 2: [1, 2, 0.0]}
            state_map = line[len("statemap:") :].strip()
            # print(state_map)

            matches = _re_state_map_entry.findall(state_map)

            if matches:
                for match in matches:
                    state_id, multiplicity, n_mult, unknown = match

                    state_id = int(state_id)
                    multiplicity = int(multiplicity)
                    n_mult = int(n_mult)
                    unknown = float(unknown)
                    state_types[state_id - 1] = multiplicity

                out.assign_coords(
                    {
                        "state_types": (
                            "state_types",
                            state_types,
                            out["state_types"].attrs,
                        )
                    }
                )

                mark_variable_assigned(out.state_types)


# Example : 12 3 ! 1 1 0 1 3 0
# Example : 12 3 ! 1 1 0
# Example : 12 3
_transition_identification_re = re.compile(
    r"(?P<num_lines>\d)\w+(?P<num_colums>\d)(\w+!\w+(?P<state_1_mult>\d)\w+(?P<state_1_n>\d)\w+(?P<state_1_misc>\d)(\w+(?P<state_2_mult>\d)\w+(?P<state_2_n>\d)\w+(?P<state_2_misc>\d))?)?"
)


def _read_dim_or_transition_identification_line(f: TextIOWrapper, main_version: int):
    """Function to optionally read the current state id or transition identification line.

    Only for versions 3.0 and up do we need to read the identification.
    Parameters
    ----------
    f : TextIOWrapper
        Input stream of text file
    main_version : int
        The main sharc version to switch the behavior


    Returns
    -------
    _type_
        _description_
    """
    # TODO: FIXME: This function could be used for debugging, currently not used.
    if main_version < 3:
        return None
    else:
        line = f.read()
        match = _transition_identification_re.match(line)
        if match:
            logging.debug(
                "Matched transition id line: %(groups)s", {'groups': match.groups()}
            )


def parse_QM_out(
    f: TextIOWrapper,
    out: xr.Dataset,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset | None:
    """Function to read all information about forces, energies, dipoles, nacs, etc. from initial condition QM.out files.

    if ``out=None`` is provided, a new Dataset is constructed and returned by the function

    Parameters
    ----------
    f : TextIOWrapper
        File input of the QM.out file to parse the data from
    out : xr.Dataset
        Target Dataset to write the loaded data into. Defaults to None.
    loading_parameters : LoadingParameters, optional
        Optional loading parameters to override variable mappings and units.

    Returns
    -------
    xr.Dataset | None
        If a new Dataset was constructed instead of being written to `out`, it will be returned.
    """
    res: xr.Dataset = out

    nstates = ConsistentValue("nstates")
    natoms = ConsistentValue("natoms")

    energy_assigned = False
    force_assigned = False
    dipole_assigned = False
    nacs_assigned = False
    phases_assigned = False

    for index, line in enumerate(f):
        line = line.strip()
        if line.startswith("SHARC_version") or line.startswith("Version"):
            version_num = line.split()[1]
            res.attrs["input_format_version"] = version_num

        if line.startswith("! 1 Hamiltonian Matrix"):
            # get number of states from dimensions of Hamiltonian
            nstates.v = int(next(f).split(" ")[0])
            energy_assigned = True

            for istate in range(nstates.v):
                energyline = re.split(" +", next(f).strip())
                res["energy"][istate] = float(energyline[2 * istate])

        elif line.startswith("! 2 Dipole Moment Matrices"):
            dim = re.split(" +", next(f).strip())
            n = int(dim[0])
            m = int(dim[1])

            dipole_assigned = True
            dip_all_tmp = nans(n, m, 3)

            dip_all_tmp[:, :, 0] = get_dipoles_per_xyz(f, n, m)
            next(f)
            dip_all_tmp[:, :, 1] = get_dipoles_per_xyz(f, n, m)
            next(f)
            dip_all_tmp[:, :, 2] = get_dipoles_per_xyz(f, n, m)

            res["dip_perm"][:], res["dip_trans"][:] = dip_sep(np.array(dip_all_tmp))

        elif line.startswith("! 3 Gradient Vectors"):
            search_res = _re_grads.search(line)
            assert search_res is not None
            get_dim = search_res.group
            nstates.v = int(get_dim("nstates"))
            natoms.v = int(get_dim("natoms"))

            force_assigned = True

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

            nacs_assigned = True

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
                phases_assigned = True
                pass

        elif line.startswith("! 8 Runtime"):
            next(f)

    if energy_assigned:
        mark_variable_assigned(res.energy)
    if force_assigned:
        mark_variable_assigned(res.forces)
    if nacs_assigned:
        mark_variable_assigned(res.nacs)
    if dipole_assigned:
        mark_variable_assigned(res.dip_perm)
        mark_variable_assigned(res.dip_trans)
    if phases_assigned:
        mark_variable_assigned(res.phases)

    # all the data has already been written to `out`
    # no need to return anything
    return None


@needs(dims={"icond"}, coords={"icond"}, not_dims={"time"})
def iconds_to_frames(iconds: xr.Dataset) -> xr.Dataset:
    """Function to convert the `icond` coordinate into a `trajid` coordinate and also build a combined `frame`+`time` multiindex as `frame`

    Will exempt atNames and atNums from the rearrangement.

    Parameters
    ----------
    iconds : xr.Dataset
        input dataset to replace the dimension in

    Raises
    ------
    ValueError
        Raised if at least one array has size 0 in one coordinate

    Returns
    -------
    xr.Dataset
        The transformed dataset
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
