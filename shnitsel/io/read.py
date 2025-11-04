import glob

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.FormatReader import FormatInformation, FormatReader
from shnitsel.io.helpers import (
    KindType,
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
import traceback
from shnitsel.io.newtonx.format_reader import NewtonXFormatReader
from shnitsel.io.pyrai2md.format_reader import PyrAI2mdFormatReader
from shnitsel.io.sharc.format_reader import SHARCFormatReader
from shnitsel.io.shnitsel.format_reader import ShnitselFormatReader
from .pyrai2md import parse_pyrai2md
from .newtonx.parse import parse_newtonx
from .sharc.parse import parse_sharc
from .shnitsel.parse import read_shnitsel_file
from tqdm.contrib.logging import logging_redirect_tqdm
from tqdm.auto import tqdm
import pandas as pd
import xarray as xr
import numpy as np
from typing import (
    Dict,
    Iterable,
    List,
    Tuple,
    TypeAlias,
    Callable,
    Literal,
    TYPE_CHECKING,
)
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import re
import logging
import os
import pathlib


# def read_trajs(
def read(
    path: PathOptionsType,
    kind: KindType | None = None,
    sub_pattern: str | None = None,
    multiple: bool = True,
    concat_method: Literal["layers", "list", "frames"] = "layers",
    parallel: bool = True,
    error_reporting: Literal["log", "raise"] = "log",
    state_names: List[str] | Callable | None = None,
    input_units: Dict[str, str] | None = None,
) -> xr.Dataset | Trajectory | List[Trajectory] | None:
    """Read all trajectories from a folder of trajectory folders

    Parameters
    ----------
    path (PathOptionsType):
        The path to the folder of folders. Can be provided as `str`, `os.PathLike` or `pathlib.Path`.
        Depending on the kind of trajectory to be loaded should denote the path of the trajectory file (``kind='shnitsel'`` or ``kind='ase'`) or a directory containing the files of the respective file format.
        Alternatively, if ``multiple=True`, this can also denote a directory containing multiple sub-directories with the actual Trajectories.
        In that case, the `concat_method` parameter should be set to specify how the .
    kind (Literal['sharc', 'nx', 'newtonx', 'pyrai2md', 'shnitsel'] | None, optional):
        The kind of trajectory, i.e. whether it was produced by SHARC, Newton-X, PyRAI2MD or Shnitsel-Tools.
        If None is provided, the function will make a best-guess effort to identify which kind of trajectory has been provided.
    sub_pattern (str|None, optional):
        If the input is a format with multiple input trajectories in different directories, this is the search pattern to append
        to the `path` (the whole thing will be read by :external:py:func:`glob.glob`).
        The default will be chosen based on `kind`, e.g., for SHARC 'TRAJ_*' or 'ICOND*' and for NewtonX 'TRAJ*'.
        If the `kind` does not support multi-folder inputs (like `shnitsel`), this will be ignored.
        If ``multiple=False``, this pattern will be ignored.
    multiple (bool, optional):
        A flag to enable loading of multiple trajectories from the subdirectories of the provided `path`.
        If set to False, only the provided path will be attempted to be loaded.
        If `sub_pattern` is provided, this parameter should not be set to `False` or the matching will be ignored.
    concat_method (Literal['layers', 'list', 'frames'])
        How to combine the loaded trajectories if multiple trajectories have been loaded.
        Defaults to ``concat_method='layers``.
        The available methods are:
        `'layers'`: Introduce a new axis `trajid` along which the different trajectories are indexed in a combined `xr.Dataset` structure.
        `'list'`: Return the multiple trajectories as a list of individually loaded data.
        `'frames'`: Concatenate the individual trajectories along the time axis ('frames') using a :external:py:class:`xarray.indexes.PandasMultiIndex`
    parallel (bool, optional):
        Whether to read multiple trajectories at the same time via parallel processing (which, in the current implementation,
        is only faster on storage that allows non-sequential reads).
        By default True.
    error_reporting (Literal['log','raise']):
        Choose whether to `log` or to `raise` errors as they occur during the import process.
        Currently, the implementation does not support `error_reporting='raise'` while `parallel=True`.
    state_names (List[str] | Callable | None, optional):
        Either a list of names to assign to states in the loaded file or a function that assigns a state name to each state id.
        If not provided or set to None, default naming will be applied, naming singlet states S0, S1,.., doublet states D0,... and triplet states T0, etc in ascending order.
    input_units: (Dict[str, str] | None, optional):
        An optional dictionary to set the units in the loaded trajectory.
        Only necessary if the units differ from that tool's default convention or if there is no default convention for the tool.
        Please refer to the names of the different unit kinds and possible values for different units in `shnitsel.units.definitions`.

    Returns
    -------
        An :external:py:class:`xarray.Dataset` containing the data of the trajectories,
        a `Trajectory` wrapper object, a list of `Trajectory` wrapper objects or `None`
        if no data could be loaded and `error_reporting='log'`.

    Raises
    ------
    FileNotFoundError
        If the `kind` does not match the provided `path` format, e.g because it does not exist or does not denote a file/directory with the required contents.
    FileNotFoundError
        If the search (``= path + pattern``) doesn't match any paths according to :external:py:func:`glob.glob`
    ValueError
        If an invalid value for ``concat_method`` is passed.
    ValueError
        If ``error_reporting`` is set to `'raise'` in combination with ``parallel=True``, the code cannot execute correctly. Only ``'log'`` is supported for parallel reading
    """
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)

    cats = {"frames": concat_trajs, "layers": layer_trajs, "list": lambda x: x}
    if concat_method not in cats:
        raise ValueError(f"`concat_method` must be one of {cats.keys()!r}")

    cat_func = cats[concat_method]

    if parallel and error_reporting != "log":
        logging.error(
            "Reading trajectories with `parallel=True` only supports `errors='log'` (the default)"
        )
        raise ValueError("parallel=True only supports errors='log' (the default)")

    loading_parameters = LoadingParameters(
        input_units=input_units,
        state_names=state_names,
        error_reporting=error_reporting,
    )

    # First check if the target path can directly be read as a Trajectory
    combined_error = None
    try:
        res = read_single(
            path, kind, error_reporting, base_loading_parameters=loading_parameters
        )

        if res is not None:
            return res

        logging.info(f"Could not read `{path}` directly as a trajectory.")
    except Exception as e:
        # Keep error in case the multiple reading also fails
        combined_error = (
            f"While trying to read as a direct trajectory: {e} [Trace:"
            + "\n".join(traceback.format_tb(e.__traceback__))
            + "]"
        )

    if multiple:
        logging.info(
            f"Attempt to read `{path}` as a directory containing multiple trajectories."
        )

        try:
            res_list = read_folder_multi(
                path,
                kind,
                sub_pattern,
                parallel,
                error_reporting,
                base_loading_parameters=loading_parameters,
            )

            if res_list is not None:
                if len(res_list) == 1:
                    return res_list[0]
                elif len(res_list) == 0:
                    message = "No trajectories could be loaded from path `{path}`."
                    if error_reporting == "log":
                        logging.error(message)
                    else:
                        raise FileNotFoundError(message)
                else:
                    return cat_func(res_list)
        except Exception as e:
            multi_error = (
                f"While trying to read as a directory containing multiple trajectories: {e} [Trace:"
                + "\n".join(traceback.format_tb(e.__traceback__))
                + "]"
            )
            combined_error = (
                multi_error
                if combined_error is None
                else combined_error + "\n" + multi_error
            )

    message = f"Could not load trajectory data from `{path}`."

    if combined_error is not None:
        message += (
            f"\nEncountered (multipe) error(s) trying to load:\n" + combined_error
        )

    if error_reporting == "log":
        logging.error(message)
        return None
    else:
        raise FileNotFoundError(message)


def read_folder_multi(
    path: PathOptionsType,
    kind: KindType | None = None,
    sub_pattern: str | None = None,
    parallel: bool = True,
    error_reporting: Literal["log", "raise"] = "log",
    base_loading_parameters: LoadingParameters | None = None,
) -> List[Trajectory] | None:
    """Function to read multiple trajectories from an input directory.

    You can either specify the kind and pattern to match relevant entries or the default pattern for `kind` will be used.
    If no `kind` is specified, all possible input formats will be checked.

    If multiple formats fit, no input will be read and either an Error will be rased or an Error will be logged and None returned.

    Otherwise, all successful reads will be returned as a list.

    Args:
        path (PathOptionsType): The path pointing to the directory where multiple trajectories may be located in the subdirectory
        kind (KindType | None,optional): The key indicating the input format.
        sub_pattern (str | None, optional): The pattern provided to "glob" to identify relevant entries in the `path` subtree. Defaults to None.
        parallel (bool, optional): A flag to enable parallel loading of trajectories. Only faster if postprocessing of read data takes up significant amounts of time. Defaults to True.
        error_reporting (Literal[&quot;log&quot;, &quot;raise&quot;], optional): Whether to raise or to log resulting errors. If errors are raised, they may also be logged. 'raise' conflicts with ``parallel=True`` setting. Defaults to "log".
        base_loading_parameters (LoadingParameters | None, optional): Base parameters to influence the loading of individual trajectories. Can be used to set default inputs and variable name mappings. Defaults to None.

    Raises:
        FileNotFoundError: If the path does not exist or Files were not founds.
        ValueError: If conflicting information of file format is detected in the target directory

    Returns:
        List[Trajectory] | None: Either a list of individual trajectories or None if loading failed.
    """

    path_obj = make_uniform_path(path)

    if not path_obj.exists() and path_obj.is_dir():
        message = f"{path} is no valid directory"
        if error_reporting == "raise":
            raise FileNotFoundError(message)
        else:
            logging.errror(message)
            return None

    relevant_kinds = [kind] if kind is not None else list(READERS.keys())

    # The kinds for which we had matches
    fitting_kinds: List[Tuple[pathlib.Path, FormatInformation]] = []
    # Entries for each kind
    matching_entries = {}

    hints_or_settings = {"kind": kind} if kind is not None else None

    for relevant_kind in relevant_kinds:
        # logging.warning(f"Considering: {relevant_kind}")
        relevant_reader = READERS[relevant_kind]

        if sub_pattern is not None:
            filter_matches = list(path_obj.glob(sub_pattern))
        else:
            filter_matches = relevant_reader.find_candidates_in_directory(path_obj)

        if filter_matches is None:
            logging.debug(f"No matches for format {relevant_kind}")
            continue

        logging.debug(
            f"Found {len(filter_matches)} matches for kind={relevant_kind}: {filter_matches}"
        )

        kind_matches = []

        kind_key = None

        for entry in filter_matches:

            # We have a match
            # logging.debug(f"Checking {entry} for format {relevant_kind}")
            try:
                res_format = relevant_reader.check_path_for_format_info(
                    entry, hints_or_settings
                )
                # res_format = identify_or_check_input_kind(entry, relevant_kind)
                if res_format is None:
                    # logging.warning(f"For {entry}, the format was None")
                    continue
                kind_key = res_format.format_name
                kind_matches.append((entry, res_format))
                # logging.info(
                #     f"Adding identified {relevant_kind}-style trajectory: {res_format}"
                # )
            except Exception as e:
                # Only consider if we hit something
                logging.debug(
                    f"Skipping {entry} for {relevant_kind} because of issue during format check: {e}"
                )
                pass

        if len(kind_matches) > 0:
            # We need to deal with the NewtonX aliases nx/newtonx
            if kind_key not in fitting_kinds:
                fitting_kinds.append(kind_key)
            matching_entries[kind_key] = kind_matches
            logging.debug(
                f"Found {len(fitting_kinds)} any appropriate matches for {relevant_kind}"
            )
        else:
            logging.debug(f"Did not find any appropriate matches for {relevant_kind}")

    if len(fitting_kinds) == 0:
        message = f"Did not detect any matching subdirectories or files for any input format in {path}"
        logging.error(message)
        if error_reporting == "raise":
            raise FileNotFoundError(message)
        else:
            return None
    elif len(fitting_kinds) > 1:
        available_formats = list(READERS.keys())
        message = f"Detected subdirectories or files of different input formats in {path} with no input format specified. Detected formats are: {fitting_kinds}. Please ensure only one format matches subdirectories in the path or denote a specific format out of {available_formats}."
        logging.error(message)
        if error_reporting == "raise":
            raise ValueError(message)
        else:
            return None
    else:
        fitting_kind = fitting_kinds[0]
        logging.debug("Opting for input format: {fitting_kind}")
        fitting_paths = matching_entries[fitting_kind]

        fitting_reader = READERS[fitting_kind]

        input_set_params = [
            (trajpath, fitting_reader, formatinfo, base_loading_parameters)
            for trajpath, formatinfo in fitting_paths
        ]
        input_paths, input_readers, input_format_info, input_loading_params = zip(
            *input_set_params
        )

        res_trajectories = []
        if parallel:
            with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
                for result in tqdm(
                    executor.map(
                        _per_traj,
                        input_paths,
                        input_readers,
                        input_format_info,
                        input_loading_params,
                    ),
                    total=len(input_set_params),
                ):
                    if result is not None and result.data is not None:
                        res_trajectories.append(result.data)
                    else:
                        logging.debug(
                            f"Reading of at least one trajectory failed. Reading routine returned value {result}."
                        )
        else:
            for params in tqdm(input_set_params, total=len(input_set_params)):
                result = _per_traj(*params)
                if result is not None and result.data is not None:
                    res_trajectories.append(result.data)
                else:
                    logging.debug(f"Failed to read trajectory from {params[1]}.")

        # TODO: FIXME: Check if trajid is actually set?
        res_trajectories.sort(
            key=lambda x: x.attrs["trajid"] if "trajid" in x.attrs else 0
        )
        return res_trajectories


def read_single(
    path: PathOptionsType,
    kind: KindType | None,
    error_reporting: Literal["log", "raise"] = "log",
    base_loading_parameters: LoadingParameters | None = None,
) -> Trajectory | None:
    try:
        res_format = identify_or_check_input_kind(path, kind)
        if res_format is not None:
            reader = READERS[res_format.format_name]
            trajectory = reader.read_from_path(
                path, res_format, base_loading_parameters
            )
            return trajectory
    except Exception as e:
        if error_reporting == "log":
            logging.exception(
                f"Caught exception while reading single trajectory input from `{path}`: \n{e}"
            )
        else:
            raise e
    return None


def identify_or_check_input_kind(
    path: PathOptionsType,
    kind_hint: KindType | None,
) -> FormatInformation | None:
    """Function to identify/guess which kind of input type the current path has if no kind was provided.
    If a kind_hint is provided, it will verify, if the path actually is of that kind

    Args:
        path (PathOptionsType): Path to a directory to be checked whether it can be read by available input readers
        kind_hint (str | None): If set, the input format specified by the user. Only that reader's result will be used eventually.

    Raises:
        FileNotFoundError: If the `path` is not valid
        ValueError: If the specified reader for `kind_hint` does not confirm validity of the directory
        ValueError: If multiple readers match and no `kind_hint` was provided.

    Returns:
        FormatInformation | None: The `FormatInformation` returned by the only successful check or None if no reader matched
    """
    # TODO: FIXME: Add ASE loading capability

    path_obj: pathlib.Path = make_uniform_path(path)

    if not path_obj.exists():
        raise FileNotFoundError("The path `{path}` is not valid.")

    # We only bother if there has been a hint to the kind of format
    # If none was specified, we take whichever fits
    is_specified_kind_satisfied = kind_hint is None
    # If the specified kind was an alias like for newtonx
    new_specified_kind = None

    resulting_format_info = {}

    hints_or_settings = {"kind": kind_hint} if kind_hint is not None else None

    for reader_kind, reader in READERS.items():
        try:
            res_format_info = reader.check_path_for_format_info(
                path_obj, hints_or_settings
            )

            if kind_hint is not None and reader_kind == kind_hint:
                is_specified_kind_satisfied = True
                new_specified_kind = res_format_info.format_name

            resulting_format_info[res_format_info.format_name] = res_format_info

        except FileNotFoundError as fn_e:
            # If required files were not found, i.e. if the path does not actually constitute input data of the denoted format
            pass
        except ValueError as v_e:
            # If the hints/settings provided by the user conflict with the requirements of the format
            pass

    if kind_hint is not None:
        if is_specified_kind_satisfied:
            return resulting_format_info[new_specified_kind]
        else:
            # The format does not fit, but another might
            message = f"The path `{path}` does not represent a directory of format `{kind_hint}`."
            possible_formats = list(resulting_format_info.keys())
            if len(possible_formats) > 0:
                joined_formats = ", ".join(possible_formats)
                message += f"\n It, however, would qualify as one of the following formats: {joined_formats}"
            else:
                message += f"\n It also didn't satisfy the conditions of any of the other known formats."

            logging.error(message)
            raise ValueError(
                f"The path `{path}` is not of the denoted format {kind_hint}."
            )
    else:
        # If there is a unique format match, use that:
        possible_formats = list(resulting_format_info.keys())
        if len(possible_formats) == 1:
            res_format = possible_formats[0]
            logging.info(
                f"Identified the path `{path}` to be of format `{res_format}`."
            )
            return resulting_format_info[res_format]
        elif len(possible_formats) > 1:
            joined_formats = ", ".join(possible_formats)
            logging.error(
                f" The path `{path}` satisfies the conditions of multiple of the known formats.: {joined_formats}. Please only provide paths containing the output data of one format."
            )
            raise ValueError(
                f"The path `{path}` is not of the denoted format {kind_hint}."
            )
        else:
            logging.error(
                f"The path `{path}` didn't satisfy the conditions of any of the known formats. Available options are: {list(READERS.keys())}"
            )

    return None


Trajid: TypeAlias = int


@dataclass
class Trajres:
    path: pathlib.Path
    misc_error: Exception | Iterable[Exception] | None
    data: Trajectory | None


# TODO: FIXME: add ASE support
_newton_reader = NewtonXFormatReader()
READERS: Dict[str, FormatReader] = {
    "nx": _newton_reader,  # parse_newtonx,
    "newtonx": _newton_reader,  # parse_newtonx,
    "sharc": SHARCFormatReader(),  # parse_sharc,
    "pyrai2md": PyrAI2mdFormatReader(),
    "shnitsel": ShnitselFormatReader(),  # read_shnitsel_file,
}


def _per_traj(
    trajdir: pathlib.Path,
    reader: FormatReader,
    format_info: FormatInformation,
    base_loading_parameters: LoadingParameters,
) -> Trajres:
    """Internal function to carry out loading of trajectories to allow for parallel processing with a ProcessExecutor.

    Args:
        trajdir (pathlib.Path): The path to read a single trajectory from
        reader (FormatReader): The reader instance to use for reading from that directory `path`.
        format_info (FormatInformation): FormatInformation obtained from previous checks of the format.
        base_loading_parameters (LoadingParameters): Settings for Loading individual trajectories like initial units and mappings of parameter names to Shnitsel variable names.

    Returns:
        Trajres|None: Either the successfully loaded trajectory in a wrapper, or the wrapper containing error information
    """

    try:
        ds = reader.read_from_path(trajdir, format_info, base_loading_parameters)
        if not ds.attrs["completed"]:
            logging.info(f"Trajectory at path {trajdir} did not complete")

        return Trajres(path=trajdir, misc_error=None, data=ds)

    except Exception as err:
        # This is fairly common and will be reported at the end
        logging.info(
            f"Reading of trajectory from path {trajdir} failed:\n"
            + str(err)
            + f"\nSkipping {trajdir}."
        )

        return Trajres(
            path=trajdir,
            misc_error=[err],
            data=None,
        )


def gather_traj_metadata(datasets: Iterable[Trajectory], time_dim="time") -> np.ndarray:
    """Function to gather metadate from a set of trajectories.

    Used to combine trajectories into one aggregate Dataset.

    Args:
        datasets (Iterable[Trajectory]): The sequence of trajctories for which metadata should be collected
        time_dim (str, optional): The name of the time dimension in the input datasets. Defaults to "time".

    Returns:
        np.ndarray: The resulting meta information
    """
    # TODO: FIXME: Rewrite such that result is a dict of conflicting settings and another one of parameters in agreement.
    # Only conflicting settings need to be indexed, others can be applied to full trajectory instead.
    # TODO: Potentially also collect atom number and other information that needs to match to be combined

    num_datasets = len(datasets)
    traj_meta = {
        "trajid": np.full((num_datasets,), -1, dtype="i4"),
        "delta_t": np.full((num_datasets,), np.nan, dtype="f8"),
        "max_ts": np.full((num_datasets,), -1, dtype="i4"),
        "t_max": np.full((num_datasets,), np.nan, dtype="f8"),
        "completed": np.full((num_datasets,), False, dtype="?"),
        "nsteps": np.full((num_datasets,), -1, dtype="i4"),
    }

    # TODO: FIXME: Check for consistency of more of the units and attributes
    for i, ds in enumerate(datasets):
        traj_meta["trajid"][i] = ds.attrs.get("trajid", -1)
        traj_meta["delta_t"][i] = ds.attrs.get("delta_t", np.nan)
        traj_meta["max_ts"][i] = ds.attrs.get("max_ts", -1)
        traj_meta["t_max"][i] = ds.attrs.get("t_max", np.nan)
        # TODO: FIXME: think about whether or not to default to False for completed parameter
        traj_meta["completed"][i] = ds.attrs.get("completed", False)
        traj_meta["nsteps"][i] = len(ds.indexes[time_dim])

    return traj_meta


def concat_trajs(datasets: Iterable[Trajectory]) -> Trajectory:
    """Function to concatenate multiple trajectories along their `time` dimension.

    Will create one continuous time dimension like an extended trajectory

    Args:
        datasets (Iterable[Trajectory]): Datasets representing the individual trajectories

    Raises:
        ValueError: Raised if there is conflicting input meta data.
        ValueError: Raised if there are no trajectories provided to this function.

    Returns:
        Trajectory: The combined and extended trajectory with a new leading `frame` dimension
    """

    # TODO:FIXME:
    datasets = list(datasets)

    if len(datasets) == 0:
        raise ValueError("No trajectories were provided.")

    if all("time" in ds.coords for ds in datasets):
        # We ensure time is always called time when loading
        time_dim = "time"
    else:
        raise ValueError("Some trajectories do not have a 'time' coordinate.")

    # TODO: Deal with trajid not being set yet
    datasets = [
        ds.expand_dims(trajid=[ds.attrs["trajid"]]).stack(frame=["trajid", time_dim])
        for ds in datasets
    ]

    # TODO: FIXME: Deal with issues arising from inconsisten meta information. E.g. ensure same number of atoms, consistent units, etc.
    frames = xr.concat(datasets, dim="frame", combine_attrs="drop_conflicts")
    traj_meta = gather_traj_metadata(datasets, time_dim=time_dim)
    frames = frames.assign_coords(trajid_=traj_meta["trajid"])
    # TODO: FIXME: Consider the naming convention of trajid and trajid_ being somewhat confusing
    frames = frames.assign(
        delta_t=("trajid_", traj_meta["delta_t"]),
        max_ts=("trajid_", traj_meta["max_ts"]),
        completed=("trajid_", traj_meta["completed"]),
        nsteps=("trajid_", traj_meta["nsteps"]),
    )

    frames = Trajectory(frames)

    if TYPE_CHECKING:
        assert isinstance(frames, Trajectory)

    return Trajectory(frames)


def layer_trajs(datasets: Iterable[Trajectory]) -> Trajectory:
    """Function to combine trajctories into one Dataset by creating a new dimension 'trajid' and indexing the different trajectories along that.

    Will create one new trajid dimension.

    Args:
        datasets (Iterable[xr.Dataset]): Datasets representing the individual trajectories

    Raises:
        ValueError: Raised if there is conflicting input meta data.
        ValueError: Raised if there are no trajectories provided to this function.


    Returns:
        xr.Dataset: The combined and extended trajectory with a new leading `trajid` dimension
    """

    meta = gather_traj_metadata(datasets)

    trajids = meta["trajid"]

    datasets = [ds.expand_dims(trajid=[id]) for ds, id in zip(datasets, trajids)]

    # trajids = pd.Index(meta["trajid"], name="trajid")
    # coords_trajids = xr.Coordinates(indexes={"trajid": trajids})
    # breakpoint()
    # TODO: FIXME: Deal with issues arising from inconsisten meta information. E.g. ensure same number of atoms, consistent units, etc.
    layers = xr.concat(datasets, dim="trajid", combine_attrs="drop_conflicts")

    layers = layers.assign_coords(trajid=trajids)

    # TODO: FIXME: All units should be converted to same unit
    # TODO: FIXME: All inconsistent meta data/attr should be stored into a meta_data object or lead to an error

    # del meta["trajid"]
    layers = layers.assign(
        {k: xr.DataArray(v, dims=["trajid"]) for k, v in meta.items() if k != "trajid"}
    )
    if TYPE_CHECKING:
        assert isinstance(layers, xr.Dataset)
    return Trajectory(layers)
