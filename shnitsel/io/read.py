from types import UnionType
from typing_extensions import TypeForm
from shnitsel.core._api_info import API, internal
from shnitsel.core.typedefs import StateTypeSpecifier
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.traj_combiner_methods import (
    concat_trajs,
    db_from_data,
    layer_trajs,
)
from shnitsel.data.dataset_containers import Trajectory, Frames

from shnitsel.data.tree import (
    ShnitselDB,
    CompoundGroup,
    DataGroup,
    DataLeaf,
)
from shnitsel.data.tree.node import TreeNode
from shnitsel.io.format_reader_base import FormatInformation, FormatReader
from shnitsel.io.format_registry import get_available_io_handlers, FormatIdentifierType
from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
import traceback
from shnitsel.io.shared.messages import (
    collect_and_clean_queue_handler,
    handle_records,
    setup_queue_handler,
)
from tqdm.contrib.logging import logging_redirect_tqdm
from tqdm.auto import tqdm
import pandas as pd
import xarray as xr
import numpy as np
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Callable,
    Literal,
    TYPE_CHECKING,
    Sequence,
    TypeVar,
    overload,
)
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import logging
import os
import pathlib

from shnitsel.data.xr_io_compatibility import SupportsFromXrConversion

DataType = TypeVar("DataType")


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["db"] = "db",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: type[DataType] | UnionType = Trajectory,
# ) -> DataType | TreeNode[Any, DataType] | None: ...


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["db"] = "db",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: None = None,
# ) -> TreeNode | None: ...


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["list"] = "list",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: type[DataType] | UnionType = Trajectory,
# ) -> DataType | List[DataType | TreeNode[Any, DataType]] | None: ...


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["list"] = "list",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: None = None,
# ) -> (
#     ShnitselDataset
#     | SupportsFromXrConversion
#     | xr.Dataset
#     | xr.DataArray
#     | TreeNode
#     | List[
#         ShnitselDataset
#         | SupportsFromXrConversion
#         | xr.Dataset
#         | xr.DataArray
#         | TreeNode
#     ]
#     | None
# ): ...


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["layers"] = "layers",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: type[DataType] | UnionType = Trajectory,
# ) -> DataType | None: ...


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["layers"] = "layers",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: None = None,
# ) -> ShnitselDataset | SupportsFromXrConversion | xr.Dataset | xr.DataArray | None: ...


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["frames"] = "frames",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: type[DataType] | UnionType = Trajectory,
# ) -> DataType | Frames | None: ...


# @overload
# def read(
#     path: PathOptionsType,
#     kind: FormatIdentifierType | None = None,
#     sub_pattern: str | None = None,
#     multiple: bool = True,
#     concat_method: Literal["frames",] = "frames",
#     parallel: bool = True,
#     error_reporting: Literal["log", "raise"] = "log",
#     input_units: Dict[str, str] | None = None,
#     input_state_types: StateTypeSpecifier
#     | List[StateTypeSpecifier]
#     | Callable[[xr.Dataset], xr.Dataset]
#     | None = None,
#     input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
#     input_trajectory_id_maps: Dict[str, int]
#     | Callable[[pathlib.Path], int]
#     | None = None,
#     expect_dtype: None = None,
# ) -> Frames | SupportsFromXrConversion | xr.Dataset | xr.DataArray | None: ...


# def read_trajs(
@API()
def read(
    path: PathOptionsType,
    kind: FormatIdentifierType | None = None,
    *,
    sub_pattern: str | None = None,
    multiple: bool = True,
    concat_method: Literal[
        "db",
        "layers",
        "list",
        "frames",
    ] = "db",
    parallel: bool = True,
    error_reporting: Literal["log", "raise"] = "log",
    input_units: Dict[str, str] | None = None,
    input_state_types: StateTypeSpecifier
    | List[StateTypeSpecifier]
    | Callable[[xr.Dataset], xr.Dataset]
    | None = None,
    input_state_names: List[str] | Callable[[xr.Dataset], xr.Dataset] | None = None,
    input_trajectory_id_maps: Dict[str, int]
    | Callable[[pathlib.Path], int]
    | None = None,
    expect_dtype: type[DataType] | UnionType | None = None,
) -> (
    xr.Dataset
    | xr.DataArray
    | ShnitselDataset
    | SupportsFromXrConversion
    | TreeNode[
        Any,
        ShnitselDataset | SupportsFromXrConversion | xr.Dataset | xr.DataArray,
    ]
    | TreeNode[Any, DataType]
    | Sequence[xr.Dataset | ShnitselDataset | SupportsFromXrConversion | xr.DataArray]
    | DataType
):
    """Read all trajectories from a folder of trajectory folder.

    The function will attempt to automatically detect the type of the trajectory if `kind` is not set.
    If `path` is a directory containing multiple trajectory sub-directories or files with `multiple=True`, this function will attempt to load all those subdirectories in parallel.
    To limit the number of considered trajectories, you can provide `sub_pattern` as a glob pattern to filter directory entries to be considered
    It will extract as much information from the trajectory as possible and return it in a standard shnitsel format.

    If multiple trajectories are loaded, they need to be combined into one return object. The method for this can be configured via `concat_method`.
    By default, `concat_method='layers'`, a new dimension `trajid` will be introduced and different trajectories can be identified by their index along this dimension.
        Please note, that additional entries along the `time` dimension in any variable will be padded by default values.
        You can either check the `max_ts` attribute for the maximum time index in the respective directory or check whether there are `np.nan` values in any of the observables.
        We recommend using the energy variable.
    `concat_method='frames'` introduces a new dimension `frame` where each tick is a combination of `trajid` and `time` in the respective trajectory. Therefore, only valid frames will be present and no padding performed.
    `concat_method='list'` simply returns the list of successfully loaded trajectories without merging them.
    `concat_method='db'` returns a Tree-structured ShnitselDB object containing all of the trajectories. Only works if all trajectories contain the same compound/molecule.
    For concatenation except `'list'`, the same number of atoms and states must be present in all individual trajectories.

    Error reporting can be configure between logging or raising exceptions via `error_reporting`.

    If `parallel=True`, multiple processes will be used to load multiple different trajectories in parallel.

    As some formats do not contain sufficient information to extract the input units of all variables, you can provide units (see `shnitsel.units.definitions.py` for unit names) of individual variables via `input_units`.
    `input_units` should be a dict mapping default variable names to the respective unit.
    The individual variable names should adhere to the shnitsel-format standard, e.g. atXYZ, force, energy, dip_perm. Unknown names or names not present in the loaded data will be ignored without warning.
    If no overrides are provided, the read function will use internal defaults for all variables.

    Similarly, as many output formats do not provide state multiplicity or state name information, we allow for the provision of state types (via `input_state_types`)
    and of state names (via `input_state_names`).
    Both can either be provided as a list of values for the states in the input in ascending index order or as a function that assigns the correct values to the coordinates `state_types` or `state_names` in the trajectory respectively.
    Types are either `1`, `2`, or `3`, whereas names are commonly of the format "S0", "D0", "T0".
    Do not modify any other variables within the respective function.
    If you modify any variable, use the `mark_variable_assigned(variable)` function, i.e. `mark_variable_assigned(dataset.state_types)` or `mark_variable_assigned(dataset.state_names)` respectively, to notify shnitsel of the respective update.
    If the notification is not applied, the coordinate may be dropped due to a supposed lack of assigned values.

    If multiple trajectories are merged, it is importand to be able to distinguish which one may be referring.
    By setting `input_trajectory_id_maps`, you can provide a mapping between input paths and the id you would like to assign to the trajectory read from that individual path as a dict.
    The key should be the absolute path as a posix-conforming string.
    The value should be the desired id. Note that ids should be pairwise distinct.
    Alternatively, `input_trajectory_id_maps` can be a function that is provided the `pathlib.Path` object of the trajectory input path and should return an associated id.
    By default, ids are exctracted from integers in the directory names of directory-based inputs.
    If no integer is found or the format does not support the directory-style input, a random id will be assigned by default.

    Parameters
    ----------
    path : PathOptionsType
        The path to the folder of folders. Can be provided as `str`, `os.PathLike` or `pathlib.Path`.
        Depending on the kind of trajectory to be loaded should denote the path of the trajectory file (``kind='shnitsel'`` or ``kind='ase'`) or a directory containing the files of the respective file format.
        Alternatively, if ``multiple=True`, this can also denote a directory containing multiple sub-directories with the actual Trajectories.
        In that case, the `concat_method` parameter should be set to specify how the .
    kind : FormatIdentifierType, optional
        The kind of trajectory, i.e. whether it was produced by SHARC, Newton-X, PyRAI2MD or Shnitsel-Tools.
        If None is provided, the function will make a best-guess effort to identify which kind of trajectory has been provided.
    sub_pattern : str, optional
        If the input is a format with multiple input trajectories in different directories, this is the search pattern to append
        to the `path` (the whole thing will be read by :external:py:func:`glob.glob`).
        The default will be chosen based on `kind`, e.g., for SHARC 'TRAJ_*' or 'ICOND*' and for NewtonX 'TRAJ*'.
        If the `kind` does not support multi-folder inputs (like `shnitsel`), this will be ignored.
        If ``multiple=False``, this pattern will be ignored.
    multiple: bool, optional
        A flag to enable loading of multiple trajectories from the subdirectories of the provided `path`.
        If set to False, only the provided path will be attempted to be loaded.
        If `sub_pattern` is provided, this parameter should not be set to `False` or the matching will be ignored.
    concat_method : Literal['db', 'layers', 'list', 'frames']
        How to combine the loaded trajectories if multiple trajectories have been loaded.
        Defaults to ``concat_method='db'``.
        The available methods are:
        `'db'` : Returns the trajectories/data points in a hierarchical tree structure to allow for easier management of complex data hierarchies.
        `'layers'`: Introduce a new axis `trajid` along which the different trajectories are indexed in a combined `xr.Dataset` structure.
        `'list'`: Return the multiple trajectories as a list of individually loaded data.
        `'frames'`: Concatenate the individual trajectories along the time axis ('frames') using a :external:py:class:`xarray.indexes.PandasMultiIndex`
    parallel : bool, optional
        Whether to read multiple trajectories at the same time via parallel processing (which, in the current implementation,
        is only faster on storage that allows non-sequential reads).
        By default True.
    error_reporting : Literal['log','raise'], optional
        Choose whether to `log` or to `raise` errors as they occur during the import process.
        Currently, the implementation does not support `error_reporting='raise'` while `parallel=True`.
    input_units : dict[str, str], optional
        An optional dictionary to set the units in the loaded trajectory.
        Only necessary if the units differ from that tool's default convention or if there is no default convention for the tool.
        Please refer to the names of the different unit kinds and possible values for different units in `shnitsel.units.definitions`.
    input_state_types : list[int] | Callable[[xr.Dataset], xr.Dataset], optional
        Either a list of state types/multiplicities to assign to states in the loaded trajectories or a function that assigns a state multiplicity to each state.
        The function may use all of the information in the trajectory if required and should return the updated Dataset.
        If not provided or set to None, default types/multipliciteis will be applied based on extracted numbers of singlets, doublets and triplets. The first num_singlet types will be set to `1`, then 2*num_doublet types will be set to `2` and then 3*num_triplets types will be set to 3.
        Will be invoked/applied before the `input_state_names` setting.
    input_state_names : list[str] | Callable[[xr.Dataset], xr.Dataset], optional
        Either a list of names to assign to states in the loaded file or a function that assigns a state name to each state.
        The function may use all of the information in the trajectory, i.e. the state_types array, and should return the updated Dataset.
        If not provided or set to None, default naming will be applied, naming singlet states S0, S1,.., doublet states D0,... and triplet states T0, etc in ascending order.
        Will be invoked/applied after the `input_state_types` setting.
    input_trajectory_id_maps : dict[str, int]| Callable[[pathlib.Path], int], optional
        A dict mapping absolut posix paths to ids to be applied or a function to convert a path into an integer id to assign to the trajectory.
        If not provided, will be chosen either based on the last integer matched from the path or at random up to `2**31-1`.
    expected_dtype: type[DataType] | UnionType, optional
        An explicit type hint to control the output type of this function where template arguments are concerned.
        Will be explicitly set on `ShnitselDB` nodes.
        If not provided, may be inferred internally.

    Returns
    -------
    Trajectory | Frames | DataType | xr.Dataset | xr.DataArray
        For simple inputs like single trajectories or non-hierarchical inputs, this function will
        return trajectory data or the data stored in the file that was attempted to be read.
        If `concat_method='frames'`, multiple data entries will be combined into a single `MultiFrames` object.
    List[Trajectory | Frames] | List[DataType]
        If `concat_method='list'` and multiple data entries were read, a list of that data may be returned.
    ShnitselDB[Trajectory | Frames]
    | ShnitselDB[DataType]
    | CompoundGroup[DataType]
    | DataGroup[DataType]
    | DataLeaf[DataType]
        If a file with hierarchical data was read or if `concat_method='db'` was set,
        a hierarchical structure will be returned.
        If such a structure was constructed, it will always complete the tree up to the `ShnitselDB` root.
        If a tree structure is read from file, completion is not automatically performed.
    xr.Dataset
        If no conversion is possible, the data is most likely returned as an xr.Dataset or a list thereof.
    None
        If no data could be loaded.

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

    cats = {
        "frames": concat_trajs,
        "layers": layer_trajs,
        "db": db_from_data,
        "list": lambda x: x,
    }
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
        error_reporting=error_reporting,
        trajectory_id=input_trajectory_id_maps,
        state_types=input_state_types,
        state_names=input_state_names,
    )

    # First check if the target path can directly be read as a Trajectory
    combined_error = None
    try:
        res = read_single(
            path,
            kind,
            error_reporting,
            base_loading_parameters=loading_parameters,
            expect_dtype=expect_dtype,
        )

        if res is not None:
            return res

        logging.debug(f"Could not read `{path}` directly as a trajectory.")
    except Exception as e:
        # Keep error in case the multiple reading also fails
        combined_error = (
            f"While trying to read as a direct trajectory: {e} [Trace:"
            + "\n".join(traceback.format_tb(e.__traceback__))
            + "]"
        )

    if multiple:
        logging.debug(
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
                expect_dtype=expect_dtype,
            )

            if res_list is not None:
                if len(res_list) == 0:
                    message = f"No trajectories could be loaded from path `{path}`."
                    if error_reporting == "log":
                        logging.error(message)
                    else:
                        raise FileNotFoundError(message)
                else:
                    return cat_func(res_list, dtype=expect_dtype)
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
        message += "\nEncountered (multipe) error(s) trying to load:\n" + combined_error

    if error_reporting == "log":
        logging.error(message)
        raise FileNotFoundError(message)
        return None
    else:
        raise FileNotFoundError(message)


@internal()
def read_folder_multi(
    path: PathOptionsType,
    kind: FormatIdentifierType | None = None,
    sub_pattern: str | None = None,
    parallel: bool = True,
    error_reporting: Literal["log", "raise"] = "log",
    base_loading_parameters: LoadingParameters | None = None,
    expect_dtype: type[DataType] | UnionType | None = None,
) -> (
    Sequence[
        xr.Dataset
        | xr.DataArray
        | ShnitselDataset
        | SupportsFromXrConversion
        | TreeNode[
            Any, ShnitselDataset | SupportsFromXrConversion | xr.Dataset | xr.DataArray
        ]
        | TreeNode[Any, DataType]
        | Sequence[
            xr.Dataset | ShnitselDataset | SupportsFromXrConversion | xr.DataArray
        ]
        | DataType
    ]
    | None
):
    """
    Function to read multiple trajectories from an input directory.

    You can either specify the kind and pattern to match relevant entries or the default pattern for `kind` will be used.
    If no `kind` is specified, all possible input formats will be checked.

    If multiple formats fit, no input will be read and either an Error will be rased or an Error will be logged and None returned.

    Otherwise, all successful reads will be returned as a list.

    Parameters
    ----------
    path : PathOptionsType, optional
        The path pointing to the directory where multiple trajectories may be located in the subdirectory, by default None,
    kind : FormatIdentifierType, optional
        The key indicating the input format, will be inferred if not provided.
    sub_pattern : str, optional
        The pattern provided to "glob" to identify relevant entries in the `path` subtree. Defaults to None.
    parallel : bool, optional
        A flag to enable parallel loading of trajectories. Only faster if postprocessing of read data takes up significant amounts of time. Defaults to True.
    error_reporting : Literal["log", "raise"], optional
        Whether to raise or to log resulting errors. If errors are raised, they may also be logged. 'raise' conflicts with ``parallel=True`` setting. Defaults to "log".
    base_loading_parameters : LoadingParameters, optional
        Base parameters to influence the loading of individual trajectories. Can be used to set default inputs and variable name mappings. Defaults to None.
    expect_dtype : type[DataType]  |  TypeForm[DataType], optional
        An explicit type hint to control the output type of this function where template arguments are concerned.
        Will be explicitly set on `ShnitselDB` nodes.
        If not provided, may be inferred internally.

    Returns
    -------
    list[Trajectory] | list[...] None
        Either a list of individual trajectories, a list of various possible result types read from file or None if loading failed.

    Raises
    ------
    FileNotFoundError
        If the path does not exist or Files were not founds.
    ValueError
        If conflicting information of file format is detected in the target directory
    """

    path_obj = make_uniform_path(path)

    if not path_obj.exists() and path_obj.is_dir():
        message = f"{path} is no valid directory"
        if error_reporting == "raise":
            raise FileNotFoundError(message)
        else:
            logging.error(message)
            return None

    if base_loading_parameters is None:
        base_loading_parameters = LoadingParameters()

    READERS = get_available_io_handlers()

    relevant_kinds = [kind] if kind is not None else list(READERS.keys())

    # The kinds for which we had matches
    fitting_kinds: List[str] = []
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
            if kind_key is not None and kind_key not in fitting_kinds:
                fitting_kinds.append(kind_key)
            matching_entries[kind_key] = kind_matches
            logging.debug(
                f"Found {len(fitting_kinds)} any appropriate matches for {relevant_kind}"
            )
        else:
            logging.debug(f"Did not find any appropriate matches for {relevant_kind}")

    if len(fitting_kinds) == 0:
        raise_msg: str
        if kind is not None:
            message = f"Did not detect any matching subdirectories or files for input format %s in `%s`."
            logging.error(message, kind, path)
            raise_msg = message % (kind, path)
        else:
            message = f"Did not detect any matching subdirectories or files for any supported input format in `%s`."
            logging.error(message, path)
            raise_msg = message % path

        if error_reporting == "raise":
            raise FileNotFoundError(raise_msg)
        else:
            return None
    elif len(fitting_kinds) > 1:
        available_formats = list(READERS.keys())
        message = (
            f"Detected subdirectories or files of different input formats in {path} with no input format specified. \n"
            f"Detected formats are: {fitting_kinds}. \n"
            f"Please ensure only one format matches subdirectories in the path or denote a specific format out of {available_formats}. \n"
            "You can import data from different import formats via multiple calls to `read()` and subsequent merging e.g. with `tree_merge()`.\n"
            "Specify your desired input format with the `kind` parameter of `read()`."
        )
        logging.error(message)
        if error_reporting == "raise":
            raise ValueError(message)
        else:
            return None
    else:
        fitting_kind = fitting_kinds[0]
        logging.debug(f"Opting for input format: {fitting_kind}")
        fitting_paths = matching_entries[fitting_kind]

        fitting_reader = READERS[fitting_kind]

        input_set_params = [
            (
                trajpath,
                fitting_reader,
                formatinfo,
                base_loading_parameters,
                expect_dtype,
            )
            for trajpath, formatinfo in fitting_paths
        ]
        (
            input_paths,
            input_readers,
            input_format_info,
            input_loading_params,
            expected_type,
        ) = zip(*input_set_params)
        log_messages: list[logging.LogRecord] = []
        res_trajectories: Sequence[
            xr.Dataset
            | xr.DataArray
            | ShnitselDataset
            | SupportsFromXrConversion
            | TreeNode[
                Any,
                ShnitselDataset | SupportsFromXrConversion | xr.Dataset | xr.DataArray,
            ]
            | TreeNode[Any, DataType]
            | Sequence[
                xr.Dataset | ShnitselDataset | SupportsFromXrConversion | xr.DataArray
            ]
            | DataType
        ] = []
        if parallel:
            with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
                for result in tqdm(
                    executor.map(
                        _per_traj,
                        input_paths,
                        input_readers,
                        input_format_info,
                        input_loading_params,
                        expected_type,
                    ),
                    total=len(input_set_params),
                ):
                    if result is not None:
                        if result.data is not None:
                            res_trajectories.append(result.data)
                        if result.log_records is not None:
                            log_messages += result.log_records
                    else:
                        logging.debug(
                            f"Reading of at least one trajectory failed. Reading routine returned value {result}."
                        )
        else:
            for params in tqdm(input_set_params, total=len(input_set_params)):
                result = _per_traj(*params)
                if result is not None:
                    if result.data is not None:
                        res_trajectories.append(result.data)
                    if result.log_records is not None:
                        log_messages += result.log_records
                else:
                    logging.debug(f"Failed to read trajectory from {params[1]}.")

        # Output collected logging messages from child processes
        handle_records(log_messages, None)

        # TODO: FIXME: Check if trajid is actually set?
        res_trajectories.sort(
            key=lambda x: x.attrs.get("trajid", x.attrs.get("trajectory_id", 0))
            if isinstance(x, xr.Dataset | ShnitselDataset)
            and ("trajid" in x.attrs or "trajectory_id" in x.attrs)
            else np.random.randint(0, high=100000)
        )
        return res_trajectories


@internal()
def read_single(
    path: PathOptionsType,
    kind: FormatIdentifierType | None,
    error_reporting: Literal["log", "raise"] = "log",
    base_loading_parameters: LoadingParameters | None = None,
    expect_dtype: type[DataType] | UnionType | None = None,
) -> (
    xr.Dataset
    | xr.DataArray
    | ShnitselDataset
    | SupportsFromXrConversion
    | TreeNode[
        Any, ShnitselDataset | SupportsFromXrConversion | xr.Dataset | xr.DataArray
    ]
    | TreeNode[Any, DataType]
    | Sequence[xr.Dataset | ShnitselDataset | SupportsFromXrConversion | xr.DataArray]
    | DataType
    | None
):
    """Helper function to read input from a single input path.

    May yield complex and iterable data structures depending on the input format.

    Parameters
    ----------
    path : PathOptionsType
        Path to a directory to be checked whether it can be read by available input readers
    kind_hint : str | None
        If set, the input format specified by the user. Only that reader's result will be used eventually.

    Raises
    ------
    FileNotFoundError
        If the `path` is not valid
    TypeError
        If the `expected_dtype` does not match the data parsed from file.

    Returns
    -------
    xr.Dataset
    | xr.DataArray
    | ShnitselDataset
    | SupportsFromXrConversion
    | TreeNode[
        Any, ShnitselDataset | SupportsFromXrConversion | xr.Dataset | xr.DataArray
    ]
    | TreeNode[Any, DataType]
    | Sequence[xr.Dataset | ShnitselDataset | SupportsFromXrConversion]
    | DataType
    | None
        The data that has been read from the `path` location.
    """
    queue, handler, logger, original_handlers = setup_queue_handler(None, 'root')

    if base_loading_parameters is None:
        base_loading_parameters = LoadingParameters()
        base_loading_parameters.error_reporting = error_reporting
        base_loading_parameters.logger = logger

    try:
        res_format = identify_or_check_input_kind(path, kind)
        if res_format is not None:
            READERS = get_available_io_handlers()
            reader = READERS[res_format.format_name.lower()]
            # TODO: FIXME: Rename to more general read_data()?
            trajectory = reader.read_data(
                path, res_format, base_loading_parameters, expect_dtype=expect_dtype
            )
            # TODO: FIXME: Deal with a full SchnitselDB being loaded from a single file in a directory and then combined with others.
            return trajectory
    except Exception as e:
        if error_reporting == "log":
            logging.exception(
                f"Caught exception while reading single trajectory input from `{path}`: \n{e}"
            )
        else:
            records = collect_and_clean_queue_handler(
                queue,
                handler,
                logger,
                original_handlers=original_handlers,
                doCollect=True,
            )
            if records is not None:
                handle_records(records, None)
            raise
    records = collect_and_clean_queue_handler(
        queue, handler, logger, original_handlers=original_handlers, doCollect=True
    )
    if records is not None:
        handle_records(records, None)
    return None


@internal()
def identify_or_check_input_kind(
    path: PathOptionsType,
    kind_hint: FormatIdentifierType | None,
) -> FormatInformation | None:
    """Function to identify/guess which kind of input type the current path has if no kind was provided.
    If a kind_hint is provided, it will verify, if the path actually is of that kind

    Parameters
    ----------
    path : PathOptionsType
        Path to a directory to be checked whether it can be read by available input readers
    kind_hint : str | None
        If set, the input format specified by the user. Only that reader's result will be used eventually.

    Raises
    ------
    FileNotFoundError
        If the `path` is not valid
    ValueError
        If the specified reader for `kind_hint` does not confirm validity of the directory
    ValueError
        If multiple readers match and no `kind_hint` was provided.

    Returns
    -------
    FormatInformation | None
        The `FormatInformation` returned by the only successful check or None if no reader matched
    """
    # TODO: FIXME: Add ASE loading capability

    path_obj: pathlib.Path = make_uniform_path(path)  # type: ignore # will always yield a pathlib.Path

    if not path_obj.exists():
        raise FileNotFoundError(f"The path `{path}` is not valid.")

    # We only bother if there has been a hint to the kind of format
    # If none was specified, we take whichever fits
    is_specified_kind_satisfied = kind_hint is None
    # If the specified kind was an alias like for newtonx
    new_specified_kind = None

    resulting_format_info = {}

    hints_or_settings = {"kind": kind_hint} if kind_hint is not None else None

    READERS = get_available_io_handlers()

    for reader_kind, reader in READERS.items():
        try:
            res_format_info = reader.check_path_for_format_info(
                path_obj, hints_or_settings
            )

            if kind_hint is not None and reader_kind == kind_hint:
                is_specified_kind_satisfied = True
                new_specified_kind = res_format_info.format_name.lower()

            resulting_format_info[res_format_info.format_name.lower()] = res_format_info

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
            message = f"The path `{path}` does not represent a directory of requested format `{kind_hint}`."
            possible_formats = list(resulting_format_info.keys())
            if len(possible_formats) > 0:
                joined_formats = ", ".join(possible_formats)
                message += f"\n It, however, would qualify as one of the following formats: {joined_formats}"
            else:
                message += "\n It also didn't satisfy the conditions of any of the other known formats."

            logging.info(message)
            # raise ValueError(
            #     f"The path `{path}` is not of the denoted format {kind_hint}."
            # )
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
            logging.warning(
                f" The path `{path}` satisfies the conditions of multiple of the known formats.: {joined_formats}. \n Please only provide paths containing the output data of one format or specify the desired output format."
            )
            # raise ValueError(
            #     f"The path `{path}` is not of the denoted format {kind_hint}."
            # )
        else:
            logging.info(
                f"The path `{path}` didn't satisfy the conditions of any of the known formats. Available options are: {list(READERS.keys())} but none matched the specific path."
            )

    return None


@dataclass
class Trajres:
    path: pathlib.Path
    misc_error: tuple[Exception, Any] | Iterable[tuple[Exception, Any]] | None
    data: (
        xr.Dataset
        | xr.DataArray
        | ShnitselDataset
        | SupportsFromXrConversion
        | TreeNode
        | Sequence[
            xr.Dataset | ShnitselDataset | SupportsFromXrConversion | xr.DataArray
        ]
        | None
    )
    log_records: list[logging.LogRecord] | None


def _per_traj(
    trajdir: pathlib.Path,
    reader: FormatReader,
    format_info: FormatInformation,
    base_loading_parameters: LoadingParameters,
    expect_dtype: type | UnionType | None = None,
) -> Trajres:
    """Internal function to carry out loading of trajectories to allow for parallel processing with a ProcessExecutor.

    Parameters
    ----------
    trajdir : pathlib.Path
        The path to read a single trajectory from
    reader : FormatReader
        The reader instance to use for reading from that directory `path`.
    format_info : FormatInformation
        FormatInformation obtained from previous checks of the format.
    base_loading_parameters : LoadingParameters
        Settings for Loading individual trajectories like initial units and mappings of parameter names to Shnitsel variable names.

    Returns
    -------
    Trajres|None
        Either the successfully loaded trajectory in a wrapper, or the wrapper containing error information
    """
    queue, handler, logger, original_handlers = setup_queue_handler(None, 'root')
    base_loading_parameters.logger = logger

    try:
        ds = reader.read_data(
            trajdir,
            format_info=format_info,
            loading_parameters=base_loading_parameters,
            expect_dtype=expect_dtype,
        )

        if isinstance(ds, (xr.Dataset, Trajectory, Frames)):
            if ds is not None and not ds.attrs.get("completed", False):
                logging.info(f"Trajectory at path {trajdir} did not complete")

        records = collect_and_clean_queue_handler(
            queue, handler, logger, original_handlers=original_handlers, doCollect=True
        )

        return Trajres(path=trajdir, misc_error=None, data=ds, log_records=records)

    except Exception as err:
        # This is fairly common and will be reported at the end
        logging.exception(
            f"Reading of trajectory from path {trajdir} failed:\n"
            + str(err)
            + f"Trace:{traceback.format_exc()}"
            + f"\nSkipping {trajdir}."
        )

        records = collect_and_clean_queue_handler(
            queue, handler, logger, original_handlers=original_handlers, doCollect=True
        )
        return Trajres(
            path=trajdir,
            misc_error=[(err, traceback.format_exc())],
            data=None,
            log_records=records,
        )
