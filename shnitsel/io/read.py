import glob

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.FormatReader import FormatInformation, FormatReader
from shnitsel.io.helpers import KindType, LoadingParameters, PathOptionsType, make_uniform_path
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
from typing import Dict, Iterable, List, TypeAlias, Callable, Literal, TYPE_CHECKING
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
    concat_method: Literal['layers', 'list', 'frames'] = 'layers',
    parallel: bool = True,
    error_reporting: Literal['log', 'raise'] = 'log',
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

    cats = {'frames': concat_trajs, 'layers': layer_trajs, 'list': lambda x: x}
    if concat_method not in cats:
        raise ValueError(f"`concat_method` must be one of {cats.keys()!r}")

    cat_func = cats[concat_method]

    if parallel and error_reporting != 'log':
        logging.error(
            "Reading trajectories with `parallel=True` only supports `errors='log'` (the default)")
        raise ValueError(
            "parallel=True only supports errors='log' (the default)")

    loading_parameters = LoadingParameters(
        input_units=input_units,
        state_names=state_names,
        error_reporting=error_reporting)

    # First check if the target path can directly be read as a Trajectory
    combined_error = None
    try:
        res = read_single(path, kind, error_reporting,
                          base_loading_parameters=loading_parameters)

        if res is not None:
            return res

        logging.info(f"Could not read `{path}` directly as a trajectory.")
    except Exception as e:
        # Keep error in case the multiple reading also fails
        combined_error = f"While trying to read as a direct trajectory: {e}"

    if multiple:
        logging.info(
            f"Attempt to read `{path}` as a directory containing multiple trajectories.")

        try:
            res_list = read_folder_multi(
                path, kind, sub_pattern,
                parallel, error_reporting,
                base_loading_parameters=loading_parameters)

            if res_list is not None:
                if len(res_list) == 1:
                    return res_list[0]
                elif len(res_list) == 0:
                    message = "No trajectories could be loaded from path `{path}`."
                    if error_reporting == 'log':
                        logging.error(message)
                    else:
                        raise FileNotFoundError(message)
                else:
                    return cat_func(res_list)
        except Exception as e:
            multi_error = f"While trying to read as a directory containing multiple trajectories: {e}"
            combined_error = multi_error if combined_error is None else combined_error+"\n" + multi_error

    message = f"Could not load trajectory data from `{path}`."

    if combined_error is not None:
        message += f"\nEncountered multipe errors trying to load:\n"+combined_error

    if error_reporting == 'log':
        logging.error(message)
        return None
    else:
        raise FileNotFoundError(message)


def read_folder_multi(
        path: PathOptionsType,
        kind: KindType | None,
        sub_pattern: str | None = None,
        parallel: bool = True,
        error_reporting: Literal['log', 'raise'] = 'log',
        base_loading_parameters: LoadingParameters | None = None) -> List[Trajectory] | None:

    # TODO: Go over readers, find matched trajectories, deal with multiple matches, read all matches and return.
    # TODO: Read in parallel if setting allows it.
    # TODO: Make individual reader routines use the loading parameters to construct the Dataset
    if parallel:
        datasets = read_trajs_parallel(paths, kind)
    else:
        datasets = read_trajs_list(
            paths, kind, error_reporting=error_reporting)
    pass


def read_single(
        path: PathOptionsType,
        kind: KindType | None,
        error_reporting: Literal['log', 'raise'] = 'log',
        base_loading_parameters: LoadingParameters | None = None) -> Trajectory | None:
    try:
        res_format = identify_or_check_input_kind(path, kind)
        if res_format is not None:
            reader = READERS[res_format.format_name]
            trajectory = reader.read_from_path(path, res_format)
            return trajectory
    except Exception as e:
        if error_reporting == 'log':
            logging.exception(
                f"Caught exception while reading single trajectory input from `{path}`: \n{e}")
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

    hints_or_settings = {'kind': kind_hint} if kind_hint is not None else None

    for reader_kind, reader in READERS.items():
        try:
            res_format_info = reader.check_path_for_format_info(
                path_obj, hints_or_settings)

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
                f"The path `{path}` is not of the denoted format {kind_hint}.")
    else:
        # If there is a unique format match, use that:
        possible_formats = list(resulting_format_info.keys())
        if len(possible_formats) == 1:
            res_format = possible_formats[0]
            logging.info(
                f"Identified the path `{path}` to be of format `{res_format}`.")
            return resulting_format_info[res_format]
        elif len(possible_formats) > 1:
            joined_formats = ", ".join(possible_formats)
            logging.error(
                f" The path `{path}` satisfies the conditions of multiple of the known formats.: {joined_formats}. Please only provide paths containing the output data of one format.")
            raise ValueError(
                f"The path `{path}` is not of the denoted format {kind_hint}.")
        else:
            logging.error(
                f"The path `{path}` didn't satisfy the conditions of any of the known formats. Available options are: {list(READERS.keys())}")

    return None


Trajid: TypeAlias = int

_exnum = re.compile("[0-9]+")


@dataclass
class Trajres:
    trajid: int
    missing_file: str | None
    misc_error: Exception | None
    data: xr.Dataset | None


def _default_idfn(path):
    global _exnum
    res = _exnum.search(os.path.basename(path))
    if res is None:
        raise ValueError(f"Could not extract trajid from path '{path}'")

    return int(res[0])


_idfn = _default_idfn
_read_traj: Callable

# TODO: FIXME: add ASE support
_newton_reader = NewtonXFormatReader()
READERS: Dict[str, FormatReader] = {
    'nx': _newton_reader,  # parse_newtonx,
    'newtonx': _newton_reader,  # parse_newtonx,
    'sharc': SHARCFormatReader(),  # parse_sharc,
    'pyrai2md': PyrAI2mdFormatReader(),
    'shnitsel': ShnitselFormatReader(),  # read_shnitsel_file,
}


def read_trajs_list(paths, kind, idfn=None, sort=True, error_reporting='log'):
    global _default_idfn
    if idfn is None:
        idfn = _default_idfn

    if kind not in READERS:
        raise ValueError(
            f"'kind' should be one of {list(READERS)}, rather than '{kind}'"
        )
    read_traj = READERS[kind]

    if error_reporting not in {'log', 'raise'}:
        raise ValueError(
            f"'errors' should be one of ['log', 'raise'], rather than '{error_reporting}'"
        )

    datasets = []
    with logging_redirect_tqdm():
        missing_files: dict[str, list[Trajid]] = {}
        misc_errors: dict[Trajid, Exception] = {}
        incomplete: list[Trajid] = []
        for trajdir in tqdm(paths):
            trajid = idfn(trajdir)
            try:
                ds = read_traj(trajdir)
            except FileNotFoundError as err:
                # This is fairly common and will be reported at the end
                logging.info(
                    f"Missing file for trajectory {trajid} at path {trajdir}:\n"
                    + str(err)
                    + f"\nSkipping {trajid}."
                )
                missing = os.path.basename(err.filename)
                missing_files[missing] = missing_files.get(
                    missing, []) + [trajid]

                continue

            except Exception as err:
                if error_reporting == 'log':
                    # Miscellaneous exceptions could indicate a problem with the parser
                    # so they enjoy a more imposing loglevel
                    logging.error(
                        f"Error for trajectory {trajid} at path {trajdir}:\n"
                        + str(err)
                        + f"\nSkipping {trajid}."
                    )
                    misc_errors[trajid] = err
                    continue
                elif error_reporting == 'raise':
                    raise err

            if not ds.attrs['completed']:
                logging.info(
                    f"Trajectory {trajid} at path {trajdir} did not complete")
                incomplete.append(trajid)

            ds.attrs['trajid'] = trajid

            datasets.append(ds)

    if sort:
        datasets.sort(key=lambda x: x.attrs['trajid'])

    if len(misc_errors):
        print("Miscellaneous errors:")
        for trajid, merr in misc_errors.items():
            print(f"{trajid:>6}  {merr}")
    if len(missing_files):
        for fname, trajids in missing_files.items():
            trajids.sort()
            print(
                f"Skipped {len(trajids)} trajectories missing file '{fname}', IDs:",
                ' '.join([str(t) for t in trajids]),
            )

    if len(incomplete):
        incomplete.sort()
        print(
            f"Included {len(incomplete)} incomplete trajectories, IDs:",
            ' '.join([str(i) for i in incomplete]),
        )

    return datasets


def _per_traj(trajdir):
    trajid = _idfn(trajdir)
    missing_file = None
    misc_error = None

    try:
        ds = _read_traj(trajdir)

    except FileNotFoundError as err:
        # This is fairly common and will be reported at the end
        logging.info(
            f"Missing file for trajectory {trajid} at path {trajdir}:\n"
            + str(err)
            + f"\nSkipping {trajid}."
        )
        missing_file = os.path.basename(err.filename)

        return Trajres(
            trajid=trajid,
            missing_file=missing_file,
            misc_error=misc_error,
            data=None,
        )

    except Exception as err:
        # Miscellaneous exceptions could indicate a problem with the parser
        # so they enjoy a more imposing loglevel
        logging.error(
            f"Error for trajectory {trajid} at path {trajdir}:\n"
            + str(err)
            + f"\nSkipping {trajid}."
        )
        misc_error = err

        return Trajres(
            trajid=trajid,
            missing_file=missing_file,
            misc_error=misc_error,
            data=None,
        )

    if not ds.attrs['completed']:
        logging.info(f"Trajectory {trajid} at path {trajdir} did not complete")

    ds.attrs['trajid'] = trajid

    return Trajres(
        trajid=trajid, missing_file=missing_file, misc_error=misc_error, data=ds
    )


def read_trajs_parallel(
    paths: Iterable[str | os.PathLike],
    kind: KindType,
    idfn=None,
    sort=True,
):
    global _idfn
    global _read_traj

    if idfn is None:
        _idfn = _default_idfn
    else:
        _idfn = idfn

    if kind not in READERS:
        raise ValueError(
            f"'kind' should be one of {list(READERS)}, rather than '{kind}'"
        )
    _read_traj = READERS[kind]

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        res = list(tqdm(executor.map(_per_traj, paths), total=len(paths)))

    datasets = []
    missing_files: dict[str, list[Trajid]] = {}
    misc_errors: dict[str, list[Trajid]] = {}
    incomplete: list[Trajid] = []
    for t in res:
        if (ds := t.data) is not None:
            datasets.append(ds)
            if not ds.attrs['completed']:
                incomplete.append(t.trajid)

        if (mf := t.missing_file) is not None:
            if mf not in missing_files:
                missing_files[mf] = []
            missing_files[mf].append(t.trajid)

        if (me := t.misc_error) is not None:
            sme = str(me)
            if sme not in missing_files:
                missing_files[sme] = []
            missing_files[sme].append(t.trajid)

    if sort:
        datasets.sort(key=lambda x: x.attrs['trajid'])

    if len(misc_errors):
        print("Miscellaneous errors:")
        for trajid, err in misc_errors.items():
            print(f"{trajid:>6}  {err}")
    if len(missing_files):
        for fname, trajids in missing_files.items():
            trajids.sort()
            print(
                f"Skipped {len(trajids)} trajectories missing file '{fname}', IDs:",
                ' '.join([str(t) for t in trajids]),
            )

    if len(incomplete):
        incomplete.sort()
        print(
            f"Included {len(incomplete)} incomplete trajectories, IDs:",
            ' '.join([str(i) for i in incomplete]),
        )

    return datasets


def gather_traj_metadata(datasets: Iterable[xr.Dataset], time_dim='ts') -> np.ndarray:
    traj_meta = np.zeros(
        len(datasets),
        dtype=[
            ('trajid', 'i4'),
            ('delta_t', 'f8'),
            ('max_ts', 'i4'),
            ('completed', '?'),
            ('nsteps', 'i4'),
        ],
    )

    for i, ds in enumerate(datasets):
        traj_meta['trajid'][i] = ds.attrs.get('trajid', -1)
        traj_meta['delta_t'][i] = ds.attrs.get('delta_t', np.nan)
        traj_meta['max_ts'][i] = ds.attrs.get('max_ts', -1)
        traj_meta['completed'][i] = ds.attrs['completed']
        traj_meta['nsteps'][i] = len(ds.indexes[time_dim])

    return traj_meta


def concat_trajs(datasets: Iterable[xr.Dataset]) -> xr.Dataset:
    """Function to concatenate multiple trajectories along their `time` dimension.

    Will create one continuous time dimension like an extended trajectory

    Args:
        datasets (Iterable[xr.Dataset]): Datasets representing the individual trajectories

    Raises:
        ValueError: Raised if there are conflicting `time` or `ts` dimensions
        ValueError: Raised if there are no trajectories provided to this function.

    Returns:
        xr.Dataset: The combined and extended trajectory with a new leading `frame` dimension
    """

    datasets = list(datasets)

    if len(datasets) == 0:
        raise ValueError("No trajectories were provided.")

    if all('ts' in ds.coords for ds in datasets):
        time_dim = 'ts'
    elif all('time' in ds.coords for ds in datasets):
        time_dim = 'time'
    else:
        raise ValueError(
            "Some trajectories have coordinate 'ts', others 'time'. "
            "Please resolve this inconsistency manually."
        )

    datasets = [
        ds.expand_dims(trajid=[ds.attrs['trajid']]).stack(
            frame=['trajid', time_dim])
        for ds in datasets
    ]

    frames = xr.concat(datasets, dim='frame', combine_attrs='drop_conflicts')
    traj_meta = gather_traj_metadata(datasets, time_dim=time_dim)
    frames = frames.assign_coords(trajid_=traj_meta['trajid'])
    frames = frames.assign(
        delta_t=('trajid_', traj_meta['delta_t']),
        max_ts=('trajid_', traj_meta['max_ts']),
        completed=('trajid_', traj_meta['completed']),
        nsteps=('trajid_', traj_meta['nsteps']),
    )
    if TYPE_CHECKING:
        assert isinstance(frames, xr.Dataset)
    return frames


def layer_trajs(datasets: Iterable[xr.Dataset]) -> xr.Dataset:
    meta = gather_traj_metadata(datasets)

    trajids = pd.Index(meta['trajid'], name='trajid')
    coords_trajids = xr.Coordinates(indexes={'trajid': trajids})
    breakpoint()
    layers = xr.concat(datasets, dim=trajids, combine_attrs='drop_conflicts')

    del meta['trajid']
    layers = layers.assign(
        {k: xr.DataArray(v, coords_trajids) for k, v in meta.items()}
    )
    if TYPE_CHECKING:
        assert isinstance(layers, xr.Dataset)
    return layers
