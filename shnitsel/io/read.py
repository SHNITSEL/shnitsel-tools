
import glob
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


@dataclass
class LoadingParameters:
    # The path to either an input file or an input trajectory depending on the kind of trajectory being requested
    input_path: str | os.PathLike
    # An indicator as to which kind of trajectory is being loaded
    kind: Literal['sharc', 'newtonx', 'pyrai2md', 'shnitsel'] | None

    # A dict containing the information, which input observable has which unit. If not provided, the loader will guess the units either based on the default values of that simulator or the data in `path`
    input_units: Dict[str, str] | None = None
    # List of the names of states or a function to label them or None and let the trajectory loader make an educated guess
    state_names: List[str] | Callable | None = None

    # Parameter to control whether multiple trajectories will be concatenated into a continuous trajectory or layered as trajectories indexed by
    concat_function: Callable = 'frames'
    # Flag to indicate whether parallel loading is requested
    parallel: bool = True
    # Flag to set how errors during loading are reported
    error_reporting: Literal['log', 'raise'] = 'log'

    # Pattern for matching paths within `input_path` for certain loaders depending on `kind`
    sub_pattern: str | None = 'TRAJ*'


# def read_trajs(
def read(
    path: str | os.PathLike,
    kind: Literal['sharc', 'newtonx', 'pyrai2md', 'shnitsel'] | None,
    sub_pattern: str | None = 'TRAJ*',
    concat_method: Literal['frames', 'layers'] = 'frames',
    parallel: bool = True,
    error_reporting: Literal['log', 'raise'] = 'log',
) -> xr.Dataset:
    """Read all trajectories from a folder of trajectory folders

    Parameters
    ----------
    path
        The path to the folder of folders
    kind
        The kind of trajectory, i.e. whether it was produced by SHARC, Newton-X, PyRAI2MD or Shnitsel-Tools.
        If None is provided, the function will make a best-guess effort to identify, which kind of trajectory has been provided
    sub_pattern
        If the input is a format with multiple input trajectories in different directories, this is the search pattern to append 
        to the `path` (the whole thing will be read by :external:py:func:`glob.glob`).
        By default 'TRAJ*'.
        If the `kind` does not support multi-folder inputs (like `shnitsel`), this will be ignored
    concat_method
        Whether to return the trajectories concatenated along the time axis ('frames') using a
        :external:py:class:`xarray.indexes.PandasMultiIndex`
        or along a new axis ('layers'), by default 'frames'
    parallel
        Whether to read multiple trajectories at the same time via parallel processing (which, in the current implementation,
        is only faster on storage that allows non-sequential reads).
        By default True.
    error_reporting:
        Choose whether to log or to raise errors as they occur during the import process. 
        Currently, the implementation does not support `error_reporting='raise'` while `parallel=True`.

    Returns
    -------
        An :external:py:class:`xarray.Dataset` containing the data of the trajectories

    Raises
    ------
    FileNotFoundError
        If the `kind` of input requires a single file and `path` does not exist or does not denote a file.
    FileNotFoundError
        If the search (``= path + pattern``) doesn't match any paths according to :external:py:func:`glob.glob`
    ValueError
        If an invalid value for ``concat_method`` is passed.
    ValueError
        If ``error_reporting`` is set to `'raise'` in combination with ``parallel=True``, the code cannot execute correctly. Only ``'log'`` is supported for parallel reading
    """
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)

    cats = {'frames': concat_trajs, 'layers': layer_trajs}
    if concat_method not in cats:
        raise ValueError(f"`concat_method` must be one of {cats.keys()!r}")

    cat_func = cats[concat_method]

    if parallel and error_reporting != 'log':
        logging.error(
            "Reading trajectories with `parallel=True` only supports `errors='log'` (the default)")
        raise ValueError(
            "parallel=True only supports errors='log' (the default)")

    glob_expr = os.path.join(path, sub_pattern)
    paths = glob.glob(glob_expr, root_dir=path)
    if len(paths) == 0:
        msg = f"The search '{glob_expr}' didn't match any paths"
        if not os.path.isabs(path):
            msg += f", relative to working directory '{os.getcwd()}'"
        raise FileNotFoundError(msg)

    if parallel:
        datasets = read_trajs_parallel(paths, kind)
    else:
        datasets = read_trajs_list(
            paths, kind, error_reporting=error_reporting)

    return cat_func(datasets)


def guess_input_kind(path: str | os.PathLike, sub_pattern: str | None) -> Literal['sharc', 'newtonx', 'pyrai2md', 'shnitsel', 'ase'] | None:
    # TODO: FIXME: Add ASE loading capability

    path_obj = pathlib.Path(path)

    if not path_obj.exists():
        return None

    if path_obj.is_file():
        # Only shnitsel and ase are file-based
        if path_obj.parts[-1].endswith(".nc"):
            return 'shnitsel'
        elif path_obj.parts[-1].endswith(".db"):
            return 'ase'
        else:
            return None
    elif path_obj.is_dir():
        possible_formats = ['sharc', 'newtonx', 'pyrai2md']

        


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

READERS = {
    'nx': parse_newtonx,
    'sharc': parse_sharc,
    'pyrai2md': parse_pyrai2md,
    'shnitsel': read_shnitsel_file
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


def read_trajs_parallel(paths: Iterable[str | os.PathLike], kind: Literal['sharc', 'newtonx', 'pyrai2md', 'shnitsel'], idfn=None, sort=True):
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

    if len(datasets) == 0:
        raise ValueError("No trajectories were parsed successfully.")

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
