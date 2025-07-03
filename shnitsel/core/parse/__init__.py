import glob
import os
import logging
import re
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor
from typing import TypeAlias, Callable

import numpy as np
import xarray as xr
import pandas as pd

from tqdm.auto import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from . import sharc_icond as sharc_icond
from . import pyrai2md as pyrai2md
from . import nx, sharc_traj

Trajid: TypeAlias = int

_exnum = re.compile('[0-9]+')

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
    'nx': nx.read_traj,
    'sharc': sharc_traj.read_traj,
    'pyrai2md': pyrai2md.read_traj,
}


def read_trajs_list(paths, kind, idfn=None, sort=True):
    global _default_idfn
    if idfn is None:
        idfn = _default_idfn

    if kind not in READERS:
        raise ValueError(
            f"'kind' should be one of {list(READERS)}, rather than '{kind}'"
        )
    read_traj = READERS[kind]

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
                missing_files[missing] = missing_files.get(missing, []) + [trajid]

                continue
            except Exception as err:
                # Miscellaneous exceptions could indicate a problem with the parser
                # so they enjoy a more imposing loglevel
                logging.error(
                    f"Error for trajectory {trajid} at path {trajdir}:\n"
                    + str(err)
                    + f"\nSkipping {trajid}."
                )
                misc_errors[trajid] = err
                continue

            if not ds.attrs['completed']:
                logging.info(f"Trajectory {trajid} at path {trajdir} did not complete")
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


def read_trajs_parallel(paths, kind, idfn=None, sort=True):
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
        res = list(executor.map(_per_traj, paths))

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


def gather_traj_metadata(datasets, time_dim='ts'):
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
        traj_meta['trajid'][i] = ds.attrs['trajid']
        traj_meta['delta_t'][i] = ds.attrs['delta_t']
        traj_meta['max_ts'][i] = ds.attrs['max_ts']
        traj_meta['completed'][i] = ds.attrs['completed']
        traj_meta['nsteps'][i] = len(ds.indexes[time_dim])

    return traj_meta


def concat_trajs(datasets):
    if all('ts' in ds.coords for ds in datasets):
        time_dim = 'ts'
    elif all('time' in ds.coords for ds in datasets):
        time_dim = 'time'
    else:
        ValueError(
            "Some trajectories have coordinate 'ts', others 'time'. "
            "Please resolve this inconsistency manually."
        )

    datasets = [
        ds.expand_dims(trajid=[ds.attrs['trajid']]).stack(frame=['trajid', time_dim])
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
    return frames


def layer_trajs(datasets):
    meta = gather_traj_metadata(datasets)

    trajids = pd.Index(meta['trajid'], name='trajid')
    coords_trajids = xr.Coordinates(indexes={'trajid': trajids})
    breakpoint()
    layers = xr.concat(datasets, dim=trajids, combine_attrs='drop_conflicts')

    del meta['trajid']
    layers = layers.assign(
        {k: xr.DataArray(v, coords_trajids) for k, v in meta.items()}
    )
    return layers


def read_trajs(path, kind, glob_suffix='TRAJ*', format='frames', parallel=False):
    glob_expr = os.path.join(path, glob_suffix)
    paths = glob.glob(glob_expr)
    if len(paths) == 0:
        raise FileNotFoundError(
            f"The search '{glob_expr}' didn't match any paths "
            f"under path '{path}' "
            f"given working directory '{os.getcwd()}'"
        )
    if parallel:
        datasets = read_trajs_parallel(paths, kind)
    else:
        datasets = read_trajs_list(paths, kind)

    cats = {'frames': concat_trajs, 'layers': layer_trajs}

    try:
        cat_func = cats[format]
    except KeyError:
        raise ValueError(f"`format` must be one of {cats.keys()!r}")

    return cat_func(datasets)