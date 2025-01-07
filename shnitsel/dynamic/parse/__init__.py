import glob
import os
import logging
import re

import numpy as np
import xarray as xr
import pandas as pd

from . import nx, sharc_traj

_exnum = re.compile('[0-9]+')


def _default_idfn(path):
    global _exnum
    return int(os.path.basename(_exnum.search(path)[0]))


def read_traj_frames(path, kind, idfn=None):
    global _default_idfn
    if idfn is None:
        idfn = _default_idfn
    trajid = idfn(path)

    if kind == 'nx':
        read_traj = nx.read_traj
    elif kind == 'sharc':
        read_traj = sharc_traj.read_traj

    try:
        ds = read_traj(path)
    except:
        logging.error(f"Error for trajectory {trajid} at path {path}")
        raise

    if not ds.attrs['completed']:
        logging.warn(f"Trajectory {trajid} at path {path} did not complete")

    ds.attrs['trajid'] = trajid

    return ds


def read_trajs_list(path, kind):
    from tqdm import tqdm
    from tqdm.contrib.logging import logging_redirect_tqdm

    paths = glob.glob(os.path.join(path, 'TRAJ*'))
    assert len(paths) > 0

    with logging_redirect_tqdm():
        datasets = [read_traj_frames(p, kind=kind) for p in tqdm(paths)]

    return datasets


def gather_traj_metadata(datasets):
    traj_meta = np.zeros(
        len(datasets), dtype=[('trajid', 'i4'), ('delta_t', 'f8'), ('nsteps', 'i4')]
    )

    for i, ds in enumerate(datasets):
        traj_meta['trajid'][i] = ds.attrs['trajid']
        traj_meta['delta_t'][i] = ds.attrs['delta_t']
        traj_meta['nsteps'][i] = len(ds.indexes['ts'])

    return traj_meta


def concat_trajs(datasets):
    datasets = [
        ds.expand_dims(trajid=[ds.attrs['trajid']]).stack(frame=['trajid', 'ts'])
        for ds in datasets
    ]

    frames = xr.concat(datasets, dim='frame', combine_attrs='drop_conflicts')
    traj_meta = gather_traj_metadata(datasets)
    frames = frames.assign_coords(trajid_=traj_meta['trajid'])
    frames = frames.assign(
        delta_t=('trajid_', traj_meta['delta_t']),
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


def read_trajs(path, kind, format='frames'):
    datasets = read_trajs_list(path, kind)

    cats = {'frames': concat_trajs, 'layers': layer_trajs}

    try:
        cat_func = cats[format]
    except KeyError:
        raise ValueError(f"`format` must be one of {cats.keys()!r}")

    return cat_func(datasets)