# pyright: basic
import itertools
import os
import logging
import math

from collections.abc import Iterable

import xarray as xr
import numpy as np
import pandas as pd

def midx_combs(values: pd.core.indexes.base.Index|list, name: str|None =None):
    if name is None:
        try:
            # if `values` is a `pandas.core.indexes.base.Index`
            # extract its name
            name = values.name
        except AttributeError:
            raise ValueError("need to specify name if values lack name attribute")

    return xr.Coordinates.from_pandas_multiindex(
      pd.MultiIndex.from_tuples(
        itertools.combinations(values, 2),
        names=['from', 'to']
      ),
      dim=f'{name}comb'
    )


def flatten_midx(ds, idx_name, renamer=None):
    midx = ds.indexes[idx_name]
    to_drop = midx.names + [midx.name]
    fidx = midx.to_flat_index()
    if renamer is not None:
        fidx = [renamer(x, y) for x,y in fidx]
    return ds.drop_vars(to_drop).assign_coords({idx_name: fidx})

def expand_midx(ds, midx_name, level_name, value):
    midx = ds.indexes[midx_name]
    to_drop = [midx.name] + midx.names
    df = midx.to_frame()
    df.insert(0, level_name, [value]*len(midx)) # in place!
    midx = pd.MultiIndex.from_frame(df)
    coords = xr.Coordinates.from_pandas_multiindex(midx, dim=midx_name)
    return (
      ds
      .drop_vars(to_drop)
      .assign_coords(coords)
    )

# TODO! Improve this, then use it to rewrite ts_to_time.
# Only worth doing once you find a second situation in which
# MultiIndex is to be modified.
def modify_midx_level(data, midx_name, coords, replace: set | None):
    """Turns the levels of a MultiIndex into coordinates,
    modifies them as coordinates, then converts back.
    """
    # find MultiIndex name
    # for idx_name, idx in data.indexes.values():
        # if not set(coords).isdisjoint(getattr(idx, 'names', set())):

    # 'frame' in energy_rel.reset_index('frame').coords['ts'].dims

    level_names = list(data.indexes[midx_name].names)
    data = (data
      .reset_index(midx_name)
      .assign_coords(coords)
    )
    if replace is not None:
        new_levels = set(level_names).union(coords.keys()).difference(replace)
        data = (data
          .reset_index(midx_name)
          .set_xindex(list(new_levels))
        )
    return data

def _ts_to_time(data, delta_t=None, swap_index=True):
    default_delta_t = 0.5

    if delta_t is None:
        delta_t = data.attrs.get('delta_t', default_delta_t)

    return modify_midx_level(
      data,
      midx_name='frame',
      coords={'time': data.coords['ts'] * delta_t},
      replace={'ts'} if swap_index else None
    )

def open_frames(path):
    try:
        frames = xr.open_dataset(path)
    except ValueError as err:
        if not os.path.exists(path):
            raise FileNotFoundError
        else:
            raise err

    tcoord = None
    if 'time' in frames.coords:
        tcoord = 'time'
    elif 'ts' in frames.coords:
        tcoord = 'ts'

    if tcoord is not None:
        frames = frames.set_xindex(['trajid', tcoord])

    return frames.set_xindex(['from', 'to'])


def save_frames(frames, path, complevel=9):
    encoding = {
        var: {"compression": "gzip", "compression_opts": complevel} for var in frames
    }

    # TODO generalize?
    # netcdf does not support bool
    if 'completed' in frames.data_vars:
        frames = frames.assign(completed=frames.completed.astype('i1'))
    if 'completed' in frames.attrs:
        frames.attrs['completed'] = int(frames.attrs['completed'])
    # or MultiIndex
    frames.reset_index(['frame', 'statecomb']).to_netcdf(
        path, engine='h5netcdf', encoding=encoding
    )

def split_for_saving(frames, bytes_per_chunk=50e6):
    trajids = frames.get('trajid_', np.unique(frames['trajid']))
    ntrajs = len(trajids)
    nchunks = math.trunc(frames.nbytes / 50e6)
    logging.debug(f"{nchunks=}")
    indices = np.trunc(np.linspace(0, ntrajs, nchunks + 1)).astype(np.integer)
    logging.debug(f"{indices=}")
    trajidsets = [trajids[a:z].values for a, z in zip(indices[:-1], indices[1:])]
    logging.debug(f"{trajidsets=}")
    return [sel_trajs(frames, trajids[a:z]) for a, z in zip(indices[:-1], indices[1:])]


def save_split(
    frames, path_template, bytes_per_chunk=50e6, complevel=9, ignore_errors=False
):
    dss = split_for_saving(frames, bytes_per_chunk=bytes_per_chunk)
    for i, ds in enumerate(dss):
        current_path = path_template.format(i)
        try:
            save_frames(ds, current_path, complevel=complevel)
        except Exception as e:
            logging.error(f"Exception while saving to {current_path=}")
            if not ignore_errors:
                raise e


#######################################
# Functions to extend xarray selection:


def sel_trajs(
    frames: xr.Dataset, trajids_or_mask: int | Iterable[bool | int], invert=False
) -> xr.Dataset:
    """Takes boolean mask (each entry corresponding to a trajectory in MuliIndex order)
    or an Iterable of trajids. In latter case, will not generally return trajectories in order given.
    `invert=True` selects those trajectories not specified."""
    trajids_or_mask = np.atleast_1d(trajids_or_mask)
    if np.issubdtype(trajids_or_mask.dtype, np.integer):
        trajids = trajids_or_mask
    elif np.issubdtype(trajids_or_mask.dtype, bool):
        mask = trajids_or_mask
        if 'trajid_' in frames.dims:
            trajids = frames['trajid_'][mask]
        else:
            raise NotImplementedError(
                "Indexing trajids with a boolean mask is only supported when the "
                "coordinate 'trajid_' is present"
            )
    else:
        raise TypeError(
            "Only indexing using a boolean mask or integer trajectory IDs is supported; "
            f"the detected dtype was {trajids_or_mask.dtype}"
        )
    return sel_trajids(frames=frames, trajids=trajids, invert=invert)


def sel_trajids(
    frames: xr.Dataset, trajids: int | Iterable[int], invert=False
) -> xr.Dataset:
    "Will not generally return trajectories in order given"
    trajids = np.atleast_1d(trajids)
    # check that all trajids are valid, as Dataset.sel() would
    if not invert and not (np.isin(trajids, frames['trajid'])).all():
        missing = trajids[~np.isin(trajids, frames['trajid'])]
        raise KeyError(
            f"Of the supplied trajectory IDs, {len(missing)} were "
            f"not found in index 'trajid': {missing}"
        )
    mask = frames['trajid'].isin(trajids)
    if invert:
        mask = ~mask
    res = frames.sel(frame=mask)

    if 'trajid_' in frames.dims:
        actually_selected = np.unique(res['trajid'])
        res = res.sel(trajid_=actually_selected)
    return res


# def selframemask(frames, frame_mask):
#     return frames.sel(frame=frame_mask)