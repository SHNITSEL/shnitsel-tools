from __future__ import annotations

import itertools
import os
import logging
import math
import typing
from typing import Sequence

import numpy.typing as npt

from xarray.core.groupby import DatasetGroupBy, DataArrayGroupBy

import xarray as xr
import numpy as np
import pandas as pd

from .._contracts import needs


def replace_total(
    da: xr.DataArray, to_replace: np.ndarray | list, value: np.ndarray | list
):
    """Replaces each occurence of `to_replace` in `da` with the corresponding element of `value`.
    Replacement must be total, i.e. every element of `da` must be in `to_replace`.
    This permits a change of dtype between `to_replace` and `value`.
    This function is based on the snippets at https://github.com/pydata/xarray/issues/6377

    Parameters
    ----------
    da
        An xr.DataArray
    to_replace
        Values to replace
    value
        Values with which to replace them

    Returns
    -------
        An xr.DataArray with dtype matching `value`.
    """
    to_replace = np.array(to_replace)
    value = np.array(value)
    flat = da.values.ravel()

    sorter = np.argsort(to_replace)
    insertion = np.searchsorted(to_replace, flat, sorter=sorter)
    indices = np.take(sorter, insertion, mode='clip')
    replaceable = to_replace[indices] == flat

    out = value[indices[replaceable]]
    return da.copy(data=out.reshape(da.shape))


def midx_combs(values: pd.core.indexes.base.Index|list, name: str|None =None):
    if name is None:
        if hasattr(values, 'name'):
            # if `values` is a `pandas.core.indexes.base.Index`
            # extract its name
            name = values.name
        else:
            raise ValueError("need to specify name if values lack name attribute")

    return xr.Coordinates.from_pandas_multiindex(
      pd.MultiIndex.from_tuples(
        itertools.combinations(values, 2),
        names=['from', 'to']
      ),
      dim=f'{name}comb'
    )


def flatten_midx(
    obj: xr.Dataset | xr.DataArray, idx_name: str, renamer: callable | None = None
) -> xr.Dataset | xr.DataArray:
    midx = obj.indexes[idx_name]
    to_drop = midx.names + [midx.name]
    fidx = midx.to_flat_index()
    if renamer is not None:
        fidx = [renamer(x, y) for x,y in fidx]
    return obj.drop_vars(to_drop).assign_coords({idx_name: fidx})


def flatten_levels(
    obj: xr.Dataset | xr.DataArray,
    idx_name: str,
    levels: Sequence[str],
    new_name: str | None = None,
    position: int = 0,
    renamer: typing.Callable | None = None,
) -> xr.Dataset | xr.DataArray:
    dims = obj.coords[idx_name].dims
    if len(dims) != 1:
        raise ValueError(
            f"Expected index '{idx_name}' to be associated with one dimension, "
            f"but it is associated with {len(dims)} dimensions: {dims}."
        )
    dim = dims[0]
    old = obj.indexes[idx_name]
    if new_name is None:
        new_name = levels[-1]
    df = old.to_frame().drop(columns=levels)

    # Construct flat index with only the specified levels:
    for level in old.names:
        if level not in levels:
            old = old.droplevel(level)
    fidx = old.to_flat_index()

    if renamer is not None:
        fidx = [renamer(x, y) for x, y in fidx]
    df.insert(position, new_name, fidx)
    new = pd.MultiIndex.from_frame(df)
    return obj.drop_vars(idx_name).assign_coords({idx_name: (dim, new)})


def expand_midx(
    obj: xr.Dataset | xr.DataArray, midx_name, level_name, value
) -> xr.Dataset | xr.DataArray:
    midx = obj.indexes[midx_name]
    to_drop = [midx.name] + midx.names
    df = midx.to_frame()
    df.insert(0, level_name, [value]*len(midx)) # in place!
    midx = pd.MultiIndex.from_frame(df)
    coords = xr.Coordinates.from_pandas_multiindex(midx, dim=midx_name)
    return obj.drop_vars(to_drop).assign_coords(coords)


def assign_levels(
    obj: xr.Dataset | xr.DataArray,
    levels: dict[str, npt.ArrayLike] | None = None,
    **levels_kwargs: npt.ArrayLike,
) -> xr.Dataset | xr.DataArray:
    """Assign new values to levels of MultiIndexes in ``obj``

    Parameters
    ----------
    obj
        An ``xarray`` object with at least one MultiIndex
    levels, optional
        A mapping whose keys are the names of the levels and whose values are the
        levels to assign. The mapping will be passed to :py:meth:`xarray.DataArray.assign_coords`
        (or the :py:class:`xarray.Dataset` equivalent).

    Returns
    -------
        A new object with the new level values replacing the old level values.
    Raises
    ------
    ValueError
        If levels are provided in both keyword and dictionary form.
    """
    if levels_kwargs != {}:
        if levels is not None:
            raise ValueError(
                "cannot specify both keyword and positional arguments to assign_levels"
            )
        levels = levels_kwargs
    lvl_names = list(levels.keys())
    midxs = set(obj.indexes[lvl].name for lvl in lvl_names)
    # Using sum() to ravel a list of lists
    to_restore = sum([list(obj.indexes[midx].names) for midx in midxs], [])
    obj = obj.reset_index(*midxs)
    obj = obj.assign_coords(levels)
    return obj.set_xindex(to_restore)


def open_frames(path):
    """Opens a NetCDF4 file saved by shnitsel-tools, specially interpreting certain attributes.

    Parameters
    ----------
    path
        The path of the file to open.

    Returns
    -------
        An :py:class:`xarray.Dataset` with any MultiIndex restored.

    Raises
    ------
    FileNotFoundError
        If there is is nothing at ``path``, or ``path`` is not a file.
    ValueError (or other exception)
        Raised by the underlying `h5netcdf <https://h5netcdf.org/>`_ engine if the file is corrupted.
    """
    # The error raised for a missing file can be misleading
    try:
        frames = xr.open_dataset(path)
    except ValueError as err:
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        else:
            raise err

    # Restore MultiIndexes
    indicator = '_MultiIndex_levels_from_attrs'
    if frames.attrs.get(indicator, False):
        # New way: get level names from attrs
        del frames.attrs[indicator]
        for k, v in frames.attrs.items():
            if k.startswith('_MultiIndex_levels_for_'):
                frames = frames.set_xindex(v)
                del frames.attrs[k]
    else:
        # Old way: hardcoded level names
        tcoord = None
        if 'time' in frames.coords:
            tcoord = 'time'
        elif 'ts' in frames.coords:
            tcoord = 'ts'

        if tcoord is not None:
            frames = frames.set_xindex(['trajid', tcoord])
        frames = frames.set_xindex(['from', 'to'])

    return frames


def save_frames(frames, path, complevel=9):
    """Save a ``Dataset``, presumably (but not necessarily) consisting of frames of trajectories, to a file at ``path``.

    Parameters
    ----------
    frames (omit if using accessor)
        The ``Dataset`` to save
    path
        The path at which to save it
    complevel, optional
        The level of ``gzip`` compression which will be applied to all variables in the ``Dataset``, by default 9

    Notes
    -----
    This function/accessor method wraps :py:meth:`xarray.Dataset.to_netcdf` but not :py:func:`numpy.any`.
    """
    frames = frames.copy()  # Shallow copy to avoid adding attrs etc. to original
    encoding = {
        var: {"compression": "gzip", "compression_opts": complevel} for var in frames
    }

    # NetCDF does not support booleans
    for data_var in frames.data_vars:
        if np.issubdtype(frames.data_vars[data_var].dtype, np.bool_):
            frames = frames.assign({data_var: frames.data_vars[data_var].astype('i1')})
    for coord in frames.coords:
        if np.issubdtype(frames.coords[coord].dtype, np.bool_):
            frames = frames.assign_coords({coord: frames.coords[coord].astype('i1')})
    for attr in frames.attrs:
        if np.issubdtype(np.asarray(frames.attrs[attr]).dtype, np.bool_):
            frames.attrs[attr] = int(frames.attrs[attr])

    # NetCDF does not support MultiIndex
    # Keep a record of the level names in the attrs
    midx_names = []
    for name, index in frames.indexes.items():
        if index.name == name and len(index.names) > 1:
            midx_names.append(name)
            midx_levels = list(index.names)
            frames.attrs[f'_MultiIndex_levels_for_{name}'] = midx_levels
    frames.attrs['_MultiIndex_levels_from_attrs'] = 1
    frames.reset_index(midx_names).to_netcdf(path, engine='h5netcdf', encoding=encoding)


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

def mgroupby(
    obj: xr.Dataset | xr.DataArray, levels: Sequence[str]
) -> DataArrayGroupBy | DatasetGroupBy:
    """Group a Dataset or DataArray by several levels of a MultiIndex it contains.

    Parameters
    ----------
    obj
        The :py:mod:`xr` object to group
    levels
        Names of MultiIndex levels all belonging to the *same* MultiIndex

    Returns
    -------
        The grouped object, which behaves as documented at :py:meth:`xr.Dataset.groupby`
        and `xr.DataArray.groupby` with the caveat that the specified levels have been
        "flattened" into a single Multiindex level of tuples.

    Raises
    ------
    ValueError
        If no MultiIndex is found, or if the named levels belong to different MultiIndexes.

    Warnings
    --------
    The function does not currently check whether the levels specified are really levels
    of a MultiIndex, as opposed to names of non-MultiIndex indexes.
    """
    # Ensure all levels belong to the same multiindex
    midxs = set(obj.indexes[lvl].name for lvl in levels)
    if len(midxs) == 0:
        raise ValueError("No index found")
    elif len(midxs) > 1:
        raise ValueError(
            f"The levels provided belong to multiple independent MultiIndexes: {midxs}"
        )
    midx = midxs.pop()
    new_name = ','.join(levels)
    # Flatten the specified levels to tuples and group the resulting object
    return flatten_levels(obj, midx, levels, new_name=new_name).groupby(new_name)


def msel(obj: xr.Dataset | xr.DataArray, **kwargs) -> xr.Dataset | xr.DataArray:
    tuples = list(zip(*kwargs.items()))
    ks, vs = list(tuples[0]), list(tuples[1])
    # Find correct index and levels
    for coord in obj.coords:
        if set(obj.coords[coord].data) <= set(ks):
            levels = obj.indexes[coord].names
            break
    else:
        raise ValueError(f"Couldn't find a coordinate containing all keys {ks}")
    to_reset = list(set(levels) - {coord})
    # Construct selector
    selectee = xr.DataArray(vs, coords=[(coord, ks)])
    # Perform selection
    return (
        selectee.sel({coord: obj.coords[coord]})
        .reset_index(to_reset)
        .set_xindex(levels)
    )

@needs(dims={'frame'}, coords_or_vars={'trajid'})
def sel_trajs(
    frames: xr.Dataset | xr.DataArray,
    trajids_or_mask: Sequence[int] | Sequence[bool],
    invert=False,
) -> xr.Dataset | xr.DataArray:
    """Select trajectories using a list of trajectories IDs or a boolean mask

    Parameters
    ----------
    frames
        The :py:class:`xr.Dataset` from which a selection is to be drawn
    trajids_or_mask
        Either
            - A sequences of integers representing trajectory IDs to be included, in which
              case the trajectories **may not be returned in the order specified**.
            - Or a sequence of booleans, each indicating whether the trajectory with an ID
              in the corresponding entry in the ``Dataset``'s ``trajid_`` coordinate
              should be included
    invert, optional
        Whether to invert the selection, i.e. return those trajectories not specified, by default False

    Returns
    -------
        A new :py:class:`xr.Dataset` containing only the specified trajectories

    Raises
    ------
    NotImplementedError
        when an attempt is made to index an :py:class:`xr.Datset` without a
        ``trajid_`` dimension/coordinate using a boolean mask
    TypeError
        If ``trajids_or_mask`` has a dtype other than integer or boolean
    """
    trajids_or_mask = np.atleast_1d(trajids_or_mask)
    trajids: npt.NDArray | xr.DataArray
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


@needs(dims={'frame'}, coords_or_vars={'trajid'})
def sel_trajids(frames: xr.Dataset, trajids: npt.ArrayLike, invert=False) -> xr.Dataset:
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

def unstack_trajs(frames: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    """Unstack the ``frame`` MultiIndex so that ``trajid`` and ``time`` become
    separate dims. Wraps the :py:meth:`xarray.Dataset.unstack` method.

    Parameters
    ----------
    frames
        An :py:class:`xarray.Dataset` with a ``frame`` dimension associated with
        a MultiIndex coordinate with levels named ``trajid`` and ``time``. The
        Dataset may also have a ``trajid_`` dimension used for variables and coordinates
        that store information pertaining to each trajectory in aggregate; this will be
        aligned along the ``trajid`` dimension of the unstacked Dataset.

    Returns
    -------
        An :py:class:`xarray.Dataset` with independent ``trajid`` and ``time``
        dimensions.
    """
    per_traj_coords = {
        k: v.rename(trajid_='trajid')
        for k, v in dict(frames.coords).items()
        if 'trajid_' in v.dims and 'frame' not in v.dims
    }
    per_traj_vars = {
        k: v.rename(trajid_='trajid')
        for k, v in dict(frames.data_vars).items()
        if 'trajid_' in v.dims and 'frame' not in v.dims
    }
    per_time_coords = {
        k: v.rename(time_='time')
        for k, v in dict(frames.coords).items()
        if 'time_' in v.dims and 'frame' not in v.dims
    }
    per_time_vars = {
        k: v.rename(time_='time')
        for k, v in dict(frames.data_vars).items()
        if 'time_' in v.dims and 'frame' not in v.dims
    }
    to_drop = to_drop = (
        list(per_traj_coords)
        + list(per_traj_vars)
        + list(per_time_coords)
        + list(per_time_vars)
    )

    # Don't re-add to unstacked dataset
    if 'trajid_' in per_traj_coords:
        del per_traj_coords['trajid_']
    if 'time_' in per_time_coords:
        del per_time_coords['time_']

    res = (
        frames.drop_vars(to_drop)
        .assign_coords({'is_frame': ('frame', np.ones(frames.sizes['frame']))})
        .unstack('frame')
        .assign_coords(per_traj_coords)
        .assign_coords(per_time_coords)
        .assign(per_traj_vars)
        .assign(per_time_vars)
    )
    res['is_frame'] = res['is_frame'].fillna(0).astype(bool)
    return res


def stack_trajs(unstacked: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    """Stack the ``trajid`` and ``time`` dims of an unstacked Dataset
    into a MultiIndex along a new dimension called ``frame``.
    Wraps the :py:meth:`xarray.Dataset.stack` method.

    Parameters
    ----------
    frames
        An :py:class:`xarray.Dataset` with independent ``trajid`` and ``time``
        dimensions.

    Returns
    -------
        An :py:class:`xarray.Dataset` with a ``frame`` dimension associated with
        a MultiIndex coordinate with levels named ``trajid`` and ``time``. Those variables
        and coordinates which only depended on one of ``trajid``
        or ``time`` but not the other in the unstacked Dataset, will be aligned along new
        dimensions named ``trajid_`` and ``time_``. The new dimensions ``trajid_`` and
        ``time_`` will be independent of the ``frame`` dimension and its ``trajid`` and
        ``time`` levels.
    """
    per_traj_coords = {
        k: v.rename(trajid='trajid_')
        for k, v in dict(unstacked.coords).items()
        if 'trajid' in v.dims and 'time' not in v.dims and v.name != 'trajid'
    }
    per_traj_vars = {
        k: v.rename(trajid='trajid_')
        for k, v in (dict(unstacked.data_vars)).items()
        if 'trajid' in v.dims and 'time' not in v.dims
    }
    per_time_coords = {
        k: v.rename(time='time_')
        for k, v in dict(unstacked.coords).items()
        if 'time' in v.dims and 'trajid' not in v.dims and v.name != 'time'
    }
    per_time_vars = {
        k: v.rename(time='time_')
        for k, v in (dict(unstacked.data_vars)).items()
        if 'time' in v.dims and 'trajid' not in v.dims
    }
    to_drop = (
        list(per_traj_coords)
        + list(per_traj_vars)
        + list(per_time_coords)
        + list(per_time_vars)
    )
    per_traj_coords['trajid_'] = unstacked.coords['trajid'].rename(trajid='trajid_')
    per_time_coords['time_'] = unstacked.coords['time'].rename(time='time_')

    res = unstacked.drop_vars(to_drop).stack({'frame': ['trajid', 'time']})
    return (
        res.isel(frame=res.is_frame)
        .drop_vars('is_frame')
        .assign_coords(per_traj_coords)
        .assign_coords(per_time_coords)
        .assign(per_traj_vars)
        .assign(per_time_vars)
    )