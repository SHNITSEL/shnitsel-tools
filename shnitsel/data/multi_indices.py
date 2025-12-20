from __future__ import annotations

import itertools
from typing import Callable, Sequence, TypeVar

import numpy.typing as npt

from xarray.core.groupby import DatasetGroupBy, DataArrayGroupBy

import xarray as xr
import numpy as np
import pandas as pd

from shnitsel._api_info import internal

from .._contracts import needs

DatasetOrArray = TypeVar("DatasetOrArray", bound=xr.Dataset | xr.DataArray)


@internal()
def midx_combs(values: pd.core.indexes.base.Index | list, name: str | None = None):
    """Helper function to create a Multi-index based dimension coordinate for an xarray
    from all (unordered) pairwise combinations of entries in `values`

    Args:
        values (pd.core.indexes.base.Index | list): The source values to generate pairwise combinations for
        name (str | None, optional): Optionally a name for the resulting combination dimension. Defaults to None.

    Raises:
        ValueError: If no name was provided and the name could not be extracted from the `values` parameter

    Returns:
        xr.Coordinates: The resulting coordinates object.
    """
    if name is None:
        if hasattr(values, 'name'):
            # if `values` is a `pandas.core.indexes.base.Index`
            # extract its name
            name = values.name
        else:
            raise ValueError("need to specify name if values lack name attribute")

    comb_name = f'{name}comb'

    return xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(
            itertools.combinations(values, 2),
            names=[f'{comb_name}_from', f'{comb_name}_to'],
        ),
        dim=comb_name,
    )


def flatten_midx(
    obj: DatasetOrArray, idx_name: str, renamer: Callable | None = None
) -> DatasetOrArray:
    """Function to flatten a multi-index into a flat index.

    Has the option to provide a custom renaming function

    Args:
        obj (xr.Dataset | xr.DataArray): The object with the index intended to be flattened
        idx_name (str): The name of the index to flatten.
        renamer (callable | None, optional): An optional function to carry out the renaming of the combined entry from individual entries. Defaults to None.

    Returns:
        xr.Dataset | xr.DataArray: The refactored object without the original index coordinates but with a combined index instead
    """
    midx = obj.indexes[idx_name]
    to_drop = midx.names + [midx.name]
    fidx = midx.to_flat_index()
    if renamer is not None:
        fidx = [renamer(x, y) for x, y in fidx]
    return obj.drop_vars(to_drop).assign_coords({idx_name: fidx})


def flatten_levels(
    obj: DatasetOrArray,
    idx_name: str,
    levels: Sequence[str],
    new_name: str | None = None,
    position: int = 0,
    renamer: Callable | None = None,
) -> DatasetOrArray:
    """Flatten specified levels of a MultiIndex into tuples occupying
    a single MultiIndex level

    Parameters
    ----------
    obj
        A Dataset or DataArray with at least one MultiIndex
    idx_name
        The name of the MultiIndex
    levels
        Which levels to flatten
    new_name, optional
        The name of the single resulting index, by default None
    position, optional
        The position of the resulting level in the MultiIndex, by default 0
    renamer, optional
        A Callable to compute the values in the new level as a
        function of the values in the original separate levels, by default None

    Returns
    -------
        An object differing from ``obj`` only in the flattening of specified levels

    Raises
    ------
    ValueError
        If the specified index is associated with more than one dimension
        (this should not be possible for a MultiIndex anyway)
    """
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
    obj: DatasetOrArray, midx_name: str, level_name: str, value
) -> DatasetOrArray:
    """Add an outer level to an existing MultiIndex in ``obj``

    Parameters
    ----------
    obj
        A Dataset or DataArray with at least one MultiIndex
    midx_name
        The name of the MultiIndex
    level_name
        The name of the new level
    value
        Values with to populate the new level

    Returns
    -------
        An object differing from ``obj`` only in the addition of the MultiIndex level
    """
    midx = obj.indexes[midx_name]
    to_drop = [midx.name] + midx.names
    df = midx.to_frame()
    df.insert(0, level_name, [value] * len(midx))  # in place!
    midx = pd.MultiIndex.from_frame(df)
    coords = xr.Coordinates.from_pandas_multiindex(midx, dim=midx_name)
    return obj.drop_vars(to_drop).assign_coords(coords)


def assign_levels(
    obj: DatasetOrArray,
    levels: dict[str, npt.ArrayLike] | None = None,
    **levels_kwargs: npt.ArrayLike,
) -> DatasetOrArray:
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
        A new object (of the same type as `obj`) with the new level values replacing the old level values.
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
    # Assignment of DataArrays fails. Workaround:
    for lvl in levels:
        if isinstance(levels[lvl], xr.DataArray):
            lvl_dims = levels[lvl].dims
            assert len(lvl_dims) == 1
            levels[lvl] = (lvl_dims[0], levels[lvl].data)
    lvl_names = list(levels.keys())
    midxs = set(
        obj.indexes[lvl].name
        for lvl in lvl_names
        # The following filter lets this function also assign normal coords:
        if obj.indexes[lvl].name != lvl
    )
    # Using sum() to ravel a list of lists
    to_restore = sum([list(obj.indexes[midx].names) for midx in midxs], [])
    if midxs:
        obj = obj.reset_index(*midxs)
    obj = obj.assign_coords(levels)
    if to_restore:
        obj = obj.set_xindex(to_restore)
    return obj


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


def msel(obj: DatasetOrArray, **kwargs) -> DatasetOrArray:
    """Add data values along a coordinate, chosen based on coordinate values

    Parameters
    ----------
    obj
        A Dataset or DataArray with at least one coordinate containing all
        the values given by the ``kwargs`` parameter name

    Returns
    -------
        The coordinate (presumably unique) from ``obj`` that contains all the parameter
        names in ``kwargs``

    Raises
    ------
    ValueError
        If no coordinate in ``obj`` contains all the parameter names in ``kwargs``
    """
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
    frames: DatasetOrArray,
    trajids_or_mask: Sequence[int] | Sequence[bool],
    invert=False,
) -> DatasetOrArray:
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
def sel_trajids(
    frames: DatasetOrArray, trajids: npt.ArrayLike, invert=False
) -> DatasetOrArray:
    """Select trajectories using a list of trajectories IDs or a boolean mask;
    note that the trajectories may not be returned in the order specified.

    Parameters
    ----------
    frames
        The :py:class:`xr.Dataset` from which a selection is to be drawn
    trajids
        A sequences of integers representing trajectory IDs to be included,
    invert
        Whether to invert the selection, i.e. return those trajectories not specified, by default False

    Returns
    -------
        A new :py:class:`xr.Dataset` containing only the specified trajectories

    Raises
    ------
    KeyError
        If some of the supplied trajectory IDs are not present in the ``trajid`` coordinate
    """
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


@internal()
def unstack_trajs(frames: DatasetOrArray) -> DatasetOrArray:
    """Unstack the ``frame`` MultiIndex so that ``trajid`` and ``time`` become
    separate dims. Wraps the :py:meth:`xarray.Dataset.unstack` method.

    Parameters
    ----------
    frames, DatasetOrArray
        An :py:class:`xarray.Dataset` with a ``frame`` dimension associated with
        a MultiIndex coordinate with levels named ``trajid`` and ``time``. The
        Dataset may also have a ``trajid_`` dimension used for variables and coordinates
        that store information pertaining to each trajectory in aggregate; this will be
        aligned along the ``trajid`` dimension of the unstacked Dataset.

    Returns
    -------
        An :py:class:`xarray.Dataset` with independent ``trajid`` and ``time``
        dimensions. Same type as `frames`
    """
    per_traj_coords = {
        k: v.rename(trajid_='trajid')
        for k, v in dict(frames.coords).items()
        if 'trajid_' in v.dims and 'frame' not in v.dims
    }
    per_time_coords = {
        k: v.rename(time_='time')
        for k, v in dict(frames.coords).items()
        if 'time_' in v.dims and 'frame' not in v.dims
    }
    if hasattr(frames, 'data_vars'):
        has_data_vars = True
        per_traj_vars = {
            k: v.rename(trajid_='trajid')
            for k, v in dict(frames.data_vars).items()
            if 'trajid_' in v.dims and 'frame' not in v.dims
        }
        per_time_vars = {
            k: v.rename(time_='time')
            for k, v in dict(frames.data_vars).items()
            if 'time_' in v.dims and 'frame' not in v.dims
        }
    else:
        has_data_vars = False
        per_traj_vars = []
        per_time_vars = []

    to_drop = to_drop = (
        list(per_traj_coords)
        + list(per_time_coords)
        + list(per_traj_vars)
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
    )
    if has_data_vars:
        res = res.assign(per_traj_vars).assign(per_time_vars)
    res['is_frame'] = res['is_frame'].fillna(0).astype(bool)
    return res


@internal()
def stack_trajs(unstacked: DatasetOrArray) -> DatasetOrArray:
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
    per_time_coords = {
        k: v.rename(time='time_')
        for k, v in dict(unstacked.coords).items()
        if 'time' in v.dims and 'trajid' not in v.dims and v.name != 'time'
    }
    if hasattr(unstacked, 'data_vars'):
        has_data_vars = True
        per_traj_vars = {
            k: v.rename(trajid='trajid_')
            for k, v in (dict(unstacked.data_vars)).items()
            if 'trajid' in v.dims and 'time' not in v.dims
        }
        per_time_vars = {
            k: v.rename(time='time_')
            for k, v in (dict(unstacked.data_vars)).items()
            if 'time' in v.dims and 'trajid' not in v.dims
        }
    else:
        has_data_vars = False
        per_traj_vars = []
        per_time_vars = []
    to_drop = (
        list(per_traj_coords)
        + list(per_traj_vars)
        + list(per_time_coords)
        + list(per_time_vars)
    )
    per_traj_coords['trajid_'] = unstacked.coords['trajid'].rename(trajid='trajid_')
    per_time_coords['time_'] = unstacked.coords['time'].rename(time='time_')

    res = unstacked.drop_vars(to_drop).stack({'frame': ['trajid', 'time']})
    res = (
        res.isel(frame=res.is_frame)
        .drop_vars('is_frame')
        .assign_coords(per_traj_coords)
        .assign_coords(per_time_coords)
    )
    if has_data_vars:
        res = res.assign(per_traj_vars).assign(per_time_vars)
    return res


@needs(dims={'frame'})
def mdiff(da: xr.DataArray) -> xr.DataArray:
    """Take successive differences along the 'frame' dimension

    Parameters
    ----------
    da
        An ``xarray.DataArray`` with a 'frame' dimension corresponding
        to a ``pandas.MultiIndex`` of which the innermost level is 'time'.

    Returns
    -------
        An ``xarray.DataArray`` with the same shape, dimension names etc.,
        but with the data of the (i)th frame replaced by the difference between
        the original (i+1)th and (i)th frames, with zeros filling in for both the
        initial frame and any frame for which time = 0, to avoid taking differences
        between the last and first frames of successive trajectories.
    """
    res = xr.apply_ufunc(
        lambda arr: np.diff(arr, prepend=np.array(arr[..., [0]], ndmin=arr.ndim)),
        da,
        input_core_dims=[['frame']],
        output_core_dims=[['frame']],
    )
    res[{'frame': res['time'] == 0}] = 0
    return res
