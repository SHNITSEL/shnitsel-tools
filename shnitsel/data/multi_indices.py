from __future__ import annotations

import itertools
from typing import Callable, Sequence, TypeVar

import numpy.typing as npt

from xarray.core.groupby import DatasetGroupBy, DataArrayGroupBy

import xarray as xr
import numpy as np
import pandas as pd

from shnitsel.core._api_info import internal

from .._contracts import needs

DatasetOrArray = TypeVar("DatasetOrArray", bound=xr.Dataset | xr.DataArray)

# TODO: FIXME: These functions occasionally return weird errors if provided with datasets of the wrong format


@internal()
def midx_combs(
    values: pd.core.indexes.base.Index | list, name: str | None = None
) -> xr.Coordinates:
    """Helper function to create a Multi-index based dimension coordinate for an xarray
    from all (unordered) pairwise combinations of entries in `values`

    Parameters
    ----------
    values : pd.core.indexes.base.Index | list
        The source values to generate pairwise combinations for
    name : str | None, optional
        Optionally a name for the resulting combination dimension. Defaults to None.

    Raises
    ------
    ValueError
        If no name was provided and the name could not be extracted from the `values` parameter

    Returns
    -------
    xr.Coordinates
        The resulting coordinates object.
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

    Parameters
    ----------
    obj : xr.Dataset | xr.DataArray
        The object with the index intended to be flattened
    idx_name : str
        The name of the index to flatten.
    renamer : callable | None, optional
        An optional function to carry out the renaming of the combined entry from individual entries. Defaults to None.

    Returns
    -------
    xr.Dataset | xr.DataArray
        The refactored object without the original index coordinates but with a combined index instead
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
    obj : DatasetOrArray
        A Dataset or DataArray with at least one MultiIndex
    idx_name : str
        The name of the MultiIndex
    levels : Sequence[str]
        Which levels to flatten
    new_name : str, optional
        The name of the single resulting index, by default None
    position : int, optional
        The position of the resulting level in the MultiIndex, by default 0
    renamer : Callable, optional
        A Callable to compute the values in the new level as a
        function of the values in the original separate levels, by default None

    Returns
    -------
    DatasetOrArray
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
    obj : DatasetOrArray
        A Dataset or DataArray with at least one MultiIndex
    midx_name : str
        The name of the MultiIndex
    level_name : str
        The name of the new level
    value
        Values with to populate the new level

    Returns
    -------
    DatasetOrArray
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
    obj : DatasetOrArray
        An ``xarray`` object with at least one MultiIndex
    levels : dict[str, npt.ArrayLike], optional
        A mapping whose keys are the names of the levels and whose values are the
        levels to assign. The mapping will be passed to :py:meth:`xarray.DataArray.assign_coords`
        (or the :py:class:`xarray.Dataset` equivalent).
    **levels_kwargs
        Keyword arguments to define the levels by instead of providing them as a dict

    Returns
    -------
    DatasetOrArray
        A new object (of the same type as `obj`) with the new level values replacing the old level values.

    Raises
    ------
    ValueError
        If levels are provided in both keyword and dictionary form.

    Notes
    -----
    Propagates attrs irrespective of ``xarray.get_options()['keep_attrs']``
    """
    if levels_kwargs != {}:
        if levels is not None:
            raise ValueError(
                "cannot specify both keyword and positional arguments to assign_levels"
            )
        levels = levels_kwargs
    # Assignment of DataArrays fails. Workaround:
    attrs_by_level = {}
    for lvl in levels:
        if isinstance(levels[lvl], xr.DataArray):
            lvl_dims = levels[lvl].dims
            assert len(lvl_dims) == 1  # Can't have multi-dimensional MultiIndex levels
            attrs_by_level[lvl] = levels[lvl].attrs
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
    for lvl, attrs in attrs_by_level.items():
        obj.coords[lvl].attrs = attrs
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
    obj : xr.Dataset | xr.DataArray
        The :py:mod:`xr` object to group
    levels : Sequence[str]
        Names of MultiIndex levels all belonging to the *same* MultiIndex

    Returns
    -------
    DataArrayGroupBy | DatasetGroupBy
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
    obj : DatasetOrArray
        A Dataset or DataArray with at least one coordinate containing all
        the values given by the ``kwargs`` parameter name
    **kwargs
        Tuples of `key:value` pairs as keyword arguments to select from entries in a multi-index.

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
    # TODO: FIXME: This does not work with a dataset input?
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
    obj: DatasetOrArray,
    trajids_or_mask: Sequence[int] | Sequence[bool],
    invert: bool = False,
) -> DatasetOrArray:
    """Select trajectories using a list of trajectories IDs or a boolean mask

    Parameters
    ----------
    obj : DatasetOrArray
        The :py:class:`xr.Dataset` from which a selection is to be drawn
    trajids_or_mask : Sequence[int] | Sequence[bool]
        Either
            - A sequences of integers representing trajectory IDs to be included, in which
              case the trajectories **may not be returned in the order specified**.
            - Or a sequence of booleans, each indicating whether the trajectory with an ID
              in the corresponding entry in the ``Dataset``'s ``trajid_`` coordinate
              should be included
    invert:bool, optional
        Whether to invert the selection, i.e. return those trajectories not specified, by default False

    Returns
    -------
    DatasetOrArray
        A new :py:class:`xr.Dataset` containing only the specified trajectories

    Raises
    ------
    NotImplementedError
        when an attempt is made to index an :py:class:`xr.Datset` without a
        ``trajid_`` dimension/coordinate using a boolean mask
    TypeError
        If ``trajids_or_mask`` has a dtype other than integer or boolean
    """
    # TODO: FIXME: This Function is currently broken in tests
    trajids_or_mask = np.atleast_1d(trajids_or_mask)
    if not is_stacked(obj):
        return _sel_trajs_unstacked(obj, trajids_or_mask, invert)
    trajids: npt.NDArray | xr.DataArray
    if np.issubdtype(trajids_or_mask.dtype, np.integer):
        trajids = trajids_or_mask
    elif np.issubdtype(trajids_or_mask.dtype, bool):
        mask = trajids_or_mask
        if 'trajid_' in obj.dims:
            trajids = obj['trajid_'][mask]
        else:
            raise NotImplementedError(
                "Indexing trajids with a boolean mask is only supported when the "
                "coordinate 'trajid_' is present, or if `frames` has unstacked trajectories "
                "(i.e. separate dimesions for trajectory and time)"
            )
    else:
        raise TypeError(
            "Only indexing using a boolean mask or integer trajectory IDs is supported; "
            f"the detected dtype was {trajids_or_mask.dtype}"
        )
    return _sel_trajids(frames=obj, trajids=trajids, invert=invert)


@needs(dims={'frame'}, coords_or_vars={'trajid'})
def _sel_trajids(
    frames: DatasetOrArray, trajids: npt.ArrayLike, invert: bool = False
) -> DatasetOrArray:
    """Select trajectories using a list of trajectories IDs;
    note that the trajectories may not be returned in the order specified.

    Parameters
    ----------
    frames : DatasetOrArray
        The :py:class:`xr.Dataset` from which a selection is to be drawn
    trajids : npt.ArrayLike
        A sequences of integers representing trajectory IDs to be included,
    invert :bool, optional
        Whether to invert the selection, i.e. return those trajectories not specified, by default False

    Returns
    -------
    DatasetOrArray
        A new :py:class:`xr.Dataset` containing only the specified trajectories

    Raises
    ------
    KeyError
        If some of the supplied trajectory IDs are not present in the ``trajectory`` coordinate
    """
    trajids = np.atleast_1d(trajids)
    # check that all trajids are valid, as Dataset.sel() would
    if not invert and not (np.isin(trajids, frames.coords['trajectory'])).all():
        missing = trajids[~np.isin(trajids, frames.coords['trajectory'])]
        raise KeyError(
            f"Of the supplied trajectory IDs, {len(missing)} were "
            f"not found in index 'trajid': {missing}"
        )
    mask = frames['trajectory'].isin(trajids)
    if invert:
        mask = ~mask
    res = frames.sel(frame=mask)

    # TODO: FIXME: This needs to be made resilient to stacked and layered sets. Selecting from within the `frames` stack fails the test

    if 'trajectory' in frames.dims:
        actually_selected = np.unique(res['trajectory'])
        res = res.sel(trajectory=actually_selected)
    return res


def _sel_trajs_unstacked(obj, indexer, invert):
    traj_dim_name = (
        'trajectory'
        if 'trajectory' in obj.dims
        else 'trajid'
        if 'trajid' in obj.dims
        else 'atrajectory'
        if 'atrajectory' in obj.dims
        else None
    )
    assert traj_dim_name is not None
    if not invert:
        return obj.loc[{traj_dim_name: indexer}]

    if np.issubdtype(indexer.dtype, np.integer):
        full_coord = obj.coords[traj_dim_name]
        indexer = full_coord[~full_coord.isin(indexer)]
    elif np.issubdtype(indexer.dtype, bool):
        indexer = ~indexer
    else:
        raise ValueError(
            "Could not invert selection, please provide integer labels or a boolean mask"
        )

class dtype_NA:
    """A sentinel value for the ``fill_value`` param in
    :py:func:`shnitsel.data.multi_indices.unstack_trajs`"""

@internal()
def unstack_trajs(frames: DatasetOrArray, fill_value=dtype_NA) -> DatasetOrArray:
    """Unstack the ``frame`` MultiIndex so that ``trajid`` and ``time`` become
    separate dims. Wraps the :py:meth:`xarray.Dataset.unstack` method.

    Parameters
    ----------
    frames : DatasetOrArray
        An :py:class:`xarray.Dataset` with a ``frame`` dimension associated with
        a MultiIndex coordinate with levels named ``trajid`` and ``time``. The
        Dataset may also have a ``trajid_`` dimension used for variables and coordinates
        that store information pertaining to each trajectory in aggregate; this will be
        aligned along the ``trajid`` dimension of the unstacked Dataset.
    fill_value
        The value used to fill in entries that were unspecified in
        stacked format; by default, the dtype's NA value will be used.

    Returns
    -------
    DatasetOrArray
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

    # NOTE: We use this kws approach to avoid importing the default value for fill_value
    # in xr's unstack, which is their internal `xarray.core.dtypes.NA`.
    kws = {'fill_value': fill_value} if fill_value is not dtype_NA else {}
    res = (
        frames.drop_vars(to_drop)
        .assign_coords({'is_frame': ('frame', np.ones(frames.sizes['frame']))})
        .unstack('frame', **kws)
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
    frames : DatasetOrArray
        An :py:class:`xarray.Dataset` with independent ``trajid`` and ``time``
        dimensions.

    Returns
    -------
    DatasetOrArray
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


def is_stacked(obj):
    """Test whether an object has stacked trajectories

    Parameters
    ----------
    obj
        An xarray Dataset/DataArray, or a wrapper around one

    Returns
    -------
        True if ``obj`` shows signs of containing multiple
        trajectories along the same dimension as used for the
        time coordinate.
    """
    from shnitsel.data.dataset_containers.frames import Frames

    TRAJECTORY_DIM_NAMES = {'trajid', 'atrajectory'}

    is_wrapped_stacked = isinstance(obj, Frames)

    c = obj.coords
    traj_coord = c.get(
        'trajid', c.get('atrajectory', c.get('trajectory', xr.DataArray()))
    )
    time_coord = c.get('time', xr.DataArray())
    coords_share_dim = not set(traj_coord.dims).isdisjoint(time_coord.dims)
    return is_wrapped_stacked or coords_share_dim


def ensure_unstacked(obj, fill_value=dtype_NA):
    """Unstack ``obj`` if it contains stacked trajectories

    Parameters
    ----------
    obj
        An xarray Dataset/DataArray, or a wrapper around one
    fill_value
        The value used to fill in entries that were unspecified in
        stacked format; by default, the dtype's NA value will be used.

    Returns
    -------
    unstacked
        The unstacked Dataset/DataArray
    was_stacked
        Whether ``obj`` had stacked trajectories
    """
    was_stacked = is_stacked(obj)
    unstacked = unstack_trajs(obj, fill_value=fill_value) if was_stacked else obj
    return unstacked, was_stacked


@needs(dims={'frame'})
def mdiff(da: xr.DataArray, dim: str | None = None) -> xr.DataArray:
    """Take successive differences along the `dim` dimension

    Parameters
    ----------
    da : xr.DataArray
        An ``xarray.DataArray`` with a dimension `dim` corresponding
        to a ``pandas.MultiIndex`` of which the innermost level is 'time'.
    dim : str, optional
        The dimension along which the successive differences should be calculated.

    Returns
    -------
        An ``xarray.DataArray`` with the same shape, dimension names etc.,
        but with the data of the (i)th frame replaced by the difference between
        the original (i+1)th and (i)th frames, with zeros filling in for both the
        initial frame and any frame for which time = 0, to avoid taking differences
        between the last and first frames of successive trajectories.
    """
    # TODO: FIXME: Tweak documentation to actually reflect what is happening here.
    if dim is None:
        leading_dim = (
            'frame'
            if 'frame' in da.dims
            else ('time' if 'time' in da.dims else da.dims[0])
        )
    else:
        leading_dim = dim

    res = xr.apply_ufunc(
        lambda arr: np.diff(arr, prepend=np.array(arr[..., [0]], ndmin=arr.ndim)),
        da,
        input_core_dims=[[leading_dim]],
        output_core_dims=[[leading_dim]],
    )

    # If the index of `dim` is a MultiIndex, we set the output to zero at boundaries
    if hasattr(da.indexes[leading_dim], 'names'):
        mask = np.zeros_like(da.coords[leading_dim], dtype=bool)
        for level_name in set(da.indexes[leading_dim].names) - {'time'}:
            level = da.coords[level_name]
            # Using manual comparison of shifted values (rather than np.diff)
            # to allow for non-subtractable level entries:
            right_shift = level.shift({leading_dim: 1}, level[0])
            mask |= level != right_shift

        res[{leading_dim: mask}] = 0

    return res
