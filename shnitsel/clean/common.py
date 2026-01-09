from logging import warning, info
from numbers import Number
from typing import Hashable, Iterable, Literal

import numpy as np
import xarray as xr

from shnitsel.core.typedefs import Frames, Stacked, Unstacked
from shnitsel.data.multi_indices import sel_trajs, stack_trajs, unstack_trajs
from shnitsel.data.trajectory_format import Trajectory
from shnitsel._contracts import needs


def is_stacked(obj: Stacked | Unstacked):
    """Test whether ``obj`` is stacked trajectories (with a 'frames' dimension, and 'trajid' and 'time' coordinates)
    or unstacked trajectories (with 'trajid' and 'time' dimensions).

    Parameters
    ----------
    obj
        A Dataset or DataArray.

    Returns
    -------
        True if ``obj`` is stacked, False if it is unstacked.

    Raises
    ------
    ValueError
        If ``obj`` is neither stacked nor unstacked.
    """
    if 'frame' in obj.dims and {'trajid', 'time'}.issubset(obj.coords):
        return True
    elif 'trajid' in obj.dims or 'time' in obj.dims:
        return False
    else:
        raise ValueError(
            "The mask argument should be trajectories, either stacked or unstacked"
        )

def ensure_unstacked(obj: Stacked | Unstacked) -> Unstacked:
    """Unstack ``obj`` unless already unstacked

    Parameters
    ----------
    obj
        A Dataset or DataArray.

    Returns
    -------
        ``obj`` in unstacked format
    """
    if is_stacked(obj):
        return unstack_trajs(obj)
    else:
        return obj


def ensure_stacked(obj: Stacked | Unstacked) -> Stacked:
    """Stack ``obj`` unless already stacked.

    Parameters
    ----------
    obj
        A Dataset or DataArray.

    Returns
    -------
        ``obj`` in unstacked format
    """
    if not is_stacked(obj):
        return stack_trajs(obj)
    else:
        return obj


def da_or_data_var(
    ds_or_da: xr.Dataset | xr.DataArray, var_name: Hashable
) -> xr.DataArray:
    """Return ``ds_or_da`` if it is a DataArray; otherwise, extract the
    data_var named ``var_name`` from it.

    Parameters
    ----------
    ds_or_da
        A Dataset containing a data_var named ``var_name``; or a DataArray
    var_name
        The name of the data_var to get

    Returns
    -------
        A DataArray
    """
    if hasattr(ds_or_da, 'data_vars'):
        return ds_or_da[var_name]
    else:
        return ds_or_da


def get_unstacked_coord(
    obj: xr.Dataset | xr.DataArray, coord_name: Hashable
) -> xr.DataArray:
    """Get a coordinate from ``obj`` in its unstacked form

    Parameters
    ----------
    obj
        A Dataset or DataArray
    coord_name
        The name of the coord to extract

    Returns
    -------
        An unstacked coordinate
    """
    if coord_name in obj.dims:
        return obj.coords[coord_name]
    else:
        dim_name = obj.indexes[coord_name].name
        return obj.coords[dim_name].unstack(dim_name)[coord_name]


########################
# Formats we will need #
########################


def cum_mask_from_filtranda(filtranda: Stacked | Unstacked) -> Unstacked:
    """Get a mask indicating whether all frames of a trajectory upto and
    including a given frame have abided by filtration criteria.

    Parameters
    ----------
    filtranda
        A DataArray of values by which to filter along with a 'threshold`
        coord containing thresholds to compare them to; separate criteria
        lie along the 'criterion' dimension

    Returns
    -------
        A mask, which for a given trajectory switches from True to False
        exactly once, at the validity cutoff.
    """
    filtranda = ensure_unstacked(filtranda)
    res = (filtranda < filtranda.coords['thresholds']).cumprod('time').astype(bool)
    return res


def cum_mask_from_dataset(ds: Stacked | Unstacked) -> Unstacked:
    """Get a mask indicating whether all frames of a trajectory upto and
    including a given frame have abided by filtration criteria.

    Parameters
    ----------
    ds
        A Dataset containing one of

            - a coordinate called 'is_good_frame'
            - a coordinate called 'good_upto'
            - a data_var called 'filtranda'

    Returns
    -------
        A mask in unstacked form (i.e. 'trajid' and 'time' are separate dimensions)
    """
    # TODO: only consider filtranda?
    if 'is_good_frame' in ds:
        mask = ensure_unstacked(ds['is_good_frame'])
    elif 'good_upto' in ds:
        # time_coord = get_unstacked_coord(ds, 'time')
        ds = ensure_unstacked(ds)
        mask = ds.coords['time'] <= ds['good_upto']
        mask = mask.assign_coords(is_frame=ds.coords['is_frame'])
    elif 'filtranda' in ds:
        return cum_mask_from_filtranda(ds['filtranda'])

    return mask.cumprod('time').astype(bool)


# TODO: unused, consider removing.
def _assign_mask(ds):
    mask = cum_mask_from_dataset(ds)
    if is_stacked(ds):
        mask = stack_trajs(mask)
    return ds.assign(is_good_frame=mask)


def true_upto(mask: xr.DataArray, dim: Hashable) -> xr.DataArray:
    """Find the last index label along dimension 'dim'
    upto which all values of ``mask`` are ``True``;
    independently vectorized over all other dimensions.

    Parameters
    ----------
    mask
        A DataArray of booleans
    dim
        The name of the dimension along which to choose indexes

    Returns
    -------
        A DataArray similar metadata to ``mask``, but with
        dimension ``dim`` removed
    """
    if is_stacked(mask):
        mask = unstack_trajs(mask).fillna(False)
    shifted_coord = np.concat([[-1], mask.coords[dim].data])
    indexer = mask.cumprod(dim).sum(dim).astype(int)
    res = np.take(shifted_coord, indexer.data)
    return indexer.copy(data=res)


def cutoffs_from_mask(mask: Stacked | Unstacked) -> Unstacked:
    """Get cutoffs, i.e. times upto which a trajectory is valid,
    from a boolean mask

    Parameters
    ----------
    mask
        A boolean mask DataArray indicating whether each frame,
        taken independently, is valid

    Returns
    -------
    A DataArray of cutoff times, with similar dimensions to ``filtranda``,
    but missing the 'trajid' dimension
    """
    if is_stacked(mask):
        mask = unstack_trajs(mask)
    good_upto = true_upto(mask, 'time')
    good_throughout = mask.all('time')
    good_upto.name = 'good_upto'
    return good_upto.assign_coords(good_throughout=good_throughout)


def cutoffs_from_filtranda(filtranda: xr.DataArray) -> xr.DataArray:
    """Get cutoffs, i.e. times upto which a trajectory is valid,
    from thresholds and filtranda

    Parameters
    ----------
    filtranda
        A DataArray with stacked or unstacked ensemble dimensions
        as well as a 'criterion' dimension;
        with a 'thresholds' coord along the 'criterion'
        dimension

    Returns
    -------
    A DataArray of cutoff times, with similar dimensions to ``filtranda``,
    but missing the 'trajid' dimension

    Notes
    -----
    Does not change stacked status
    """
    thresholds = filtranda.coords['thresholds']
    is_good_frame = (filtranda < thresholds).astype(bool)
    return cutoffs_from_mask(is_good_frame)


def cutoffs_from_dataset(ds: xr.Dataset) -> Unstacked:
    """Get cutoffs, i.e. times upto which a trajectory is valid,
    from an ensemble Dataset containing filtration information

    Parameters
    ----------
    ds
        A Dataset containing either:

            - a 'good_upto' data_var and a 'good_throughout' coordinate
            - a 'filtranda' data_var with a 'threshold' coordinate

    Returns
    -------
        A DataArray containing cutoff times (the same as the good_upto data_var)
        and with a coord called good_throughout
        The returned object has dimensions {'trajid', 'criterion'}.

    Raises
    ------
    ValueError
        If there is no filtration information in the Dataset
    """
    if 'good_upto' in ds.data_vars:
        res = ds.data_vars['good_upto']
        if 'good_throughout' in res.coords:
            if 'filtranda' in ds.data_vars:
                warning(
                    "data_vars 'filtranda' and 'good_upto' present in "
                    "the same dataset; ignoring 'filtranda'"
                )
            if 'trajid_' in res.dims and 'trajid' not in res.dims:
                res = res.rename(trajid_='trajid')
            return res
        else:
            warning(
                "data_var 'good_upto' is missing expected coord "
                "'good_throughout'; will recalculate."
            )
    elif 'filtranda' in ds.data_vars:
        return cutoffs_from_filtranda(ds.data_vars['filtranda'])
    else:
        raise ValueError(
            "Please set data_vars 'filtranda' and 'thresholds', "
            "or alternatively supply cutoffs directly using data_var 'good_upto'"
        )

# TODO: unused, remove?
def assign_cutoffs(ds):
    """Assign cutoffs to a dataset

    Parameters
    ----------
    ds
        A Dataset acceptable to
        :py:func:`shnitsel.clean.common.cutoffs_from_dataset`

    Returns
    -------
        ``ds`` with an additional Variable 'good_upto', which
        also adds a coord 'good_throughout'
    """
    cutoffs = cutoffs_from_dataset(ds)
    if is_stacked(ds):
        cutoffs = cutoffs.rename(trajid='trajid_')
    return ds.assign(good_upto=cutoffs)


####################
# Action functions #
####################

# All the action functions take a dataset
# They can use the functions above to get the info they need


def _log_omit(before, after):
    kept = set(after.trajid.values.tolist())
    omitted = set(before.trajid.values.tolist()).difference(kept)
    info(
        f"Kept {len(kept)} trajectories, IDs: {kept}; \n"
        f"Dropped {len(omitted)} trajectories, IDs: {omitted}"
    )

@needs(data_vars={'filtranda'})
def omit(ds:Trajectory):
    """Remove all trajectories containing even one bad frames from an ensemble.
    Badness is determined based on the information considered by
    :py:func:`shnitsel.clean.cutoffs_from_dataset`.

    Parameters
    ----------
    ds
        An ensemble of trajectories

    Returns
    -------
        The filtered ensemble

    See also
    --------
    :py:func:`shnitsel.clean.truncate`
    :py:func:`shnitsel.clean.transect`
    """
    cutoffs = cutoffs_from_dataset(ds)
    good_throughout = cutoffs['good_throughout']
    selection = good_throughout.all('criterion')
    res = sel_trajs(ds, selection)
    _log_omit(before=ds, after=res)
    return res


@needs(data_vars={'filtranda'})
def truncate(ds:Trajectory):
    """Cut off each trajectory in an ensemble just before the first bad frame.
    Badness is determined based on the information considered by
    :py:func:`shnitsel.clean.cutoffs_from_dataset`.

    Parameters
    ----------
    ds
        An ensemble of trajectories

    Returns
    -------
        The filtered ensemble

    See also
    --------

        - :py:func:`shnitsel.clean.omit`
        - :py:func:`shnitsel.clean.transect`
    """
    if is_stacked(ds):
        cutoffs = cutoffs_from_dataset(ds).min('criterion')
        sel = cutoffs.sel(trajid=ds.coords['trajid'])
        stacked_mask = ds.coords['time'].data < sel.data
        res = ds.isel(frame=stacked_mask)
        nbefore = ds.sizes['frame']
        nafter = res.sizes['frame']
    else:
        unstacked_mask = cum_mask_from_dataset(ds).all('criterion')
        res = ds.assign_coords(is_frame=ds.coords['is_frame'] & unstacked_mask)
        nbefore = ds.coords['is_frame'].sum().item()
        nafter = res.coords['is_frame'].sum().item()
    info(f"Kept {nafter} frames; dropped {nbefore - nafter}")
    return res


@needs(data_vars={'filtranda'})
def transect(ds: Trajectory|Frames, cutoff: float):
    """Cut off all trajectories at a specified time, discarding trajectories
    which have bad frames earlier.
    Badness is determined based on the information considered by
    :py:func:`shnitsel.clean.cutoffs_from_dataset`.

    Parameters
    ----------
    ds
        An ensemble of trajectories
    cutoff
        The time unitl which all trajectories should last

    Returns
    -------
        The filtered ensemble

    See also
    --------

        - :py:func:`shnitsel.clean.truncate`
        - :py:func:`shnitsel.clean.omit`
    """
    if is_stacked(ds):
        ds = unstack_trajs(ds)
    nbefore = ds.coords['is_frame'].sum().item()
    ds = ds.loc[{'time': slice(float(cutoff))}]
    nafter1 = ds.coords['is_frame'].sum().item()
    good_upto = cutoffs_from_dataset(ds)
    traj_selection = (
        (good_upto >= cutoff).all('criterion')
        &
        # the following must be calculated after time-slicing.
        ds.coords['is_frame'].all('time')
    )
    res = ds.isel({'trajid': traj_selection})
    trajs_dropped = set(ds['trajid'].values.tolist()).difference(
        res['trajid'].values.tolist()
    )
    nafter2 = res.coords['is_frame'].sum().item()
    info(
        f"Kept {nbefore} frames; dropped {nafter1} frames past cutoff;\n"
        f"dropped {len(trajs_dropped)} trajectories "
        f"({nafter1 - nafter2} additional frames) "
        f"that did not reach the cutoff, IDs: {trajs_dropped}"
    )
    return res


def dispatch_cut(
    frames,
    cut: Literal['truncate', 'omit', False] | Number = 'truncate',
):
    """Filter ``frames`` in a manner depending on ``cut``

    Parameters
    ----------
    frames
        An ensemble of trajectories
    cut
        One of

            - 'truncate', in which case trajectories will be cut off
                just before the first bad frame using
                :py:func:`shnitsel.clean.truncate`
            - 'omit', in which case trajectories containing a bad frame
                will be removed entirely using
                :py:func:`shnitsel.clean.omit`
            - a number, in which case this number is interpreted as the
                desired duration of all trajectories; trajectories in
                which a bad frame occurs earlier are discarded using
                :py:func:`shnitsel.clean.transect`
    Returns
    -------
        The filtered ensemble

    Raises
    ------
    ValueError
        If ``cut`` is not one of the accepted values
    """
    if not cut:
        return frames.assign(good_upto=cutoffs_from_dataset(frames))
    elif cut == 'truncate':
        return truncate(frames)
    elif cut == 'omit':
        return omit(frames)
    elif isinstance(cut, Number):
        return transect(frames, cut)
    else:
        raise ValueError(
            "`cut` should be one of {'truncate', 'omit'}, or a number, " f"not {cut}"
        )


###########################################
# Formats directly prerequisite for plots #
###########################################


def cum_max_quantiles(
    obj: xr.Dataset | xr.DataArray, quantiles: Iterable[float] | None = None
) -> xr.DataArray:
    """Quantiles of cumulative maxima

    Parameters
    ----------
    obj
        A DataArray, or a Dataset with a data_var 'filtranda';
        either way, the Variable should have dimensions and
        coordinates corresponding to a
        (stacked or unstacked) ensemble of trajectories.
    quantiles, optional
        Which quantiles to calculate,
        by default ``[0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]``

    Returns
    -------
        A DataArray with 'quantile' and 'time' dimensions;
        'frame' or 'trajid' dimensions will have been removed;
        other dimensions remain unaffected.

    See also
    --------
    The data returned by this function is intended for consumption by
    :py:func:`shnitsel.vis.plot.filtration.check_thresholds`
    """
    if quantiles is None:
        quantiles = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]

    da = da_or_data_var(obj, 'filtranda')
    da = ensure_unstacked(da)
    da = da.fillna(0)
    time_axis = da.get_axis_num('time')

    cum_max = da.copy(data=np.maximum.accumulate(da.data, axis=time_axis))
    return cum_max.quantile(quantiles, 'trajid')


