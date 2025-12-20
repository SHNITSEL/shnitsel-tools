from logging import warning, info
from numbers import Number
from typing import Literal

import numpy as np
import xarray as xr

from shnitsel.core.typedefs import Frames, Stacked, Unstacked
from shnitsel.data.multi_indices import sel_trajs, stack_trajs, unstack_trajs
from shnitsel.data.trajectory_format import Trajectory


def is_stacked(obj):
    if 'frame' in obj.dims and {'trajid', 'time'}.issubset(obj.coords):
        return True
    elif 'trajid' in obj.dims or 'time' in obj.dims:
        return False
    else:
        raise ValueError(
            "The mask argument should be trajectories, either stacked or unstacked"
        )

def ensure_unstacked(obj):
    if is_stacked(obj):
        return unstack_trajs(obj)
    else:
        return obj


def ensure_stacked(obj):
    if not is_stacked(obj):
        return stack_trajs(obj)
    else:
        return obj


def da_or_data_var(ds_or_da, var_name):
    if hasattr(ds_or_da, 'data_vars'):
        return ds_or_da[var_name]
    else:
        return ds_or_da


def get_unstacked_coord(obj, coord_name):
    if coord_name in obj.dims:
        return obj.coords[coord_name]
    else:
        dim_name = obj.indexes[coord_name].name
        return obj.coords[dim_name].unstack(dim_name)[coord_name]


########################
# Formats we will need #
########################


def cum_mask_from_filtranda(filtranda: Stacked | Unstacked) -> Unstacked:
    filtranda = ensure_unstacked(filtranda)
    res = (filtranda < filtranda.coords['thresholds']).cumprod('time').astype(bool)
    return res


def cum_mask_from_dataset(ds: Stacked | Unstacked) -> Unstacked:
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


def assign_mask(ds):
    "Does not change stacked status"
    mask = cum_mask_from_dataset(ds)
    if is_stacked(ds):
        mask = stack_trajs(mask)
    return ds.assign(is_good_frame=mask)


def true_upto(mask, dim):
    if is_stacked(mask):
        mask = unstack_trajs(mask).fillna(False)
    shifted_coord = np.concat([[-1], mask.coords[dim].data])
    indexer = mask.cumprod(dim).sum(dim).astype(int)
    res = np.take(shifted_coord, indexer.data)
    return indexer.copy(data=res)


def cutoffs_from_mask(mask: Stacked | Unstacked) -> Unstacked:
    if is_stacked(mask):
        mask = unstack_trajs(mask)
    good_upto = true_upto(mask, 'time')
    good_throughout = mask.all('time')
    good_upto.name = 'good_upto'
    return good_upto.assign_coords(good_throughout=good_throughout)


def cutoffs_from_filtranda(filtranda: xr.DataArray) -> xr.DataArray:
    """
    Does not change stacked status
    """
    thresholds = filtranda.coords['thresholds']
    is_good_frame = (filtranda < thresholds).astype(bool)
    return cutoffs_from_mask(is_good_frame)


def cutoffs_from_dataset(ds) -> Unstacked:
    """
    Returns a da containing cutoff times (the same as the good_upto data_var)
    and with a coord called good_throughout

    The returned object has dimensions {'trajid', 'criterion'}.
    In that sense it is unstacked.
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

def assign_cutoffs(ds):
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


def omit(ds:Trajectory):
    cutoffs = cutoffs_from_dataset(ds)
    good_throughout = cutoffs['good_throughout']
    selection = good_throughout.all('criterion')
    res = sel_trajs(ds, selection)
    _log_omit(before=ds, after=res)
    return res


def truncate(ds:Trajectory):
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


def transect(ds: Trajectory|Frames, cutoff: float):
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


def cum_max_quantiles(obj, quantiles=None):
    if quantiles is None:
        quantiles = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]

    da = da_or_data_var(obj, 'filtranda')
    da = ensure_unstacked(da)
    da = da.fillna(0)
    time_axis = da.get_axis_num('time')

    cum_max = da.copy(data=np.maximum.accumulate(da.data, axis=time_axis))
    return cum_max.quantile(quantiles, 'trajid')


