from logging import warning
from numbers import Number

import numpy as np
import xarray as xr

from shnitsel.core.midx import sel_trajs, unstack_trajs


def is_stacked(obj):
    if 'frame' in obj.dims and {'trajid', 'time'}.issubset(obj.coords):
        return True
    elif {'trajid', 'time'}.issubset(obj.dims):
        return False
    else:
        raise ValueError(
            "The mask argument should be trajectories, either stacked or unstacked"
        )

########################
# Formats we will need #
########################

def cum_mask_from_filtranda(filtranda):
    # isn't this literally just
    return (filtranda < filtranda.coords['thresholds']).cumprod('time').astype(bool)
    # why is this a function
    # also is "cumulative mask" really your best attempt to explain this concept?


def cum_mask_from_cutoffs(good_upto, thresholds):
    return good_upto < thresholds


def cum_mask_from_dataset():
    # dataset might contain filtranda, the mask we're after, or the cutoffs
    # it doesn't matter whether the mask is cumulative, because our
    # accumulation is idempotent
    
    # two versions currently, for stacked and for unstacked
    # use `filtranda` and `thresholds` if available
    # otherwise use `good_upto`
    ... # TODO

def true_upto(mask, dim):
    if is_stacked(mask):
        mask = unstack_trajs(mask).fillna(False)
    shifted_coord = xr.DataArray(np.concat([[-1], mask.coords[dim].data]), dims=dim)
    res = shifted_coord[mask.cumprod(dim).sum(dim)]
    return res


def cutoffs_from_mask(is_good_frame):
    good_upto = true_upto(is_good_frame, 'time')
    good_throughout = is_good_frame.groupby('time').all()
    good_upto.name = 'good_upto'
    return good_upto.assign_coords(good_throughout=good_throughout)
    

def cutoffs_from_filtranda(filtranda):
    thresholds = filtranda.coords['thresholds']
    is_good_frame = (filtranda < thresholds).astype(bool)
    return cutoffs_from_mask(is_good_frame)


def cutoffs_from_dataset(ds):
    """
    Returns a da containing cutoff times (the same as the good_upto data_var)
    and with a coord called good_throughout
    -- question: is good_upto expensive enough to justify separating out
    calculation of good_throughout?

    Another question -- shall we prefer filtranda, good_upto, or check consistency?
    I think we should ignore filtranda if good_upto is available.
    """
    if 'good_upto' in ds.data_vars:
        res = ds.data_vars['good_upto']
        if 'good_throughout' in res.coords:
            if 'filtranda' in ds.data_vars:
                warning(
                    "data_vars 'filtranda' and 'good_upto' present in "
                    "the same dataset; ignoring 'filtranda'"
                )
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


####################
# Action functions #
####################

# All the action functions take a dataset
# They can use the functions above to get the info they need


def omit(ds):
    cutoffs = cutoffs_from_dataset(ds)
    good_throughout = cutoffs['good_throughout']
    selection = good_throughout.all('criterion')
    return sel_trajs(ds, selection)


def truncate(ds):
    mask = cum_mask_from_dataset(ds)
    # So does `mask_from_dataset` return stacked or unstacked?
    assert is_stacked(mask)
    selection = mask.all('criterion')
    return ds.isel(frame=selection)


def transect(ds, cutoff: float):
    if is_stacked(ds):
        ds = unstack_trajs(ds)
    ds = ds.loc[{'time': slice(float(cutoff))}]
    good_upto = cutoffs_from_dataset(ds)
    # NB. the second mask must be calculated after time-slicing.
    traj_selection = (good_upto >= cutoff).all('criterion') & ds.coords['is_frame'].all('time')
    return ds.isel[{'trajid': traj_selection}]


#
# Filtranda derivat
#


#########################
# Convenience functions #
#########################


def energy_filtranda(): ...


def sanity_check(
    frames,
    cut=False,
    *,
    units='eV',  # TODO: FIXME: Actually implement!
    etot_drift=0.2,
    etot_step=0.1,
    epot_step=0.7,
    ekin_step=0.7,
    hop_epot_step=1.0,
):
    settings = {
        k: locals()[k]
        for k in ['etot_drift', 'etot_step', 'epot_step', 'ekin_step', 'hop_epot_step']
    }
    frames = frames.assign(filtranda=energy_filtranda(frames, **settings))
    frames = get_cutoffs(frames)
    if not cut:
        return frames
    elif cut == 'truncate':
        return truncate(frames)
    elif cut == 'omit':
        return omit(frames)
    else:
        assert isinstance(cut, Number)
        return transect(frames, cut)


def filter_cleavages(
    atom_pair_identifiers: list,
    atom_pair_thresholds: dict,
): ...


###########################################
# Formats directly prerequisite for plots #
###########################################


# For check_thresholds aka plot_thresholds
def cum_max_quantiles(ds_or_da, quantiles=None):
    ...
    # Current implementation expects filtranda
    # or dataset containing filtranda


# For validity_populations
def validity_populations(ds, intersections=False): ...