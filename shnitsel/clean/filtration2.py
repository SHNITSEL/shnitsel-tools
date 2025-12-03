from logging import warning
from numbers import Number

import numpy as np
import xarray as xr

from shnitsel.core.typedefs import Frames
from shnitsel.data.multi_indices import sel_trajs, stack_trajs, unstack_trajs
from shnitsel.data.trajectory_format import Trajectory


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
    if is_stacked(filtranda):
        restack = True
        filtranda = unstack_trajs(filtranda)
    else:
        restack = False

    res = (filtranda < filtranda.coords['thresholds']).cumprod('time').astype(bool)

    # Is unstack-restack meaningfully less efficient than groupby?
    if restack:
        return stack_trajs(res)
    else:
        return res


# def cum_mask_from_dataset_cutoffs(ds):
#     return (ds.coords['time'] <= ds.coords['good_upto']).cumprod('time').astype(bool)


def cum_mask_from_dataset(ds):
    # dataset might contain filtranda, the mask we're after, or the cutoffs
    # it doesn't matter whether the mask is cumulative, because our
    # accumulation is idempotent

    # two versions currently, for stacked and for unstacked
    # use `filtranda` and `thresholds` if available
    # otherwise use `good_upto`
    if 'is_good_frame' in ds:
        return ds['is_good_frame'].cumprod('time').astype(bool)
    elif 'good_upto' in ds.data_vars:
        return (
            (ds.coords['time'] <= ds.coords['good_upto']).cumprod('time').astype(bool)
        )
    elif 'filtranda' in ds.data_vars:
        return cum_mask_from_filtranda(ds)


def true_upto(mask, dim):
    if is_stacked(mask):
        mask = unstack_trajs(mask).fillna(False)
    # shifted_coord = xr.DataArray(np.concat([[-1], mask.coords[dim].data]), dims=dim)
    shifted_coord = np.concat([[-1], mask.coords[dim].data])
    indexer = mask.cumprod(dim).sum(dim).astype(int)
    # res = shifted_coord[{dim: indexer}]
    # res = np.choose(indexer.data, shifted_coord)
    res = np.take(shifted_coord, indexer.data)
    return indexer.copy(data=res)


def cutoffs_from_mask(mask):
    if is_stacked(mask):
        mask = unstack_trajs(mask)
    good_upto = true_upto(mask, 'time')
    good_throughout = mask.all('time')
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


def assign_cutoffs(ds): ...


####################
# Action functions #
####################

# All the action functions take a dataset
# They can use the functions above to get the info they need


def omit(ds:Trajectory):
    cutoffs = cutoffs_from_dataset(ds)
    good_throughout = cutoffs['good_throughout']
    selection = good_throughout.all('criterion')
    return sel_trajs(ds, selection)


def truncate(ds:Trajectory):
    mask = cum_mask_from_dataset(ds)
    # So does `mask_from_dataset` return stacked or unstacked?
    assert is_stacked(mask)
    selection = mask.all('criterion')
    return ds.isel(frame=selection)


def transect(ds: Trajectory|Frames, cutoff: float):
    if is_stacked(ds):
        ds = unstack_trajs(ds)
    ds = ds.loc[{'time': slice(float(cutoff))}]
    good_upto = cutoffs_from_dataset(ds)
    # NB. the second mask must be calculated after time-slicing.
    traj_selection = (good_upto >= cutoff).all('criterion') & ds.coords['is_frame'].all(
        'time'
    )
    return ds.isel[{'trajid': traj_selection}]


#
# Filtranda derivat
#


#########################
# Convenience functions #
#########################


def energy_filtranda(
    frames,
    *,
    etot_drift=0.2,
    etot_step=0.1,
    epot_step=0.7,
    ekin_step=0.7,
    hop_epot_step=1.0,
):
    from logging import warning
    from shnitsel.data.multi_indices import mdiff
    from shnitsel.units.conversion import convert_energy

    default_thresholds = {
        'etot_drift': etot_drift,
        'etot_step': etot_step,
        'epot_step': epot_step,
        'ekin_step': ekin_step,
        'hop_epot_step': hop_epot_step,
    }

    res = xr.Dataset()
    is_hop = mdiff(frames['astate']) != 0
    e_pot = frames.energy.sel(state=frames.astate).drop_vars('state')
    e_pot.attrs['units'] = frames['energy'].attrs.get('units', 'unknown')
    e_pot = convert_energy(e_pot, to='eV')

    res['epot_step'] = mdiff(e_pot).where(~is_hop, 0)
    res['hop_epot_step'] = mdiff(e_pot).where(is_hop, 0)

    if 'e_kin' in frames.data_vars:
        e_kin = frames['e_kin']
        e_kin.attrs['units'] = frames['e_kin'].attrs.get('units', 'unknown')
        e_kin = convert_energy(e_kin, to='eV')

        e_tot = e_pot + e_kin
        res['etot_drift'] = e_tot.groupby('trajid').map(
            lambda traj: abs(traj - traj.item(0))
        )
        res['ekin_step'] = mdiff(e_kin).where(~is_hop, 0)
        res['etot_step'] = mdiff(e_tot)
    else:
        e_kin = None
        warning("data does not contain kinetic energy variable ('e_kin')")

    da: xr.DataArray = np.abs(res.to_dataarray('criterion'))  # type: ignore # numpy on DataArray yields DataArray.
    thresholds = [default_thresholds[x] for x in da.coords['criterion'].data]
    return da.assign_coords(thresholds=('criterion', thresholds))


def sanity_check(
    frames,
    cut='truncate',
    *,
    units='eV',  # TODO: FIXME: Actually implement!
    etot_drift: float = 0.2,
    etot_step: float = 0.1,
    epot_step: float = 0.7,
    ekin_step: float = 0.7,
    hop_epot_step: float = 1.0,
):
    settings = {
        k: locals()[k]
        for k in ['etot_drift', 'etot_step', 'epot_step', 'ekin_step', 'hop_epot_step']
    }
    # filtranda = energy_filtranda(frames, **settings)
    frames = frames.assign(filtranda=energy_filtranda(frames, **settings))
    # frames = get_cutoffs(frames)
    if not cut:
        return frames.assign(good_upto=cutoffs_from_dataset(frames))
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
