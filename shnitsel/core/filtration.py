from logging import warning

import numpy as np
import xarray as xr

from shnitsel.core.midx import mdiff
from .._contracts import needs

# link functions that have moved:
from .geom import get_bond_lengths as get_bond_lengths

@needs(data_vars={'energy', 'astate'})
def energy_filtranda(frames: xr.Dataset) -> xr.Dataset:
    if 'e_kin' in frames.data_vars:
        has_e_kin = True
        res = frames[['e_kin']]
    else:
        has_e_kin = False
        res = xr.Dataset()
        warning("data does not contain kinetic energy variable ('e_kin')")

    res['e_pot'] = frames.energy.sel(state=frames.astate).drop_vars('state')
    if has_e_kin:
        res['e_tot'] = res['e_pot'] + res['e_kin']

        res['etot_drift'] = (
            res['e_tot'].groupby('trajid').map(lambda traj: abs(traj - traj.isel(frame=0)))
        )
        res['ekin_step'] = mdiff(res['e_kin'])
        res['etot_step'] = mdiff(res['e_tot'])

    res['epot_step'] = mdiff(res['e_pot'])
    res['is_hop'] = mdiff(frames['astate']) != 0

    return res


@needs(dims={'frame'}, coords={'trajid', 'time'})
def last_time_where(mask):
    # If tree:
    if isinstance(mask, xr.DataTree):
        raise NotImplementedError("This function does not yet support DataTree")
    # If stacked:
    elif 'frame' in mask.dims and {'trajid', 'time'}.issubset(mask.coords):
        mask = mask.unstack('frame', fill_value=False).transpose('trajid', 'time', ...)
    # Otherwise, had better be unstacked
    elif {'trajid', 'time'}.issubset(mask.dims):
        mask = mask.transpose('trajid', 'time', ...)
    else:
        raise ValueError(
            "The mask argument should be trajectories, "
            "either stacked or unstacked"
        )
    idxs = np.logical_not((~mask.values).cumsum(axis=1)).sum(axis=1)
    times = np.concat([[-1], mask.time.values])
    return mask[:, 0].copy(data=times[idxs]).drop_vars('time').rename('time')


@needs(dims={'frame'}, coords={'trajid', 'time'})
def get_cutoffs(masks_ds):
    ds = masks_ds.map(last_time_where)
    names = list(ds.data_vars)
    ds['original'] = (
        masks_ds.coords['time'].groupby('trajid').max().rename(trajid='trajid_')
    )
    # Put 'original' first, so that min() chooses it in cases of ambiguity
    ds = ds[['original'] + names]
    names = list(ds.data_vars)

    da = ds.to_dataarray('cutoff')
    ds['earliest'] = da.min('cutoff')
    ds['reason'] = da.argmin('cutoff')
    ds.attrs['reasons'] = names
    return ds

@needs(dims={'frame'}, coords={'trajid', 'time'})
def truncate(frames, cutoffs):
    if 'trajid_' not in cutoffs.coords and 'trajid' in cutoffs.coords:
        cutoffs = cutoffs.rename(trajid='trajid_')
    expansion = cutoffs.sel(trajid_=frames.coords['trajid']).drop_vars('trajid_')
    mask = frames['time'] <= expansion
    return frames.sel(frame=mask)


def trajs_where(mask_da): ...