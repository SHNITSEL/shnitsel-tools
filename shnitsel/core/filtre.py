import numpy as np
import xarray as xr

# link functions that have moved:
from .geom import get_bond_lengths as get_bond_lengths


def energy_filtranda(frames: xr.Dataset) -> xr.Dataset:
    res = frames[['e_kin']]
    res['e_pot'] = frames.energy.sel(state=frames.astate).drop_vars('state')
    res['e_tot'] = res['e_pot'] + res['e_kin']

    res['etot_drift'] = (
        res['e_tot'].groupby('trajid').map(lambda traj: abs(traj - traj.isel(frame=0)))
    )
    res['ekin_step'] = res['e_kin'].sh.sudi()
    res['epot_step'] = res['e_pot'].sh.sudi()
    res['etot_step'] = res['e_tot'].sh.sudi()
    res['is_hop'] = frames['astate'].sh.sudi() != 0

    return res


def last_time_where(mask):
    mask = mask.unstack('frame', fill_value=False).transpose('trajid', 'time', ...)
    idxs = np.logical_not((~mask.values).cumsum(axis=1)).sum(axis=1)
    times = np.concat([[-1], mask.time.values])
    return mask[:, 0].copy(data=times[idxs]).drop_vars('time').rename('time')


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


def truncate(frames, cutoffs):
    if 'trajid_' not in cutoffs.coords and 'trajid' in cutoffs.coords:
        cutoffs = cutoffs.rename(trajid='trajid_')
    expansion = cutoffs.sel(trajid_=frames.coords['trajid']).drop_vars('trajid_')
    mask = frames['time'] <= expansion
    return frames.sel(frame=mask)


def trajs_where(mask_da): ...