import numpy as np
import xarray as xr


def bonds(): ...


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


# def last_time_where(mask_da):
#     return (
#         mask_da.groupby('trajid')
#         .map(
#             lambda traj: traj.coords['time'][-1]
#             if traj.all()
#             else traj.sel(frame=~traj).coords['time'][0]
#         )
#         .rename(trajid='trajid_')
#     )


def last_time_where(mask):
    before_first_false = ~((~mask).groupby('trajid').cumsum().astype(bool))
    upto_first_false = mask.coords['time'][before_first_false].groupby('trajid').last()
    fallback = (~before_first_false).groupby('trajid').all()
    if fallback.any():
        fallback = fallback.copy(
            data=np.full((len(fallback),), -1)  # -1 indicates first ts fails test
        )
        res = upto_first_false.combine_first(fallback)
    else:
        res = upto_first_false
    return res.rename(trajid='trajid_')


# def earliest_false(mask):
#     times_for_trajs_with_false = mask[~mask].coords['time'].groupby('trajid').min()
#     fallback = mask.coords['time'].groupby('trajid').max()
#     return xr.merge(
#         [times_for_trajs_with_false, fallback], join='right', compat='override'
#     ).fillna(fallback)


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
    expansion = cutoffs.sel(trajid_=frames.coords['trajid']).drop_vars('trajid_')
    mask = frames['time'] <= expansion
    return frames.sel(frame=mask)


def trajs_where(mask_da): ...