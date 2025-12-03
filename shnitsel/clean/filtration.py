from logging import warning

import numpy as np
import xarray as xr

from shnitsel.data.multi_indices import mdiff
from shnitsel.data.trajectory_format import Trajectory
from shnitsel.units.conversion import convert_energy
from .._contracts import needs

# link functions that have moved:
from ..geo.geom import get_bond_lengths as get_bond_lengths


@needs(data_vars={'energy', 'astate'})
def energy_filtranda(
    frames: Trajectory,
    *,
    etot_drift=0.2,
    etot_step=0.1,
    epot_step=0.7,
    ekin_step=0.7,
    hop_epot_step=1.0,
):
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

    da: xr.DataArray = np.abs(res.to_dataarray('criterion'))  # type: ignore # Result of numpy application to DataArray is a DA.
    thresholds = [default_thresholds[x] for x in da.coords['criterion'].data]
    return da.assign_coords(thresholds=('criterion', thresholds))


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
            "The mask argument should be trajectories, either stacked or unstacked"
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
