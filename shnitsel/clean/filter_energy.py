from logging import warning
from numbers import Number
from typing import Literal

import numpy as np
import xarray as xr

from shnitsel.data.multi_indices import mdiff
from shnitsel.clean.common import dispatch_cut
from shnitsel.units.conversion import convert_energy

_default_energy_thresholds_eV = {
    'etot_drift': 0.2,
    'etot_step': 0.1,
    'epot_step': 0.7,
    'ekin_step': 0.7,
    'hop_epot_step': 1.0,
}


def energy_filtranda(
    frames,
    *,
    etot_drift: float | None = None,
    etot_step: float | None = None,
    epot_step: float | None = None,
    ekin_step: float | None = None,
    hop_epot_step: float | None = None,
    units='eV',
):
    res = xr.Dataset()
    is_hop = mdiff(frames['astate']) != 0
    e_pot = frames.energy.sel(state=frames.astate).drop_vars('state')
    e_pot.attrs['units'] = frames['energy'].attrs['units']
    e_pot = convert_energy(e_pot, to=units)

    res['epot_step'] = mdiff(e_pot).where(~is_hop, 0)
    res['hop_epot_step'] = mdiff(e_pot).where(is_hop, 0)

    if 'e_kin' in frames.data_vars:
        e_kin = frames['e_kin']
        e_kin.attrs['units'] = frames['e_kin'].attrs['units']
        e_kin = convert_energy(e_kin, to=units)

        e_tot = e_pot + e_kin
        res['etot_drift'] = e_tot.groupby('trajid').map(
            lambda traj: abs(traj - traj.item(0))
        )
        res['ekin_step'] = mdiff(e_kin).where(~is_hop, 0)
        res['etot_step'] = mdiff(e_tot)
    else:
        e_kin = None
        warning("data does not contain kinetic energy variable ('e_kin')")

    da = np.abs(res.to_dataarray('criterion')).assign_attrs(units=units)

    # Make threshold coordinates

    def dict_to_thresholds(d: dict, units: str) -> xr.DataArray:
        criteria = da.coords['criterion'].data
        data = [d[c] for c in criteria]
        res = xr.DataArray(
            list(data), coords={'criterion': criteria}, attrs={'units': units}
        )
        return res.astype(float)

    default_thresholds = dict_to_thresholds(_default_energy_thresholds_eV, units='eV')
    default_thresholds = convert_energy(default_thresholds, to=units)
    user_thresholds = dict_to_thresholds(locals(), units=units)
    thresholds = user_thresholds.where(~np.isnan(user_thresholds), default_thresholds)

    da = da.assign_coords(thresholds=thresholds)
    return da


def sanity_check(
    frames,
    cut: Literal['truncate', 'omit', False] | Number = 'truncate',
    *,
    units='eV',
    etot_drift: float = np.nan,
    etot_step: float = np.nan,
    epot_step: float = np.nan,
    ekin_step: float = np.nan,
    hop_epot_step: float = np.nan,
):
    settings = {
        'etot_drift': etot_drift,
        'etot_step': etot_step,
        'epot_step': epot_step,
        'ekin_step': ekin_step,
        'hop_epot_step': hop_epot_step,
        'units': units,
    }
    frames = frames.assign(filtranda=energy_filtranda(frames, **settings))
    return dispatch_cut(frames, cut)