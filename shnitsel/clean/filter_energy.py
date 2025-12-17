from logging import warning
from numbers import Number
from typing import Literal, Sequence

import numpy as np
import xarray as xr

from shnitsel.data.multi_indices import mdiff
from shnitsel.clean.common import dispatch_cut
from shnitsel.clean.dispatch_plots import dispatch_plots
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
    """Derive energetic filtration targets from an xr.Dataset

    Parameters
    ----------
    frames
        A xr.Dataset with ``astate``, ``energy``, and ideally ``e_kin`` variables
    etot_drift, optional
        Threshold for drift of total energy over an entire trajectory, by default 0.2 eV
    etot_step, optional
        Threshold for difference in total energy from one frame to the next, ignoring hops
        , by default 0.1 eV
    epot_step, optional
        Threshold for difference in potential energy from one frame to the next, ignoring hops, by default 0.7 eV
    ekin_step, optional
        Threshold for difference in kinetic energy from one frame to the next, ignoring hops, by default 0.7 eV
    hop_epot_step, optional
        Threshold for difference in potential energy across hops, by default 1.0 eV
    units, optional
        Units in which custom thresholds are given, and to which defaults and data will be converted, by default 'eV'

    Returns
    -------
        An xr.DataArray of filtration targets stacked along the ``criterion`` dimension;
        criteria comprise epot_step and hop_epot_step, as well as
        etot_drift, etot_step and ekin_step if the input contains an e_kin variable
    """
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
    plot_thresholds: bool | Sequence[float] = False,
    plot_populations: bool | Literal['independent', 'intersections'] = False,
):
    """Filter trajectories according to energy to exclude unphysical (insane) behaviour

    Parameters
    ----------
    frames
        A xr.Dataset with ``astate``, ``energy``, and ideally ``e_kin`` variables
    cut, optional
        Specifies the manner in which to remove data;

            - if 'omit', drop trajectories unless all frames meet criteria (:py:func:`shnitsel.clean.omit`)
            - if 'truncate', cut each trajectory off just before the first frame that doesn't meet criteria
              (:py:func:`shnitsel.clean.truncate`)
            - if a number, interpret this number as a time, and cut all trajectories off at this time,
              discarding those which violate criteria before reaching the given limit,
              (:py:func:`shnitsel.clean.transect`)
            - if ``False``, merely annotate the data;
        see :py:func:`shnitsel.clean.dispatch_cut`.
    units, optional
        Units in which custom thresholds are given, and to which defaults and data will be converted, by default 'eV'
    etot_drift, optional
        Threshold for drift of total energy over an entire trajectory, by default 0.2 eV
    etot_step, optional
        Threshold for difference in total energy from one frame to the next, ignoring hops
        , by default 0.1 eV
    epot_step, optional
        Threshold for difference in potential energy from one frame to the next, ignoring hops, by default 0.7 eV
    ekin_step, optional
        Threshold for difference in kinetic energy from one frame to the next, ignoring hops, by default 0.7 eV
    hop_epot_step, optional
        Threshold for difference in potential energy across hops, by default 1.0 eV
    plot_thresholds
        See :py:func:`shnitsel.vis.plot.filtration.check_thresholds`.

        - If ``True``, will plot using ``check_thresholds`` with
        default quantiles
        - If a ``Sequence``, will plot using ``check_thresholds``
        with specified quantiles
        - If ``False``, will not plot threshold plot
    plot_populations
        See :py:func:`shnitsel.vis.plot.filtration.validity_populations`.

        - If ``True`` or ``'intersections'``, will plot populations of
        trajectories satisfying intersecting conditions
        - If ``'independent'``, will plot populations of
        trajectories satisfying conditions taken independently
        - If ``False``, will not plot populations plot
    Returns
    -------
        The sanitized xr.Dataset

    Notes
    -----
    The resulting object has a ``filtranda`` data_var, representing the values by which the data were filtered.
    If the input has a ``filtranda`` data_var, it is overwritten.
    """
    settings = {
        'etot_drift': etot_drift,
        'etot_step': etot_step,
        'epot_step': epot_step,
        'ekin_step': ekin_step,
        'hop_epot_step': hop_epot_step,
        'units': units,
    }
    filtranda = energy_filtranda(frames, **settings)
    dispatch_plots(filtranda, plot_thresholds, plot_populations)
    frames = frames.drop_dims(['criterion'], errors='ignore').assign(
        filtranda=filtranda
    )
    return dispatch_cut(frames, cut)