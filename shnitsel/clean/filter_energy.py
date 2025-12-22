from dataclasses import dataclass
from logging import warning
from numbers import Number
from typing import Literal, Sequence

import numpy as np
import xarray as xr

from shnitsel.data.multi_indices import mdiff
from shnitsel.clean.common import dispatch_cut
from shnitsel.clean.dispatch_plots import dispatch_plots
from shnitsel.units.conversion import convert_energy
from shnitsel.units.definitions import energy


@dataclass
class EnergyThresholds:
    etot_drift: float = 0.2
    etot_step: float = 0.1
    epot_step: float = 0.7
    ekin_step: float = 0.7
    hop_epot_step: float = 1.0
    energy_unit: str = energy.eV


def energy_filtranda(
    frames,
    *,
    energy_thresholds: EnergyThresholds | None = None,
):
    """Derive energetic filtration targets from an xr.Dataset

    Parameters
    ----------
    frames
        A xr.Dataset with ``astate``, ``energy``, and ideally ``e_kin`` variables
    energy_thresholds, optional
        Threshold for total, potential and kinetic energy of the system.
        Can specify thresholds for overall drift and individual time step changes.
        Can also specify thresholds for energy steps at hops.
        Unit should be specified as a member variable.
        If not provided will default to some reasonable default values as seen in `EnergyThresholds` definition.

    Returns
    -------
        An xr.DataArray of filtration targets stacked along the ``criterion`` dimension;
        criteria comprise epot_step and hop_epot_step, as well as
        etot_drift, etot_step and ekin_step if the input contains an e_kin variable
    """

    if energy_thresholds is None:
        energy_thresholds = EnergyThresholds()

    res = xr.Dataset()
    is_hop = mdiff(frames['astate']) != 0
    # TODO: FIXME: Shouldn't we drop coords instead?
    e_pot_active = frames.energy.sel(state=frames.astate).drop_vars('state')
    e_pot_active.attrs['units'] = frames['energy'].attrs['units']
    # e_pot_active = convert_energy(e_pot_active, to=units)

    res['epot_step'] = mdiff(e_pot_active).where(~is_hop, 0)
    res['hop_epot_step'] = mdiff(e_pot_active).where(is_hop, 0)

    if 'e_kin' in frames.data_vars:
        e_kin = frames['e_kin']
        e_kin.attrs['units'] = frames['e_kin'].attrs['units']
        # e_kin = convert_energy(e_kin, to=units)

        e_tot = e_pot_active + e_kin
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
    energy_thresholds: EnergyThresholds | None = None,
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
    energy_thresholds, optional
        Threshold for total, potential and kinetic energy of the system.
        Can specify thresholds for overall drift and individual time step changes.
        Can also specify thresholds for energy steps at hops.
        Unit should be specified as a member variable.
        If not provided will default to some reasonable default values as seen in `EnergyThresholds` definition.
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
    if energy_thresholds is None:
        energy_thresholds = EnergyThresholds()
    filtranda = energy_filtranda(frames, energy_thresholds=energy_thresholds)
    dispatch_plots(filtranda, plot_thresholds, plot_populations)
    filtered_frames = frames.drop_dims(['criterion'], errors='ignore').assign(
        filtranda=filtranda
    )
    return dispatch_cut(filtered_frames, cut)
