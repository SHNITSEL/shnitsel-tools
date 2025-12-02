from itertools import product
import logging
from typing import Iterable

import numpy as np
import xarray as xr

from shnitsel.units.definitions import energy, dipole

from .typedefs import InterState, DimName, SpectraDictType

from .generic import keep_norming, subtract_combinations
from .._contracts import needs
from ..units import convert_energy, convert_dipole


def _get_fosc(energy_interstate: xr.DataArray, dip_trans_norm: xr.DataArray) -> xr.DataArray:
    """Internal function to actually calculate the oscillator frequency for energies and transition dipoles.

    Args:
        energy_interstate (DataArray): The Array of Energies in the system.
        dip_trans_norm (DataArray): The array of associated norm of transition dipoles in the system.


    Returns:
        DataArray: The resulting oscillation frequency (f_osc) array.
    """

    return 2 / 3 * energy_interstate * dip_trans_norm**2


def get_fosc(
    energy_per_or_interstate: xr.DataArray, dip_trans_norm: xr.DataArray
) -> xr.DataArray:
    """Function to obtain a dataarray containing the oscillator strength as a dataarray.

    Args:
        energy_interstate (DataArray): The array of per- or inter-state energies in the system.
            If provided as a per-state energy, inter-state barriers will automatically be calculated.
        dip_trans_norm (DataArray): The array of associated transition dipoles in the system with their norm calculated across the direction dimension.

    Returns:
        DataArray: The resulting datarray of oscillator strength f_osc
    """
    if 'state' in energy_per_or_interstate.dims:
        assert 'statecomb' not in energy_per_or_interstate.dims
        energy_interstate = subtract_combinations(energy_per_or_interstate, 'state')
    else:
        energy_interstate = energy_per_or_interstate

    assert (
        energy_interstate.values.shape == dip_trans_norm.values.shape
    ), f"Energy and dip_trans do not have the same shapes: {energy_interstate.values.shape} <-> {dip_trans_norm.values.shape}"

    da = _get_fosc(convert_energy(energy_interstate, to=energy.Hartree), convert_dipole(dip_trans_norm, to=dipole.au))
    da.name = 'fosc'
    da.attrs.update(
        {
            "long_name": r"$f_{\mathrm{osc}}$",
            "description": "derived from 'energy_interstate' and 'dip_trans' variables",
        }
    )
    return da


# TODO: deprecate (made redundant by DerivedProperties)
@needs(data_vars={'energy_interstate', 'dip_trans'}, coords={'statecomb'})
def assign_fosc(ds: xr.Dataset) -> xr.Dataset:
    """Function to calculate oscillator strength fosc and create a new dataset with this variable assigned.

    Args:
        ds (xr.Dataset): Dataset from which to calculate fosc

    Returns:
        xr.Dataset: Dataset with the member variable fosc set
    """
    if 'dip_trans_norm' not in ds:
        ds['dip_trans_norm'] = keep_norming(ds['dip_trans'])

    da = get_fosc(ds['energy_interstate'], ds['dip_trans_norm'])
    res = ds.assign(fosc=da)
    return res


@needs(data_vars={'energy', 'fosc'})
def broaden_gauss(
    E: xr.DataArray,
    fosc: xr.DataArray,
    agg_dim: DimName = 'frame',
    *,
    width: float = 0.5, # in eV
    nsamples: int = 1000,
    xmin: float = 0,
    xmax: float | None = None,
) -> xr.DataArray:
    r"""
    Applies a gaussian smoothing kernel to the fosc data.

    Aggregation is performed along the `agg_dim` dimension.

    Parameters
    ----------
    E
        values used for the x-axis, presumably $E_i$
    fosc
        values used for the y-axis, presumably $f_\mathrm{osc}$
    agg_dim, optional
        dimension along which to aggregate the many Gaussian distributions,
        by default 'frame'
    width, optional
        the width (i.e. 2 standard deviations) of the Gaussian distributions
        used, by default 0.001
    nsamples, optional
        number of evenly spaced x-values over which to sample the distribution,
        by default 1000
    xmax, optional
        the maximum x-value, by default 3 standard deviations
        beyond the pre-broadened maximum
    """
    assert (
        agg_dim in E.sizes
    ), f"E does not have required dimension {agg_dim} for aggregation"

    stdev = width / 2

    E_eV = convert_energy(E, to=energy.eV)

    def g(x):
        nonlocal stdev
        return 1 / (np.sqrt(2 * np.pi) * stdev) * np.exp(-(x**2) / (2 * stdev**2))

    xname = getattr(E_eV, 'name', 'energy') or 'xdim'
    yname = getattr(fosc, 'name', 'fosc') or 'ydim'

    if xmax is None:
        # TODO: FIXME: The calculation does not fit the statement of the comment above.
        # broadening could visibly overshoot the former maximum by 3 standard deviations
        xmax = E_eV.max().item() * (1 + 1.5 * width)

        assert xmax is not None, "Could not calculate maximum of the provided energy"

    xs = np.linspace(0, xmax, num=nsamples)
    Espace = xr.DataArray(xs, dims=[xname], attrs=E_eV.attrs)
    res: xr.DataArray = (g(Espace - E_eV) * fosc).mean(dim=agg_dim)
    # print(res)
    res.name = yname
    res.attrs = fosc.attrs
    for cname, coord in res.coords.items():
        if cname in fosc.coords:
            coord.attrs = fosc.coords[cname].attrs
    return res.assign_coords({xname: Espace})


@needs(data_vars={'energy', 'fosc'})
def ds_broaden_gauss(
    ds: xr.Dataset, width: float = 0.5, nsamples: int = 1000, xmax: float | None = None
) -> xr.DataArray:
    return broaden_gauss(
        ds['energy'], ds['fosc'], width=width, nsamples=nsamples, xmax=None
    )


@needs(data_vars={'energy_interstate', 'fosc'}, coords={"statecomb", "time"})
def get_spectrum(
    data: InterState, t: float, sc: tuple[int, int], rel_cutoff: float = 0.01
) -> xr.DataArray:
    """Function to calculate a gaussian-smoothed spectrum of an interstate dataset

    Args:
        data (InterState): An InterState dataset with fosc and energy data
        t (float): The time at which to evaluate the spectrum
        sc (tuple[int, int]): State combination identifier. Possibly an index or a tuple (from, to) of states.
        rel_cutoff (float, optional): Relative cutoff threshold. Values below the max of the resulting spectrum times this scale will be ignored. Defaults to 0.01.

    Returns:
        xr.DataArray: The Gauss-broadened spectrum of the provided `data` system.
            If broadening across trajectories could not be performed, just returns the fosc array.
    """
    # following required because `method='nearest'` doesn't work for MultiIndex
    if t not in data.coords['time']:
        times = np.unique(data.coords['time'])
        diffs = np.abs(times - t)
        t = times[np.argmin(diffs)]

    # Only take one timestep and one state combination
    data = data.sel(time=t, statecomb=sc)

    # Figure out how the trajectory is indexed across multiple trajectories or whether a single trajectory is provided
    trajid_dim = None
    if "trajid_" in data.energy_interstate.sizes:
        trajid_dim = "trajid_"
    elif "trajid" in data.energy_interstate.sizes:
        trajid_dim = "trajid"

    if trajid_dim is not None:
        res: xr.DataArray = broaden_gauss(
            data.energy_interstate, data.fosc, agg_dim=trajid_dim
        )
        max_ = res.max().item()
        non_negligible = res.where(res > rel_cutoff * max_, drop=True).energy_interstate
        if len(non_negligible) == 0:
            return res.sel(energy_interstate=non_negligible)
        return res.sel(
            energy_interstate=slice(non_negligible.min(), non_negligible.max())
        )
    else:
        logging.warning(
            "A single trajectory was provided. No gaussian smoothing across multiple trajectories could be performed."
        )
        raise NotImplementedError(
            "Cannot perform spectra calculation for single trajectory."
        )


@needs(data_vars={'energy', 'fosc'}, coords={"statecomb", "time"})
def calc_spectra(
    interstate: InterState,
    times: Iterable[float] | None = None,
    rel_cutoff: float = 0.01,
) -> SpectraDictType:
    """Function to

    Args:
        interstate (InterState): An InterState transformed Dataset.
        times (Iterable[float]|None, optional): The times at which the spectrum should be calculated. Defaults to None. If None, will be initialized as [0,10,20,30]
        rel_cutoff (float, optional): Factor for the cutoff of broadened/smoothened spectrum relative to maximum. Defaults to 0.01.

    Returns:
        SpectraDictType: Returns a `dict` of DataArrays indexed by `(time, statecomb)` tuples.
    """
    if times is None:
        times = [0, 10, 20, 30]

    sc_values: Iterable[tuple[int, int]] = interstate.statecomb.values

    res: SpectraDictType = {
        (t, sc): get_spectrum(interstate, t, sc, rel_cutoff=rel_cutoff)
        for t, sc in product(times, sc_values)
    }
    return res


def get_spectra_groups(
    spectra: SpectraDictType,
) -> tuple[
    SpectraDictType,
    SpectraDictType,
]:
    """Group spectra results into spectra involving the ground state or only excited states.

    Args:
        spectra (SpectraDictType): The Spectral calculation results, e.g. from `calc_spectra()`. Indexed by (timestep, state_combination) and yielding the associated spectrum.

    Returns:
        SpectraDictType: First the spectra involving the ground state
        SpectraDictType: Second the spectra involving only excited states.
    """

    ground, excited = {}, {}
    for (t, (sc_from, sc_to)), v in spectra.items():
        if sc_from > 1 and sc_to > 1:
            excited[t, (sc_from, sc_to)] = v
        else:
            ground[t, (sc_from, sc_to)] = v

    sgroups = (ground, excited)
    return sgroups


# TODO: FIXME: This looks like the more general version of get_spectra_groups?
def sep_ground_excited_spectra(
    spectra: SpectraDictType, excited_transitions: set[tuple[int, int]] | None = None
) -> tuple[
    SpectraDictType,
    SpectraDictType,
]:
    """Function to split Spectra into two groups.

    Can specify which state combinations should be grouped into `excited` transitions.
    If not provided a `excited_transitions` set, will assume all transitions not involving the ground state to be 'excited'.

    Args:
        spectra (SpectraDictType): The spectral results, e.g. from `calc_spectra()`
        excited_transitions (set[tuple[int, int]] | None, optional): An optional set specifying all state transitions that should be filtered as 'excited. Defaults to None.

    Returns:
        SpectraDictType: First the spectra not involving excited state transitions
        SpectraDictType: Second the spectra involving only excited state transitions
    """
    # This is too specific and does not fit the integer data model anymore.
    # if excited_transitions is None:
    #     excited_transitions = {'$S_2 - S_1$'}

    ground, excited = {}, {}

    for (t, sc), v in spectra.items():
        if (excited_transitions is not None and sc in excited_transitions) or (
            excited_transitions is None and sc[0] > 1 and sc[1] > 1
        ):
            excited[t, sc] = v
        else:
            ground[t, sc] = v

    return ground, excited


@needs(data_vars={'energy', 'fosc'}, coords={'frame', 'trajid_'})
def spectra_all_times(inter_state: xr.Dataset) -> xr.DataArray:
    """Function to calculate the spectra at all times.

    Does not return a dict with only the relevant (t,sc) combinations as above but instead a full
    xr.DataArray with a time dimension that has spectrum data for all times within the dataset averaged across trajectories.

    Args:
        inter_state (xr.Dataset): The InterState transformed Dataset.

    Raises:
        ValueError: If required variables or dimensions are missing

    Returns:
        xr.DataArray: The resulting spectra across all times.
    """
    assert isinstance(inter_state, xr.Dataset)
    if 'energy' not in inter_state.data_vars:
        raise ValueError("Missing required variable 'energy'")
    if 'fosc' not in inter_state.data_vars:
        raise ValueError("Missing required variable 'fosc'")
    assert (
        'frame' in inter_state and 'trajid_' in inter_state
    ), "Missing required dimensions"

    data = inter_state.unstack('frame')
    return broaden_gauss(data.energy, data.fosc, agg_dim='trajid_')
