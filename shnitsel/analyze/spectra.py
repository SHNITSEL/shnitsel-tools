from itertools import product
from typing import Hashable, TypeAlias

import numpy as np
import xarray as xr

from .generic import subtract_combinations
from .._contracts import needs
from ..units import convert_energy

DimName: TypeAlias = Hashable


def _get_fosc(energy, dip_trans):
    return 2 / 3 * energy * dip_trans**2


def get_fosc(energy, dip_trans):
    """Function to obtain a dataarray containing the oscillator strength as a dataarray.

    Args:
        energy (DataArray): _description_
        dip_trans (DataArray): _description_

    Returns:
        DataArray: The resulting datarray of oscillator strength f_osc
    """
    if 'state' in energy.dims:
        assert 'statecomb' not in energy.dims
        energy = subtract_combinations(energy, 'state')

    da = _get_fosc(convert_energy(energy, to='hartree'), dip_trans)
    da.name = 'fosc'
    da.attrs['long_name'] = r"$f_{\mathrm{osc}}$"
    return da


# TODO: deprecate (made redundant by DerivedProperties)
@needs(data_vars={'energy', 'dip_trans'})
def assign_fosc(ds: xr.Dataset) -> xr.Dataset:
    """Function to calculate oscillator strength fosc and create a new dataset with this variable assigned.

    Args:
        ds (xr.Dataset): Dataset from which to calculate fosc

    Returns:
        xr.Dataset: Dataset with the member variable fosc set
    """
    da = get_fosc(ds['energy'], ds['dip_trans'])
    return ds.assign(fosc=da)


@needs(data_vars={'energy', 'fosc'})
def broaden_gauss(
    E: xr.DataArray,
    fosc: xr.DataArray,
    agg_dim: DimName = 'frame',
    *,
    width: float = 0.5,
    nsamples: int = 1000,
    xmin: float = 0,
    xmax: float | None = None,
) -> xr.DataArray:
    r"""
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

    stdev = width / 2

    def g(x):
        nonlocal stdev
        return 1 / (np.sqrt(2 * np.pi) * stdev) * np.exp(-(x**2) / (2 * stdev**2))

    xname = getattr(E, 'name', 'energy') or 'xdim'
    yname = getattr(fosc, 'name', 'fosc') or 'ydim'

    if xmax is None:
        # broadening could visibly overshoot the former maximum by 3 standard deviations
        xmax = E.max().item() * (1 + 1.5 * width)
    xs = np.linspace(0, xmax, num=nsamples)
    Espace = xr.DataArray(xs, dims=[xname], attrs=E.attrs)
    res: xr.DataArray = (g(Espace - E) * fosc).mean(dim=agg_dim)
    res.name = yname
    res.attrs = fosc.attrs
    for cname, coord in res.coords.items():
        if cname in fosc.coords:
            coord.attrs = fosc.coords[cname].attrs
    return res.assign_coords({xname: Espace})


def ds_broaden_gauss(
    ds: xr.Dataset, width: float = 0.5, nsamples: int = 1000, xmax: float | None = None
) -> xr.DataArray:
    return broaden_gauss(
        ds['energy'], ds['fosc'], width=width, nsamples=nsamples, xmax=None
    )


def get_spectrum(data, t, sc, cutoff=0.01):
    # following required because `meathod='nearest'` doesn't work for MultiIndex
    if t not in data.coords['time']:
        times = np.unique(data.coords['time'])
        diffs = np.abs(times - t)
        t = times[np.argmin(diffs)]
    data = data.sel(time=t, statecomb=sc)
    res = broaden_gauss(data.energy, data.fosc, agg_dim='trajid')

    max_ = res.max().item()
    non_negligible = res.where(res > cutoff * max_, drop=True).energy
    if len(non_negligible) == 0:
        return res.sel(energy=non_negligible)
    return res.sel(energy=slice(non_negligible.min(), non_negligible.max()))


def calc_spectra(spectral, times=None, cutoff=0.01):
    """Returns a `dict` of DataArrays indexed by `(time, statecomb)` tuples."""
    if times is None:
        times = [0, 10, 20, 30]
    return {
        (t, sc): get_spectrum(spectral, t, sc, cutoff=cutoff)
        for t, sc in product(times, spectral.statecomb.values)
    }


def get_sgroups(spectra):
    ground, excited = {}, {}
    for (t, sc), v in spectra.items():
        if sc == '$S_2 - S_1$':
            excited[t, sc] = v
        else:
            ground[t, sc] = v

    sgroups = (ground, excited)
    return sgroups


def sep_ground_excited_spectra(spectra, excited_transitions=None):
    if excited_transitions is None:
        excited_transitions = {'$S_2 - S_1$'}

    ground, excited = {}, {}

    for (t, sc), v in spectra.items():
        if sc in excited_transitions:
            excited[t, sc] = v
        else:
            ground[t, sc] = v

    return ground, excited

@needs(data_vars={'energy', 'fosc'}, coords={'frame', 'trajid'})
def spectra_all_times(inter_state: xr.Dataset):
    assert isinstance(inter_state, xr.Dataset)
    if 'energy' not in inter_state.data_vars:
        raise ValueError("Missing required variable 'energy'")
    if 'fosc' not in inter_state.data_vars:
        raise ValueError("Missing required variable 'fosc'")
    assert (
        'frame' in inter_state and 'trajid' in inter_state
    ), "Missing required dimensions"

    data = inter_state.unstack('frame')
    return broaden_gauss(data.energy, data.fosc, agg_dim='trajid')