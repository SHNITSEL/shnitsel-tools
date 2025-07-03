from itertools import product
import numpy as np
from . import postprocess as P


def get_spectrum(data, t, sc, cutoff=0.01):
    # following required because `meathod='nearest'` doesn't work for MultiIndex
    if t not in data.coords['time']:
        times = np.unique(data.coords['time'])
        diffs = np.abs(times - t)
        t = times[np.argmin(diffs)]
    data = data.sel(time=t, statecomb=sc)
    res = P.broaden_gauss(data.energy, data.fosc, agg_dim='trajid')

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
