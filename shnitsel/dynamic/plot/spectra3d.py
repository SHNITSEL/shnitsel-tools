import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

from shnitsel.dynamic import postprocess as P


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
    return P.broaden_gauss(data.energy, data.fosc, width=0.1, agg_dim='trajid')


def inlabel(s, ax, ha='center', va='center'):
    return ax.text(
        0.05,
        0.95,
        s,
        fontweight='bold',
        transform=ax.transAxes,
        ha=ha,
        va=va,
    )


def ski_plots(spectra: xr.DataArray) -> mpl.figure.Figure:
    assert 'time' in spectra.coords, "Missing 'time' coordinate"
    assert 'statecomb' in spectra.coords, "Missing 'statecomb' coordinate"
    assert 'energy' in spectra.coords, "Missing 'energy' coordinate"

    nstatecombs = spectra.sizes['statecomb']
    fig, axs = plt.subplots(nstatecombs, 1, layout='constrained', sharex=True)
    fig.set_size_inches(6, 10)

    cnorm = mpl.colors.Normalize(spectra.time.min(), spectra.time.max())
    cmap = plt.get_cmap('viridis')
    for ax, (sc, scdata) in zip(axs, spectra.groupby('statecomb')):
        for t, tdata in scdata.groupby('time'):
            ax.plot(tdata.energy, tdata.squeeze(), c=cmap(cnorm(t)), linewidth=0.2)
        maxes = scdata[scdata.argmax('energy')]
        ax.plot(
            maxes.energy.squeeze(),
            maxes.squeeze(),
            c='black',
            linewidth=1,
            linestyle='--',
        )

        inlabel(sc, ax)
        ax.set_ylabel(r'$f_\mathrm{osc}$')
    ax.set_xlabel(r'$E$ / eV')
    return fig


def pcm_plots(spectra: xr.DataArray) -> mpl.figure.Figure:
    assert 'time' in spectra.coords, "Missing 'time' coordinate"
    assert 'statecomb' in spectra.coords, "Missing 'statecomb' coordinate"
    assert 'energy' in spectra.coords, "Missing 'energy' coordinate"

    nstatecombs = spectra.sizes['statecomb']
    fig, axs = plt.subplots(1, nstatecombs, layout='constrained')

    cnorm = mpl.colors.LogNorm(0.0005, spectra.max())
    for ax, (sc, scdata) in zip(axs, spectra.groupby('statecomb')):
        qm = scdata.squeeze().plot.pcolormesh(x='energy', y='time', ax=ax, norm=cnorm)
        qm.axes.invert_yaxis()
    return fig