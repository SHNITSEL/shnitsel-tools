import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

from shnitsel.core import postprocess as P


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
    return P.broaden_gauss(data.energy, data.fosc, agg_dim='trajid')


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
    """Plot spectra for different times on top of each other,
    along with a dashed line that tracks the maximum.
    One plot per statecomb; plots stacked vertically.
    Expected to be used on data produced by `spectra3d.spectra_all_times`.

    Parameters
    ----------
    spectra
        DataArray containing fosc values organized along 'energy', 'time' and
        'statecomb' dimensions.

    Returns
    -------
        Figure object corresponding to plot.

    Examples
    --------
        >>> from shnitsel.dynamic import xrhelpers as xh, postprocess as P
        >>> from shnitsel.dynamic.plot import spectra3d
        >>> spectra_data = (
                xh.open_frames(path)
                .pipe(P.get_inter_state)
                .pipe(P.assign_fosc)
                .pipe(spectra3d.spectra_all_times))
        >>> spectra3d.ski_plots(spectra_data)
    """
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
    """Represent fosc as colour in a plot of fosc against time and energy.
    The colour scale is logarithmic.
    One plot per statecomb; plots stacked horizontally.
    Expected to be used on data produced by `spectra3d.spectra_all_times`.

    Parameters
    ----------
    spectra
        DataArray containing fosc values organized along 'energy', 'time' and
        'statecomb' dimensions.

    Returns
    -------
        Figure object corresponding to plot.

    Examples
    --------
        >>> from shnitsel.dynamic import xrhelpers as xh, postprocess as P
        >>> from shnitsel.dynamic.plot import spectra3d
        >>> spectra_data = (
                xh.open_frames(path)
                .pipe(P.get_inter_state)
                .pipe(P.assign_fosc)
                .pipe(spectra3d.spectra_all_times))
        >>> spectra3d.pcm_plots(spectra_data)
    """
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