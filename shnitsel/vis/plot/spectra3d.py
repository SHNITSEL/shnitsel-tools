import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

# TODO: Use plot.common.inlabel instead?
def _inlabel(s, ax, ha='center', va='center'):
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
    Expected to be used on data produced by ``spectra.spectra_all_times``.

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
        >>> import shnitsel as st
        >>> from shnitsel.core.plot import spectra3d
        >>> spectra_data = (
                st.io.read(path)
                .st.get_inter_state()
                .st.assign_fosc()
                .st.spectra_all_times())
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

    if nstatecombs == 1:
        axs = [axs]

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

        _inlabel(sc, ax)
        ax.set_ylabel(r'$f_\mathrm{osc}$')
    ax.set_xlabel(r'$E$ / eV')
    return fig


def pcm_plots(spectra: xr.DataArray) -> mpl.figure.Figure:
    """Represent fosc as colour in a plot of fosc against time and energy.
    The colour scale is logarithmic.
    One plot per statecomb; plots stacked horizontally.
    Expected to be used on data produced by `spectra.spectra_all_times`.

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
        >>> import shnitsel as st
        >>> from shnitsel.core.plot import spectra3d
        >>> spectra_data = (
                st.io.read(path)
                .st.get_inter_state()
                .st.assign_fosc()
                .st.spectra_all_times())
        >>> spectra3d.pcm_plots(spectra_data)
    """
    assert 'time' in spectra.coords, "Missing 'time' coordinate"
    assert 'statecomb' in spectra.coords, "Missing 'statecomb' coordinate"
    assert 'energy' in spectra.coords, "Missing 'energy' coordinate"

    nstatecombs = spectra.sizes['statecomb']
    fig, axs = plt.subplots(1, nstatecombs, layout='constrained')

    cnorm = mpl.colors.LogNorm(0.0005, spectra.max())
    
    if nstatecombs == 1:
        axs = [axs]
    for ax, (sc, scdata) in zip(axs, spectra.groupby('statecomb')):
        qm = scdata.squeeze().plot.pcolormesh(x='energy', y='time', ax=ax, norm=cnorm)
        qm.axes.invert_yaxis()
    return fig