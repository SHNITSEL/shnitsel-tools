import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


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

def ski_plots(spectra: xr.DataArray, threshold: float = np.inf) -> mpl.figure.Figure:
    """Plot spectra for different times on top of each other,
    along with a dashed "ski"-line that tracks the maximum.
    One plot per statecomb; plots stacked vertically.
    Expected to be used on data produced by ``spectra.spectra_all_times``.

    Parameters
    ----------
    spectra
        DataArray containing fosc values organized along 'energy', 'time' and
        'statecomb' dimensions.
    threshold
        The "ski" line will not be drawn between successive
        points for which the change in energy-coordinate is greater than this.

    Returns
    -------
        Figure object corresponding to plot.

    .. !the example is commented out for now

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
    assert 'energy_interstate' in spectra.coords, (
        "Missing 'energy_interstate' coordinate"
    )

    nstatecombs = spectra.sizes['statecomb']
    fig, axs = plt.subplots(nstatecombs, 1, layout='constrained', sharex=True)
    fig.set_size_inches(6, 10)

    cnorm = mpl.colors.Normalize(spectra.time.min(), spectra.time.max())
    cmap = plt.get_cmap('viridis')

    if nstatecombs == 1:
        axs = [axs]

    for ax, (sc, scdata) in zip(axs, spectra.groupby('statecomb')):
        scdata = scdata.squeeze('statecomb')
        for t, tdata in scdata.groupby('time'):
            ax.plot(
                tdata['energy_interstate'],
                tdata.squeeze('time'),
                c=cmap(cnorm(t)),
                linewidth=0.2,
            )
        maxes = scdata[scdata.argmax('energy_interstate')]
        xs = maxes['energy_interstate'].data
        ys = maxes.data

        segments = np.c_[xs[:-1], ys[:-1], xs[1:], ys[1:]].reshape(-1, 2, 2)
        mask = np.abs(segments[:, 0, 0] - segments[:, 1, 0]) < threshold
        segments = segments[mask]
        lc = mpl.collections.LineCollection(
            segments,
            color='k',
            linewidths=1,
            linestyles='--',
            # cmap=cmap, norm=cnorm,
        )
        lc.set_array(maxes['time'].data)
        ax.add_collection(lc)

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
    assert 'energy_interstate' in spectra.coords, (
        "Missing 'energy_interstate' coordinate"
    )

    nstatecombs = spectra.sizes['statecomb']
    fig, axs = plt.subplots(1, nstatecombs, layout='constrained', sharey=True)

    cnorm = mpl.colors.LogNorm(5e-4, spectra.max())
    
    if nstatecombs == 1:
        axs = [axs]
    for ax, (sc, scdata) in zip(axs, spectra.groupby('statecomb')):
        qm = scdata.squeeze('statecomb').plot.pcolormesh(
            x='energy_interstate', y='time', ax=ax, norm=cnorm
        )
        qm.axes.invert_yaxis()
        ax.set_title(str(sc))  # TODO (thevro): Use TeX state names
    return fig