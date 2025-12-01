from .common import figaxs_defaults, centertext
from .hist import create_marginals, trunc_max

from matplotlib.axes import Axes


@figaxs_defaults(mosaic=[['ntd'], ['nde']], scale_factors=(1 / 3, 1 / 3))
def plot_nacs_histograms(inter_state, hop_idxs, fig=None, axs=None) -> dict[str, Axes]:
    """Plot 2D histograms of NACS vs delta_E or dip_trans"""

    nacs_data = inter_state.sel(frame=hop_idxs)
    axs['nde'].set_ylabel(r'$\Delta E$ / eV')
    axs['nde'].minorticks_on()
    axs['nde'].set_xlabel(r"$\|\mathrm{NAC}_{i,j}\|_2$")
    axs['ntd'].tick_params(axis="x", labelbottom=False)
    if 'dip_trans' in inter_state:
        axs['ntd'].set_ylabel(r"$\|\mathbf{\mu}_{i,j}\|_2$")
        axs['ntd'].minorticks_on()

    def plot(label, yname):
        ax = axs[label]
        axx, axy = create_marginals(ax)
        bins = 100

        for i, (sc, data) in enumerate(nacs_data.groupby('statecomb')):
            ydata = data[yname].squeeze()
            xdata = data['nacs'].squeeze()
            # xmax = trunc_max(xdata)
            ymax = trunc_max(ydata)
            color = data['_color'].item()
            axx.hist(
                xdata,
                # range=(0, xmax),
                color=color,
                bins=bins,
            )
            axy.hist(
                ydata, range=(0, ymax), orientation='horizontal', color=color, bins=bins
            )

            ax.scatter(xdata, ydata, color=color, s=0.2, alpha=0.5)

    plot('nde', 'energy')
    if 'dip_trans' in inter_state:
        plot('ntd', 'dip_trans')
    else:
        centertext(r"No $\mathbf{\mu}_{ij}$ data", axs['ntd'])
        # axs['ntd'].tick_params(axis='y', labelleft=False)
        axs['ntd'].get_yaxis().set_visible(False)
        axs['ntd'].get_xaxis().set_visible(False)

    return axs