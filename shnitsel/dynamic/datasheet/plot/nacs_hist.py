import matplotlib.pyplot as plt

from .hist import create_marginals, trunc_max


def plot_nacs_histograms(inter_state, hop_idxs, col_inter, axs=None):
    """Plot 2D histograms of NACS vs delta_E or dip_trans"""

    if axs is None:
        fig, axs = plt.subplot_mosaic(
            [
                ['energies', 'forces', 'dip_perm'],
                ['fde', 'pca', 'pca'],
                ['t0', 'pca', 'pca'],
                ['t1', '.', 'pop'],
                ['t2', 'ntd', 'de'],
                ['.', 'nde', 'ft'],
            ],
            layout='constrained',
        )
        fig.set_size_inches(8.27, 11.69)  # portrait A4
    else:
        fig = axs['nde'].figure

    nacs_data = inter_state.sel(frame=hop_idxs)
    axs['nde'].set_ylabel(r'$\Delta E$ / eV')
    axs['nde'].minorticks_on()
    axs['nde'].tick_params(axis="x", labelbottom=False)
    if 'dip_trans' in inter_state:
        axs['ntd'].set_ylabel(r"$\|\mathbf{\mu}_{i,j}\|_2$")
        axs['ntd'].set_xlabel(r"$\|\mathrm{NAC}_{i,j}\|_2$")
        axs['ntd'].minorticks_on()

    def plot(label, yname):
        ax = axs[label]
        axx, axy = create_marginals(ax)
        bins = 100
        # for sc, data in inter_state.groupby('statecomb'):
        for i, (sc, data) in enumerate(nacs_data.groupby('statecomb')):
            if i != 0:
                continue
            ydata = data[yname].squeeze()
            xdata = data['nacs'].squeeze()
            xmax = trunc_max(xdata)
            ymax = trunc_max(ydata)
            # ymax = trunc_max(ydata) # if you truncate nacs, there's nothing left
            color = col_inter[i]
            axx.hist(xdata, range=(0, xmax), color=color, bins=bins)
            axy.hist(
                ydata, range=(0, ymax), orientation='horizontal', color=color, bins=bins
            )

            ax.scatter(xdata, ydata, color=color, s=0.2, alpha=0.5)
            # ax.set_xlim(0, xmax)

    plot('nde', 'energies')
    if 'dip_trans' in inter_state:
        plot('ntd', 'dip_trans')

    return axs