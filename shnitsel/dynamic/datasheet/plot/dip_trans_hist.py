import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from .hist import trunc_max, create_marginals
from .colormaps import magma_rw, custom_ylgnr


def single_hist(data, shi, slo, color, bins=100, ax=None, cmap=None, cnorm=None):
    if ax is None:
        _, ax = plt.subplots(1, 1)
    if cmap is None:
        cmap = magma_rw

    axx, axy = create_marginals(ax)
    xdata = data['energies'].squeeze()
    ydata = data['dip_trans'].squeeze()
    xmax = trunc_max(xdata)
    ymax = trunc_max(ydata)
    axx.hist(xdata, range=(0, xmax), color=color, bins=bins)
    axy.hist(ydata, range=(0, ymax), orientation='horizontal', color=color, bins=bins)
    hist2d_output = ax.hist2d(
        xdata, ydata, range=[(0, xmax), (0, ymax)], bins=bins, cmap=cmap, norm=cnorm
    )

    ax.set_ylabel(r"$\|\mathbf{\mu}_{%d,%d}\|_2$" % (shi, slo))
    ax.text(
        1.05,
        1.05,
        "$S_%d/S_%d$" % (shi, slo),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        color=color,
        #   fontweight='bold',
    )

    return hist2d_output


def plot_dip_trans_histograms(inter_state, dcol_inter, axs=None, cnorm=None):
    if axs is None:
        nplots = len(inter_state.coords['statecomb'])
        _, axs = plt.subplots(nplots, 1, layout='constrained')

    # TODO obviate following cludge:
    sclabels = [(int(x[3]), int(x[9])) for x in inter_state.statecomb.values]

    hist2d_outputs = []
    for i, (sc, data) in enumerate(inter_state.groupby('statecomb')):
        # label = f't{i}'
        shi, slo = sclabels[i]
        ax = axs[i]

        color = dcol_inter[sc]
        hist2d_outputs.append(
            single_hist(data, shi, slo, color=color, ax=ax, cnorm=cnorm)
        )
    return hist2d_outputs


def plot_spectra(spectra, dcol_inter, ax=None, cmap=None, cnorm=None):
    if ax is None:
        _, ax = plt.subplots(1, 1)
    cmap = plt.get_cmap(cmap) if cmap else custom_ylgnr
    times = [t for (t, sc) in spectra]
    cnorm = cnorm if cnorm else plt.Normalize(min(times), max(times))
    ax.set_ylabel(r'$f_\mathrm{osc}$')
    ax.invert_xaxis()
    # linestyles = {t: ['-', '--', '-.', ':'][i]
    #               for i, t in enumerate(np.unique(list(zip(*spectra.keys()))[0]))}
    for (t, sc), data in spectra.items():
        c = cmap(cnorm(t))
        # ax.fill_between(data['energy'], data, alpha=0.5, color=c)
        ax.plot(
            data['energy'],
            data,
            # linestyle=linestyles[t], c=dcol_inter[sc],
            c=c,
            linewidth=0.5,
        )

    return ax


def _hist_min_max(inter_state):
    minmax = []
    debug = []
    for sc, data in inter_state.groupby('statecomb'):
        data = data.squeeze()
        H, xe, ye = np.histogram2d(
            data.energies.values, data.dip_trans.values, bins=100
        )
        minmax += [H.min(), H.max()]
        debug.append((H, xe, ye))

    return min(minmax), max(minmax), debug


def plot_separated_spectra_and_hists(inter_state, sgroups, dcol_inter, axs=None):
    if axs is None:
        mosaic = [['sg'], ['t0'], ['t1'], ['se'], ['t2'], ['cb_spec'], ['cb_hist']]
        fig, axs = plt.subplot_mosaic(
            mosaic, layout='constrained', height_ratios=([1] * 5) + ([0.1] * 2)
        )
        fig.set_size_inches(8.27 / 3, 11.69 * 0.8)  # portrait A4

    ground, excited = sgroups
    times = [tup[0] for lst in sgroups for tup in lst]
    scnorm = plt.Normalize(min(times), max(times))
    scmap = plt.get_cmap('viridis_r')
    scscale = mpl.cm.ScalarMappable(norm=scnorm, cmap=scmap)

    hist2d_outputs = []
    # ground-state spectra and histograms
    plot_spectra(ground, dcol_inter, ax=axs['sg'], cnorm=scnorm, cmap=scmap)
    hist2d_outputs += plot_dip_trans_histograms(
        inter_state.isel(statecomb=[0, 1]),
        dcol_inter,
        axs=[axs[k] for k in ['t1', 't0']],
    )
    # excited-state spectra and histograms
    plot_spectra(excited, dcol_inter, ax=axs['se'], cnorm=scnorm, cmap=scmap)
    hist2d_outputs += plot_dip_trans_histograms(
        inter_state.isel(statecomb=[2]), dcol_inter, axs=[axs['t2']]
    )

    hists = np.array([tup[0] for tup in hist2d_outputs])
    hcnorm = plt.Normalize(hists.min(), hists.max())

    quadmeshes = [tup[3] for tup in hist2d_outputs]
    for quadmesh in quadmeshes:
        quadmesh.set_norm(hcnorm)

    def ev2nm(ev):
        return 4.135667696 * 2.99792458 * 100 / np.where(ev != 0, ev, 1)

    lims = [l for ax in axs.values() for l in ax.get_xlim()]
    new_lims = (min(lims), max(lims))
    for lax, ax in axs.items():
        if lax.startswith('cb'):
            continue
        ax.set_xlim(*new_lims)
        ax.invert_xaxis()

    for ax in list(axs.values()):
        ax.tick_params(axis="x", labelbottom=False)
    axs['t2'].tick_params(axis="x", labelbottom=True)

    secax = axs['sg'].secondary_xaxis('top', functions=(ev2nm, ev2nm))
    secax.set_xticks([50, 75, 100, 125, 150, 200, 300, 500, 1000])
    secax.tick_params(axis='x', rotation=45, labelsize='small')
    for l in secax.get_xticklabels():
        l.set_horizontalalignment('left')
        l.set_verticalalignment('bottom')
    secax.set_xlabel(r'$\Delta E$ / nm')

    for lax in ['cb_spec', 'cb_hist']:
        axs[lax].get_yaxis().set_visible(False)

    axs['cb_spec'].figure.colorbar(scscale, cax=axs['cb_spec'], location='bottom')
    axs['cb_spec'].set_xlabel('time / fs')
    hcscale = mpl.cm.ScalarMappable(norm=hcnorm, cmap=magma_rw)
    axs['cb_hist'].figure.colorbar(hcscale, cax=axs['cb_hist'], location='bottom')
    axs['cb_hist'].set_xlabel('# data points')

    axs['se'].set_title(
        r"$\uparrow$ground state" + "\n" + r"$\downarrow$excited state absorption"
    )
    axs['t2'].set_xlabel(r'$\Delta E$ / eV')

    return axs