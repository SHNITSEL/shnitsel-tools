import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from .common import figaxs_defaults
from .hist import trunc_max, create_marginals
from .colormaps import magma_rw, custom_ylgnr


def single_hist(data, shi, slo, color, bins=100, ax=None, cmap=None, cnorm=None):
    if ax is None:
        _, ax = plt.subplots(1, 1)
    if cmap is None:
        cmap = magma_rw

    axx, axy = create_marginals(ax)
    xdata = data['energy'].squeeze()
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


def plot_dip_trans_histograms(inter_state, axs=None, cnorm=None):
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

        color = data['_color'].item()
        hist2d_outputs.append(
            single_hist(data, shi, slo, color=color, ax=ax, cnorm=cnorm)
        )
    return hist2d_outputs


def plot_spectra(spectra, ax=None, cmap=None, cnorm=None, mark_peaks=False):
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
        # special casing for now
        if sc == '$S_2 - S_0$':
            linestyle = '--'
        else:
            linestyle = '-'
        c = cmap(cnorm(t))
        # ax.fill_between(data['energy'], data, alpha=0.5, color=c)
        ax.plot(
            data['energy'],
            data,
            # linestyle=linestyles[t], c=dcol_inter[sc],
            linestyle=linestyle,
            c=c,
            linewidth=0.5,
        )
        if mark_peaks:
            try:
                peak = data[data.argmax('energy')]
                ax.text(peak['energy'], peak, f"{t:.2f}:{sc}", fontsize='xx-small')
            except Exception as e:
                print(e)
    _, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax)

    return ax


@figaxs_defaults(
    mosaic=[['sg'], ['t0'], ['t1'], ['se'], ['t2'], ['cb_spec'], ['cb_hist']],
    scale_factors=(1 / 3, 4 / 5),
    height_ratios=([1] * 5) + ([0.1] * 2),
)
def plot_separated_spectra_and_hists(
    inter_state, sgroups, fig=None, axs=None, cb_spec_vlines=True
):
    ground, excited = sgroups
    times = [tup[0] for lst in sgroups for tup in lst]
    scnorm = plt.Normalize(inter_state.time.min(), inter_state.time.max())
    scmap = plt.get_cmap('turbo')
    scscale = mpl.cm.ScalarMappable(norm=scnorm, cmap=scmap)

    hist2d_outputs = []
    # ground-state spectra and histograms
    plot_spectra(ground, ax=axs['sg'], cnorm=scnorm, cmap=scmap)

    # We show at most the first two statecombs
    if inter_state.sizes['statecomb'] >= 2:
        selsc = [0, 1]
        selaxs = [axs['t1'], axs['t0']]
    elif inter_state.sizes['statecomb'] == 1:
        selsc = [0]
        selaxs = [axs['t1']]
    else:
        raise ValueError(
            "Too few statecombs (expecting at least 2 states => 1 statecomb)"
        )
    hist2d_outputs += plot_dip_trans_histograms(
        inter_state.isel(statecomb=selsc),
        axs=selaxs,
    )

    # excited-state spectra and histograms
    if inter_state.sizes['statecomb'] >= 2:
        plot_spectra(excited, ax=axs['se'], cnorm=scnorm, cmap=scmap)
        hist2d_outputs += plot_dip_trans_histograms(
            inter_state.isel(statecomb=[2]), axs=[axs['t2']]
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

    cb_spec = axs['cb_spec'].figure.colorbar(
        scscale,
        cax=axs['cb_spec'],
        location='bottom',
        extend='both',
        extendrect=True,
    )
    axs['cb_spec'].set_xlabel('time / fs')
    if cb_spec_vlines:
        for t in times:
            lo, hi = scscale.get_clim()  # lines at these points don't show
            if t == lo:
                t += t / 100  # so we shift them slightly
            elif t == hi:
                t -= t / 100

            cb_spec.ax.axvline(t, c='white', linewidth=0.5)

    hcscale = mpl.cm.ScalarMappable(norm=hcnorm, cmap=magma_rw)
    axs['cb_hist'].figure.colorbar(hcscale, cax=axs['cb_hist'], location='bottom')
    axs['cb_hist'].set_xlabel('# data points')

    axs['se'].set_title(
        r"$\uparrow$ground state" + "\n" + r"$\downarrow$excited state absorption"
    )
    axs['t2'].set_xlabel(r'$\Delta E$ / eV')

    legend_lines, legend_labels = zip(
        *[
            (Line2D([0], [0], color='k', linestyle='-', linewidth=0.5), "$S_1/S_0$"),
            (Line2D([0], [0], color='k', linestyle='--', linewidth=0.5), "$S_2/S_0$"),
        ]
    )
    axs['sg'].legend(legend_lines, legend_labels, fontsize='x-small')

    return axs