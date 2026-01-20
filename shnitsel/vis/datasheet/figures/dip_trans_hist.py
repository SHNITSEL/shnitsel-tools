import logging
import matplotlib as mpl
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, Normalize
from matplotlib.figure import Figure, SubFigure
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from shnitsel.data.dataset_containers.inter_state import InterState
from shnitsel.filtering.state_selection import StateSelection

from ....core.typedefs import SpectraDictType
from ....units.definitions import energy

from .common import figaxs_defaults
from .hist import calc_truncation_maximum, create_marginals
from ...colormaps import magma_rw, custom_ylgnr
from ....units.conversion import convert_energy


def single_dip_trans_hist(
    interstate: InterState,
    sc_label: str,
    state_labels: tuple[str, str],
    color: str,
    bins: int = 100,
    ax: Axes | None = None,
    cmap=None,
    cnorm=None,
):
    """Function to plot a single histogram of interstate dip_trans data.

    Parameters
    ----------
    interstate : InterState
        Inter-state Dataset
    sc_label : str
        Label to use for the state combination.
    state_labels : tuple[str,str]
        Labels for the individual states.
    color : str): Color for the histogram of this state combination.
    bins : int, optional
        Number of bins for the histogram. Defaults to 100.
    ax : Axes, optional
        Axes object to plot into. Defaults to None.
    cmap : str, optional
        Colormap to use. Defaults to None.
    cnorm : str, optional
        Norming method to apply to the colormap. Defaults to None.

    Returns
    -------
        ?: The result of ax.hist2d will be returned
    """
    if ax is None:
        _, ax = plt.subplots(1, 1)
    if cmap is None:
        cmap = magma_rw

    axx, axy = create_marginals(ax)
    # We expect energies in eV for the plot
    xdata = interstate['energy_interstate'].squeeze()
    xdata = convert_energy(xdata, to=energy.eV)

    # We need the normed transition dipole
    ydata = interstate['dip_trans_norm'].squeeze()

    xmax = calc_truncation_maximum(xdata)
    ymax = calc_truncation_maximum(ydata)
    axx.hist(xdata, range=(0, xmax), color=color, bins=bins)
    axy.hist(ydata, range=(0, ymax), orientation='horizontal', color=color, bins=bins)
    hist2d_output = ax.hist2d(
        xdata, ydata, range=[(0, xmax), (0, ymax)], bins=bins, cmap=cmap, norm=cnorm
    )

    ax.set_ylabel(
        r"$\|\mathbf{\mu}_{%s,%s}\|_2$" % (state_labels[0], state_labels[1])
    )  # r"$\|\mathbf{\mu}_{%d,%d}\|_2$" % (sc[0], sc[1]))
    ax.text(
        1.05,
        1.05,
        "$%s/%s$"
        % (
            state_labels[0],
            state_labels[1],
        ),  # f"{sc_label}",  # $S_%d/S_%d$" % (sc[0], sc[1]),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        color=color,
        #   fontweight='bold',
    )

    return hist2d_output


def plot_dip_trans_histograms(
    inter_state: InterState,
    state_selection: StateSelection,
    axs: list[Axes] | None = None,
    cnorm: str | None = None,
) -> list:
    """function to plot all relevant histograms for the provided inter_state data

    Parameters
    ----------
    inter_state : InterState
        Inter-state data to get the transitional dipole data from.
    state_selection : StateSelection
        State selection object to limit the states included in plotting and to provide state names.
    axs : Axes, optional
        Axes objects to plot into. If not provided, will be created.
    cnorm : str, optional
        Optional specification of a colormap norm method. Defaults to None.

    Returns
    -------
    list
        The list of the results of hist2d() calls for the provided data in inter_state
    """
    if axs is None:
        state_combs = list(state_selection.combination_info())
        nplots = len(state_combs)
        _, axs = plt.subplots(nplots, 1, layout='constrained')

    assert axs is not None, "Could not create subplot axes."

    hist2d_outputs = []
    selected_scs = 0
    for i, (sc_, data) in enumerate(inter_state.groupby('statecomb')):
        if not state_selection.has_state_combination(sc_):
            continue

        if selected_scs >= len(axs):
            break

        # label = f't{i}'
        sc_label = state_selection.get_state_combination_tex_label(sc_)  # = sclabels[i]
        state_labels = (
            state_selection.get_state_tex_label(sc_[0]),
            state_selection.get_state_tex_label(sc_[1]),
        )
        ax = axs[selected_scs]

        color = state_selection.get_state_combination_color(sc_)
        hist2d_outputs.append(
            single_dip_trans_hist(
                data, sc_label, state_labels, color=color, ax=ax, cnorm=cnorm
            )
        )
        selected_scs += 1
    return hist2d_outputs


def plot_spectra(
    spectra: SpectraDictType,
    state_selection: StateSelection,
    ax: Axes | None = None,
    lim_num_sc: int = -1,
    cmap: str | Colormap | None = None,
    cnorm: str | Normalize | None = None,
    mark_peaks: bool = False,
) -> Axes:
    """Create the spectra plot of the system denoted by the results in spectra.

    Parameters
    ----------
    spectra : SpectraDictType
        The spectra (t, state combination) -> fosc data to plot. t is expected in fs.
    ax : Axes, optional
        Axis object to plot into. If not provided, will be created.
    lim_num_sc : int, optional
        Maximum number of state combinations to consider for the plot.
    state_selection : StateSelection
        State selection object to limit the states included in plotting and to provide state names.
    cmap : str | Colormap, optional
        Optional specification of a desired colormap. Defaults to None.
    cnorm : str | Normalize, optional
        Optional specification of a colormap norm method. Defaults to None.
    mark_peaks : bool, optional
        Flag whether peaks should be clearly marked. Defaults to False.

    Returns
    -------
    Axes
        The axes object into which the graph was plotted
    """
    if ax is None:
        _, ax = plt.subplots(1, 1)
    # cmap = plt.get_cmap(cmap) if cmap else custom_ylgnr
    # times = [t for (t, sc) in spectra]
    # cnorm = cnorm if cnorm else plt.Normalize(min(times), max(times))
    ax.set_ylabel(r'$f_\mathrm{osc}$')
    ax.invert_xaxis()
    # linestyles = {t: ['-', '--', '-.', ':'][i]
    #               for i, t in enumerate(np.unique(list(zip(*spectra.keys()))[0]))}

    times = list(set([tup[0] for tup in spectra]))
    times.sort()

    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']

    times_styles = {t: linestyles[i % len(linestyles)] for i, t in enumerate(times)}
    sc_count = {}
    for i, ((t, sc), data) in enumerate(spectra.items()):
        if not state_selection.has_state_combination(sc):
            continue

        if data.isnull().all():
            continue

        if sc not in sc_count:
            if lim_num_sc > 0 and len(sc_count) >= lim_num_sc:
                continue
            sc_count[sc] = 0

        # curr_count = sc_count[sc]
        sc_count[sc] += 1

        sc_color = state_selection.get_state_combination_color(sc)
        # sc_label = state_selection.get_state_combination_tex_label(sc)
        # special casing for now
        linestyle = times_styles[t]

        # c = cmap(cnorm(t))
        # ax.fill_between(data['energy'], data, alpha=0.5, color=c)
        converted_energy = convert_energy(data['energy_interstate'], to=energy.eV)
        converted_energy = np.abs(converted_energy)
        ax.plot(
            converted_energy,
            data,
            # linestyle=linestyles[t], c=dcol_inter[sc],
            linestyle=linestyle,
            c=sc_color,
            linewidth=0.8,
        )
        if mark_peaks:
            try:
                peak = data[data.argmax('energy_interstate')]
                ax.text(
                    float(
                        convert_energy(peak['energy_interstate'], to=energy.eV).values
                    ),
                    float(peak),
                    f"{t:.2f}:{sc}",
                    fontsize='xx-small',
                )
            except Exception as e:
                logging.warning(f"{e}")

    handles = []
    labels = []
    for t, style in times_styles.items():
        handles.append(Line2D([0], [0], color='k', ls=style))
        labels.append(f'$t={t:.1f}\\,fs$')
    _, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax)
    # ax.legend(
    #     handles=handles,
    #     labels=labels,
    #     loc='center left',
    #     bbox_to_anchor=(1.1, 0.5),
    #     prop={'size': 6},
    # )

    return ax


@figaxs_defaults(
    mosaic=[['sg'], ['t0'], ['t1'], ['se'], ['t2'], ['cb_spec'], ['cb_hist']],
    scale_factors=(1 / 3, 4 / 5),
    height_ratios=([1] * 5) + ([0.1] * 2),
)
def plot_separated_spectra_and_hists(
    inter_state: InterState,
    spectra_groups: tuple[SpectraDictType, SpectraDictType],
    state_selection: StateSelection,
    fig: Figure | SubFigure | None = None,
    axs: dict[str, Axes] | None = None,
    cb_spec_vlines: bool = True,
):
    """Create separate spectra plots for ground and excited states.

    Parameters
    ----------
    inter_state : InterState
        The interstate data to use for the spectra plots
    spectra_groups : tuple[SpectraDictType, SpectraDictType]
        Spectra separated into `ground` state spectra and `excited` spectra.
    state_selection : StateSelection
        State selection object to limit the states included in plotting and to provide state names.
    fig : Figure| SubFigure
        A figure, consumed by the automatic axes generation. Not used by the function itself.
    axs : dict[str,Axes], optional
        Axis dictionary object to plot into.
    cb_spec_vlines : bool, optional
        Whether to mark spectral lines in the energy spectrum. Defaults to True.

    Returns
    -------
    dict[str, Axes]
        The axes dict loc='upper right', bbox_to_anchor=(1.1, 1.1))after plotting to it.
    """
    assert axs is not None, "Could not acquire axes for plotting"
    ground, excited = spectra_groups
    times = [tup[0] for lst in spectra_groups for tup in lst]
    scnorm = plt.Normalize(inter_state.time.min(), inter_state.time.max())
    scmap = plt.get_cmap('turbo')
    scscale = mpl.cm.ScalarMappable(norm=scnorm, cmap=scmap)

    hist2d_outputs = []
    # ground-state spectra and histograms
    plot_spectra(
        ground, ax=axs['sg'], state_selection=state_selection, cnorm=scnorm, cmap=scmap
    )

    # # We show at most the first two statecombs
    # if inter_state.sizes['statecomb'] >= 2:
    #     selsc = [0, 1]
    #     selaxs = [axs['t1'], axs['t0']]
    # elif inter_state.sizes['statecomb'] == 1:
    #     selsc = [0]
    #     selaxs = [axs['t1']]
    # else:
    #     raise ValueError(
    #         "Too few statecombs (expecting at least 2 states => 1 statecomb)"
    #     )
    # inter_state_sel = inter_state.isel(statecomb=selsc)

    if len(state_selection.state_combinations) > 1:
        selaxs = [axs['t1'], axs['t0']]
    else:
        selaxs = [axs['t1']]

    hist2d_outputs += plot_dip_trans_histograms(
        inter_state,  # inter_state_sel,
        axs=selaxs,
        state_selection=state_selection,
    )

    # excited-state spectra and histograms
    if len(excited) >= 2:
        plot_spectra(
            excited,
            ax=axs['se'],
            state_selection=state_selection,
            cnorm=scnorm,
            cmap=scmap,
        )
        hist2d_outputs += plot_dip_trans_histograms(
            inter_state,  # inter_state.isel(statecomb=[2]),
            axs=[axs['t2']],
            state_selection=state_selection,
        )

    hists = np.array([tup[0] for tup in hist2d_outputs])

    if len(hists) > 0:
        hcnorm = plt.Normalize(hists.min(), hists.max())

        quadmeshes = [tup[3] for tup in hist2d_outputs]
        for quadmesh in quadmeshes:
            quadmesh.set_norm(hcnorm)
        hcscale = mpl.cm.ScalarMappable(norm=hcnorm, cmap=magma_rw)
        axs['cb_hist'].figure.colorbar(hcscale, cax=axs['cb_hist'], location='bottom')

    def ev2nm(delta_E):
        """Helper function to convert delta energy of a transition in eV to a nanometer wavelength

        Args:
            delta_E (float-like): Energy delta in electron Volts (eV) units.

        Returns:
            float-like: Resulting wavelength in nm.
        """
        return 4.135667696 * 2.99792458 * 100 / np.where(delta_E != 0, delta_E, 1)

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
    secax.set_xlabel(r'$\lambda$ / nm')

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


@figaxs_defaults(
    # mosaic=[['sg'], ['t0'], ['t1'], ['se'], ['t2'], ['cb_spec'], ['cb_hist']],
    mosaic=[['sg', 't1'], ['cb_spec', 'cb_hist']],
    scale_factors=(4 / 5, 4 / 5),
    # height_ratios=([1] * 5) + ([0.1] * 2),
    height_ratios=([1]) + ([0.1]),
)
def plot_separated_spectra_and_hists_groundstate(
    inter_state: InterState,
    spectra_groups: tuple[SpectraDictType, SpectraDictType],
    state_selection: StateSelection,
    fig: Figure | SubFigure | None = None,
    axs: dict[str, Axes] | None = None,
    cb_spec_vlines: bool = True,
    scmap: Colormap = plt.get_cmap('turbo'),
) -> dict[str, Axes]:
    """Function to plot separated spectra and histograms of ground state data only.


    Parameters
    ----------
    inter_state : InterState
        Inter-State dataset containing energy differences
    spectra_groups : tuple[SpectraDictType, SpectraDictType]
        Tuple holding the spectra groups of ground-state transitions and excited-state transitions.
    state_selection : StateSelection
        State selection object to limit the states included in plotting and to provide state names.
    fig : Figure | SubFigure | None, optional
        Figure to plot the graphs to. Defaults to None.
    axs : dict[str, Axes] | None, optional
        Dict of named axes to plot to. Defaults to None.
    cb_spec_vlines : bool, optional
        Flag to enable vertical lines in the time-dependent spectra. Defaults to True.
    scmap : Colormap, optional
        State combination colormap. Defaults to plt.get_cmap('turbo').

    Raises
    ------
    AssertionError
        Could not acquire axes to plot to

    Returns
    -------
        dict[str, Axes]: The named axes after plotting
    """
    assert axs is not None, "Could not acquire axes to plot to."
    ground, excited = spectra_groups
    times = [tup[0] for lst in spectra_groups for tup in lst]
    scnorm = plt.Normalize(inter_state.time.min(), inter_state.time.max() + 50)
    scscale = mpl.cm.ScalarMappable(norm=scnorm, cmap=scmap)

    hist2d_outputs = []
    # ground-state spectra and histograms
    plot_spectra(
        ground, ax=axs['sg'], state_selection=state_selection, cnorm=scnorm, cmap=scmap
    )

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
        inter_state.isel(statecomb=selsc), axs=selaxs, state_selection=state_selection
    )

    # excited-state spectra and histograms
    # if inter_state.sizes['statecomb'] >= 2:
    #    plot_spectra(excited, ax=axs['se'], cnorm=scnorm, cmap=scmap)
    #    hist2d_outputs += plot_dip_trans_histograms(
    #        inter_state.isel(statecomb=[2]), axs=[axs['t2']]
    #    )

    hists = np.array([tup[0] for tup in hist2d_outputs])
    hcnorm = plt.Normalize(hists.min(), hists.max())

    quadmeshes = [tup[3] for tup in hist2d_outputs]
    for quadmesh in quadmeshes:
        quadmesh.set_norm(hcnorm)

    def ev2nm(ev):
        """Helper function to convert eV to nm of wavelength

        Args:
            ev (ArrayLike): Float data of energy transitions in units of eV

        Returns:
            ArrayLike: The associated nm wafelength
        """
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
    axs['t1'].tick_params(axis="x", labelbottom=True)
    axs['sg'].tick_params(axis="x", labelbottom=True)

    secax = axs['sg'].secondary_xaxis('top', functions=(ev2nm, ev2nm))
    secax.set_xticks(
        [10, 25, 35, 40, 50, 75, 100, 125, 150, 200, 250, 300, 350, 400, 500, 750, 1000]
    )
    secax.tick_params(axis='x', rotation=45, labelsize='small')
    for l in secax.get_xticklabels():
        l.set_horizontalalignment('left')
        l.set_verticalalignment('bottom')
    secax.set_xlabel(r'$\lambda$ / nm')

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

            cb_spec.ax.axvline(t, c='white', linewidth=2)

    hcscale = mpl.cm.ScalarMappable(norm=hcnorm, cmap=magma_rw)
    axs['cb_hist'].figure.colorbar(hcscale, cax=axs['cb_hist'], location='bottom')
    axs['cb_hist'].set_xlabel('# data points')

    # axs['se'].set_title(
    #    r"$\uparrow$ground state" + "\n" + r"$\downarrow$excited state absorption"
    # )
    axs['t1'].set_xlabel(r'$\Delta E$ / eV')
    axs['sg'].set_xlabel(r'$\Delta E$ / eV')

    legend_lines, legend_labels = zip(
        *[
            (Line2D([0], [0], color='k', linestyle='-', linewidth=2), "$S_1/S_0$"),
            # (Line2D([0], [0], color='k', linestyle='--', linewidth=0.5), "$S_2/S_0$"),
        ]
    )
    axs['sg'].legend(legend_lines, legend_labels, fontsize='small')

    return axs
