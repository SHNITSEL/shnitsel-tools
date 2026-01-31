import logging
import matplotlib as mpl
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, Normalize
from matplotlib.figure import Figure, SubFigure
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import LinAlgError

from shnitsel.filtering.state_selection import StateSelection
from shnitsel.vis.datasheet.figures.dip_trans_hist import (
    plot_spectra,
    # single_dip_trans_hist,
)

from shnitsel.data.dataset_containers import InterState

from ....core.typedefs import SpectraDictType
from ....units.definitions import energy

from .common import centertext, figaxs_defaults
from .hist import calc_truncation_maximum, create_marginals
from ...colormaps import magma_rw
from ....units.conversion import convert_energy
from scipy.stats import gaussian_kde


def plot_energy_histogram(
    inter_state: InterState,
    state_selection: StateSelection,
    ax: Axes | None = None,
    bins: int = 100,
    cmap: str | Colormap | None = None,
    cnorm: Normalize | None = None,
    mark_peaks: bool = False,
    rasterized: bool = False,
) -> Axes:
    """Create the alternative to spectra plots where the oscillator strength is not available.

    Parameters
    ----------
    inter_state : InterState
        Interstate data to get delta energy histograms from
    state_selection : StateSelection
        State selection object to limit the states included in plotting and to provide state names.
    ax : Axes, optional
        Axis object to plot into. If not provided, will be created.
    bins : int, optional
        Optional number of bins for histogram creation. Defaults to 100.
    cmap : str | Colormap, optional
        Optional specification of a desired colormap. Defaults to None.
    cnorm : Normalize, optional
        Optional specification of a colormap norm method. Defaults to None.
    mark_peaks : bool, optional
        Flag whether peaks should be clearly marked. Defaults to False.
    rasterized : bool, optional
        Flag to control whether the histogram plot should be rasterized to cut down on loading times and file sizes in complex plot environments.

    Returns
    -------
    Axes
        The axes object into which the graph was plotted
    """
    if ax is None:
        _, ax = plt.subplots(1, 1)
    ax.set_ylabel(r'$\rho_{\Delta E}$')
    ax.invert_xaxis()
    # linestyles = {t: ['-', '--', '-.', ':'][i]
    #               for i, t in enumerate(np.unique(list(zip(*spectra.keys()))[0]))}

    if not inter_state.has_variable('energy_interstate'):
        centertext(r"No $\Delta E$ data provided.", ax=ax)
        return ax

    for i, sc in enumerate(state_selection.state_combinations):
        sc_data = inter_state.dataset.sel(statecomb=sc)
        sc_color = state_selection.get_state_combination_color(sc)
        sc_label = state_selection.get_state_combination_tex_label(sc)

        xdata = sc_data['energy_interstate'].squeeze()
        xdata = np.abs(convert_energy(xdata, to=energy.eV))

        xmax = calc_truncation_maximum(xdata)
        xmin = -calc_truncation_maximum(-xdata)

        n, x = np.histogram(xdata, range=(xmin, xmax), bins=100, density=True)
        bin_centers = (x[:-1] + x[1:]) / 2.0
        try:
            distr = gaussian_kde(xdata)
        except LinAlgError as e:
            logging.debug(f"Caught LinAlg error during histogram creation :{e}")
            continue

        # ax.fill_between(data['energy'], data, alpha=0.5, color=c)
        ax.plot(
            x,
            distr(x),
            # linestyle=linestyles[t], c=dcol_inter[sc],
            # linestyle=linestyle,
            c=sc_color,
            linewidth=0.8,
        )
        if mark_peaks:
            try:
                peak_pos = n.argmax()
                peak_dE = bin_centers[peak_pos]
                peak_dens = n[peak_pos]

                ax.text(
                    peak_dE,
                    peak_dens,
                    f"{sc_label}",
                    fontsize='xx-small',
                )
            except Exception as e:
                logging.warning(f"{e}")
    _, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax)

    return ax


def single_trans_hist(
    interstate: InterState,
    xvariable: str,
    yvariable: str,
    state_labels: tuple[str, str],
    color: str,
    xlabel: str | None = None,
    ylabel: str | None = None,
    bins: int = 100,
    ax: Axes | None = None,
    plot_marginals: bool = True,
    cmap=None,
    cnorm=None,
    rasterized: bool = True,
):
    """Function to plot a single histogram of interstate soc data vs. energy gaps.

    Parameters
    ----------
    interstate : InterState
        Inter-state Dataset
    xvariable : str
        Name of the interstate variable to use for the x-axis.
    yvariable : str
        Name of the interstate variable to use for the y-axis.
    xlabel : str, optional
        Label for the x-axis. Defaults to None.
    ylabel : str, optional
        Label for the y-axis. Defaults to None.
    sc_label : str
        Label to use for the state combination.
    state_labels : tuple[str,str]
        Labels for the individual states.
    color : str
        Color for the histogram of this state combination.
    bins : int, optional
        Number of bins for the histogram. Defaults to 100.
    ax : Axes, optional
        Axes object to plot into. Defaults to None.
    plot_marginals : bool, optional
        Flag to include marginal histograms. Defaults to True, meaning marginal histograms will be included.
    cmap : str, optional
        Colormap to use. Defaults to None.
    cnorm : str, optional
        Norming method to apply to the colormap. Defaults to None.
    rasterized : bool, optional
        Flag to control whether the histogram plot should be rasterized to cut down on loading times and file sizes in complex plot environments.

    Returns
    -------
        The result of ax.hist2d will be returned. 
    None 
        is returned if data is missing
    """
    if ax is None:
        _, ax = plt.subplots(1, 1)
    if cmap is None:
        cmap = magma_rw
    # TODO: FIXME: Merge with dip_trans_plot and generalize
    if xvariable not in interstate.data_vars or yvariable not in interstate.data_vars:
        # print("No SOC plot for missing data")
        # print(interstate.data_vars.keys())
        return None
    xdata = interstate.dataset[xvariable].squeeze()
    if xvariable == 'energy_interstate':
        # We expect energies in eV for the energy delta plot
        xdata = convert_energy(xdata, to=energy.eV)
        xdata = np.abs(convert_energy(xdata, to=energy.eV))

    # We need the normed transition soc
    ydata = interstate.dataset[yvariable].squeeze()
    # print(f"{xdata=}")
    # print(f"{ydata=}")

    xmax = calc_truncation_maximum(xdata)
    ymax = calc_truncation_maximum(ydata)
    # Get lower bounds dynamically as well.
    xmin = -calc_truncation_maximum(-xdata)
    ymin = -calc_truncation_maximum(-ydata)

    if plot_marginals:
        axx, axy = create_marginals(ax)
        axx.hist(
            xdata, range=(xmin, xmax), color=color, bins=bins, rasterized=rasterized
        )
        axy.hist(
            ydata,
            range=(ymin, ymax),
            orientation='horizontal',
            color=color,
            bins=bins,
            rasterized=rasterized,
        )

    hist2d_output = ax.hist2d(
        xdata,
        ydata,
        range=[(xmin, xmax), (ymin, ymax)],
        bins=bins,
        cmap=cmap,
        norm=cnorm,
        rasterized=rasterized,
    )

    if ylabel is not None:
        ax.set_ylabel(
            ylabel % (state_labels[0], state_labels[1])
        )  # r"$\|\mathbf{\mu}_{%d,%d}\|_2$" % (sc[0], sc[1]))

    if xlabel is not None:
        ax.set_xlabel(xlabel)

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


def single_soc_trans_hist(
    interstate: InterState,
    sc_label: str,
    state_labels: tuple[str, str],
    color: str,
    bins: int = 100,
    ax: Axes | None = None,
    plot_marginals: bool = True,
    cmap=None,
    cnorm=None,
    rasterized: bool = True,
):
    """Function to plot a single histogram of interstate data.

    Used for both delta_E vs. \\mu and delta_E vs SOC plots.

    Parameters
    ----------
    interstate : InterState
            Inter-state Dataset
    sc_label : str
            Label to use for the state combination.
    state_labels : tuple[str,str]
            Labels for the individual states.
    color : str
            Color for the histogram of this state combination.
    bins : int, optional
            Number of bins for the histogram. Defaults to 100.
    ax : Axes, optional
            Axes object to plot into. Defaults to None.
    plot_marginals : bool, optional
            Flag to include marginal histograms. Defaults to True, meaning marginal histograms will be included.
    cmap : str, optional
            Colormap to use. Defaults to None.
    cnorm : str, optional
            Norming method to apply to the colormap. Defaults to None.
    rasterized : bool, optional
            Flag to control whether the histogram plot should be rasterized to cut down on loading times and file sizes in complex plot environments.

    Returns
    -------
    The result of ax.hist2d will be returned. 
    None is returned if data is missing
    """
    return single_trans_hist(
        interstate=interstate,
        xvariable='energy_interstate',
        yvariable='socs_norm',
        ylabel=r"$\|SOC_{%s,%s}\|_2$",
        state_labels=state_labels,
        color=color,
        bins=bins,
        ax=ax,
        cmap=cmap,
        cnorm=cnorm,
        rasterized=rasterized,
    )


def single_dip_trans_hist(
    interstate: InterState,
    sc_label: str,
    state_labels: tuple[str, str],
    color: str,
    bins: int = 100,
    ax: Axes | None = None,
    plot_marginals: bool = True,
    cmap=None,
    cnorm=None,
    rasterized: bool = True,
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
    color : str
        Color for the histogram of this state combination.
    bins : int, optional
        Number of bins for the histogram. Defaults to 100.
    ax : Axes, optional
        Axes object to plot into. Defaults to None.
    plot_marginals : bool, optional
        Flag to include marginal histograms. Defaults to True, meaning marginal histograms will be included.
    cmap : str, optional
        Colormap to use. Defaults to None.
    cnorm : str, optional
        Norming method to apply to the colormap. Defaults to None.
    rasterized : bool, optional
        Flag to control whether the histogram plot should be rasterized to cut down on loading times and file sizes in complex plot environments.

    Returns
    -------
    The result of ax.hist2d will be returned
    """
    return single_trans_hist(
        interstate=interstate,
        xvariable='energy_interstate',
        yvariable='dip_trans_norm',
        ylabel=r"$\|\mathbf{\mu}_{%s,%s}\|_2$",
        state_labels=state_labels,
        color=color,
        bins=bins,
        ax=ax,
        cmap=cmap,
        cnorm=cnorm,
        rasterized=rasterized,
    )


def plot_soc_or_dip_trans_histograms(
    inter_state: InterState,
    state_selection: StateSelection,
    axs: list[Axes] | None = None,
    cnorm: str | None = None,
) -> list:
    """Function to plot all relevant histograms for the provided inter_state data

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
    for i, (sc_, data) in enumerate(inter_state.dataset.groupby('statecomb')):
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
        # print(data)
        # print(data.dip_trans_norm)
        # print(data.socs_norm)
        color = state_selection.get_state_combination_color(sc_)
        if 'dip_trans_norm' in data and data.dip_trans_norm.max() > 1e-9:
            # print("Opting for dip trans norm plot.")
            tmp_res = single_dip_trans_hist(
                InterState(direct_interstate_data=data),
                sc_label,
                state_labels,
                color=color,
                ax=ax,
                cnorm=cnorm,
            )
            if tmp_res is not None:
                hist2d_outputs.append(tmp_res)
        elif 'socs_norm' in data and data.socs_norm.max() > 1e-9:
            # print("Opting for socs norm plot.")
            tmp_res = single_soc_trans_hist(
                InterState(direct_interstate_data=data),
                sc_label,
                state_labels,
                color=color,
                ax=ax,
                cnorm=cnorm,
            )
            if tmp_res is not None:
                hist2d_outputs.append(tmp_res)
        else:
            centertext(r"No $\mathbf{\mu}_{ij}$ or $SOC$ data", ax, clearticks='xy')
            continue

        selected_scs += 1
    for s in range(selected_scs, len(axs)):
        centertext(r"No $\mathbf{\mu}_{ij}$ or $SOC$ data", axs[s], clearticks='xy')

    return hist2d_outputs


@figaxs_defaults(
    mosaic=[['sg'], ['t0'], ['t1'], ['se'], ['t2'], ['legend_spec'], ['cb_hist']],
    scale_factors=(1 / 3, 4 / 5),
    height_ratios=([1] * 5) + ([0.2] * 1) + ([0.1] * 1),
)
def plot_separated_spectra_and_soc_dip_hists(
    inter_state: InterState,
    spectra_groups: tuple[SpectraDictType, SpectraDictType],
    state_selection: StateSelection,
    fig: Figure | SubFigure | None = None,
    axs: dict[str, Axes] | None = None,
    cb_spec_vlines: bool = True,
    current_multiplicity: int | None = None,
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
    current_multiplicity : int, optional
        Can denote the current multiplicity to change the way we plot the ground state and/or excited state transitions.

    Returns
    -------
    dict[str, Axes]
        The axes dict after plotting to it.
    """
    assert axs is not None, "Could not acquire axes for plotting"
    ground, excited = spectra_groups
    scnorm = plt.Normalize(
        inter_state.dataset.time.min(), inter_state.dataset.time.max()
    )
    scmap = plt.get_cmap('turbo')
    # scscale = mpl.cm.ScalarMappable(norm=scnorm, cmap=scmap)
    non_degenerate_selection = state_selection.non_degenerate()

    # print(state_selection.states, state_selection.state_combinations)
    # print(non_degenerate_selection.states, non_degenerate_selection.state_combinations)

    time_unit = inter_state.dataset.time.attrs.get('units', 'fs')

    times = list(set([tup[0] for tup in ground]))
    times.sort()

    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']

    times_styles = {t: linestyles[i % len(linestyles)] for i, t in enumerate(times)}

    selection_ground_states = state_selection.ground_state_transitions()

    hist2d_outputs = []
    if current_multiplicity is None or current_multiplicity == 1 and len(ground) > 0:
        # Only plot ground state spectra in singlet mode
        # ground-state spectra and histograms
        plot_spectra(
            ground,
            ax=axs['sg'],
            lim_num_sc=2,
            state_selection=non_degenerate_selection,
            cnorm=scnorm,
            cmap=scmap,
        )
        # TODO: FIXME: Think about how to make the state transitions identifyable

        # if current_multiplicity == 1:
        #     legend_lines, legend_labels = zip(
        #         *[
        #             (
        #                 Line2D([0], [0], color='k', linestyle='-', linewidth=0.5),
        #                 "$S_1/S_0$",
        #             ),
        #             (
        #                 Line2D([0], [0], color='k', linestyle='--', linewidth=0.5),
        #                 "$S_2/S_0$",
        #             ),
        #         ]
        #     )
        #     axs['sg'].legend(legend_lines, legend_labels, fontsize='x-small')
    else:
        plot_energy_histogram(
            inter_state=inter_state,
            state_selection=selection_ground_states,
            ax=axs['sg'],
        )

    if len(selection_ground_states.state_combinations) > 1:
        selaxs = [axs['t1'], axs['t0']]
    else:
        selaxs = [axs['t1']]
        centertext(r"No $\mathbf{\mu}_{ij}$ or $SOC$ data", axs['t0'], clearticks='xy')

    res = plot_soc_or_dip_trans_histograms(
        inter_state,  # inter_state_sel,
        axs=selaxs,
        state_selection=selection_ground_states,
    )
    if res is not None:
        hist2d_outputs += res

    # excited-state spectra and histograms
    if len(excited) >= 1:
        plot_spectra(
            excited,
            ax=axs['se'],
            lim_num_sc=1,
            state_selection=non_degenerate_selection,
            cnorm=scnorm,
            cmap=scmap,
        )
        if current_multiplicity is None or current_multiplicity == 1:
            # Plot an excited state transition in the singlet case
            res = plot_soc_or_dip_trans_histograms(
                inter_state,  # inter_state.isel(statecomb=[2]),
                axs=[axs['t2']],
                state_selection=state_selection.excited_state_transitions(),
            )
            if res is not None:
                hist2d_outputs += res
        else:
            # Plot an energy histogram of non-permitted transitions in the higher-order case.
            plot_energy_histogram(
                inter_state=inter_state,
                state_selection=non_degenerate_selection.excited_state_transitions(),
                ax=axs['t2'],
            )
    else:
        plot_energy_histogram(
            inter_state=inter_state,
            state_selection=non_degenerate_selection.excited_state_transitions(),
            ax=axs['se'],
        )
        centertext(r"No $\mathbf{\mu}_{ij}$ or $SOC$ data", axs['t2'], clearticks='xy')

    hists = np.array([tup[0] for tup in hist2d_outputs])

    if len(hists) > 0:
        hcnorm = plt.Normalize(hists.min(), hists.max())

        quadmeshes = [tup[3] for tup in hist2d_outputs]
        for quadmesh in quadmeshes:
            quadmesh.set_norm(hcnorm)
        hcscale = mpl.cm.ScalarMappable(norm=hcnorm, cmap=magma_rw)
        axs['cb_hist'].figure.colorbar(hcscale, cax=axs['cb_hist'], location='bottom')
    else:
        axs['cb_hist'].set_axis_off()

    def ev2nm(delta_E):
        """Helper function to convert delta energy of a transition in eV to a nanometer wavelength

        Parameters
        ----------
        delta_E : float-like
            Energy delta in electron Volts (eV) units.

        Returns
        -------
        float-like
            Resulting wavelength in nm.
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

    # for lax in ['legend_spec', 'cb_hist']:
    #     axs[lax].get_yaxis().set_visible(False)

    # axs['cb_hist'].get_yaxis().set_visible(False)

    x_pos = list(range(0, len(times)))

    x_ticks = []
    x_labels = []
    for x, time in zip(x_pos, times):
        x_ticks.append(x)
        x_labels.append(f"{time:.1f}")

        axs['legend_spec'].plot(
            [x - 0.4, x + 0.4], [1, 1], color='k', ls=times_styles[time]
        )

    if len(x_ticks) > 0:
        axs['legend_spec'].get_yaxis().set_visible(False)
        # axs['legend_spec'].set_yticks([], [])
        axs['legend_spec'].set_xlabel(f'$t$/{time_unit}')
        axs['legend_spec'].set_xlim((-0.5, len(times) - 0.5))
        axs['legend_spec'].spines['top'].set_visible(False)
        axs['legend_spec'].spines['right'].set_visible(False)
        axs['legend_spec'].spines['bottom'].set_visible(False)
        axs['legend_spec'].spines['left'].set_visible(False)
        axs['legend_spec'].set_xticks(x_ticks, x_labels)
        # axs['legend_spec'].get_xaxis().set_visible(True)
    else:
        axs['legend_spec'].set_axis_off()

    # cb_spec = axs['cb_spec'].figure.colorbar(
    #     scscale,
    #     cax=axs['cb_spec'],
    #     location='bottom',
    #     extend='both',
    #     extendrect=True,
    # )
    # axs['cb_spec'].set_xlabel('time / fs')
    # if cb_spec_vlines:
    #     for t in times:
    #         lo, hi = scscale.get_clim()  # lines at these points don't show
    #         if t == lo:
    #             t += t / 100  # so we shift them slightly
    #         elif t == hi:
    #             t -= t / 100

    #         cb_spec.ax.axvline(t, c='white', linewidth=0.5)

    axs['cb_hist'].get_yaxis().set_visible(False)
    axs['cb_hist'].set_xlabel('# data points')

    axs['legend_spec'].tick_params(axis="x", labelbottom=True)
    axs['cb_hist'].tick_params(axis="x", labelbottom=True)

    axs['se'].set_title(
        r"$\uparrow$ground state" + "\n" + r"$\downarrow$excited state absorption"
    )
    axs['t2'].set_xlabel(r'$\Delta E$ / eV')

    return axs


@figaxs_defaults(
    # mosaic=[['sg'], ['t0'], ['t1'], ['se'], ['t2'], ['cb_spec'], ['cb_hist']],
    mosaic=[['sg', 't1'], ['cb_spec', 'cb_hist']],
    scale_factors=(4 / 5, 4 / 5),
    # height_ratios=([1] * 5) + ([0.1] * 2),
    height_ratios=([1]) + ([0.1]),
)
def plot_separated_spectra_and_soc_dip_hists_groundstate(
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
    ValueError
        Too few statecombs (expecting at least 2 states => 1 statecomb)

    Returns
    -------
        dict[str, Axes]: The labeled axes after plotting
    """
    assert axs is not None, "Could not acquire axes to plot to."
    ground, excited = spectra_groups
    times = [tup[0] for lst in spectra_groups for tup in lst]
    scnorm = plt.Normalize(
        inter_state.dataset.time.min(), inter_state.dataset.time.max() + 50
    )
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

    res = plot_soc_or_dip_trans_histograms(
        InterState(direct_interstate_data=inter_state.dataset.isel(statecomb=selsc)),
        axs=selaxs,
        state_selection=state_selection,
    )
    if res is not None:
        hist2d_outputs += res

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

        Parameters
        ----------
        ev :ArrayLike
            Float data of energy transitions in units of eV

        Returns
        -------
        ArrayLike
            The associated nm wafelength
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
