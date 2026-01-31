from matplotlib.figure import Figure, SubFigure
from shnitsel._contracts import needs
from shnitsel.analyze.generic import keep_norming
from shnitsel.core.typedefs import InterState
from shnitsel.filtering.state_selection import StateSelection
from ....units.definitions import energy
from ....units.conversion import convert_energy

from .common import figaxs_defaults, centertext
from .hist import create_marginals, calc_truncation_maximum

from matplotlib.axes import Axes
import xarray as xr


@figaxs_defaults(mosaic=[['ntd'], ['nde']], scale_factors=(1 / 3, 1 / 3))
@needs(data_vars={'socs'})
def plot_socs_histograms(
    inter_state: InterState,
    hops_mask: xr.DataArray,
    state_selection: StateSelection,
    fig: Figure | SubFigure | None = None,
    axs: dict[str, Axes] | None = None,
) -> dict[str, Axes]:
    """Plot 2D histograms of NACS vs delta_E or dip_trans

    Parameters
    ----------
    inter_state : InterState
            The dataset containing inter-state data including NACs
    hop_idxs: 
        Argument to specify, which frames should be selected for the histograms.
    state_selection : StateSelection
        State selection object to limit the states included in plotting and to provide state names.
    fig : Figure | SubFigure, optional
        Unused figure provided to the plot. Consumed by the figaxs_defaults decorator.
    axs : dict[str, Axes]
        Axes objects to plot to with the respective keys of the plot. Defaults to None.

    Returns
    -------
    dict[str, Axes]
        The axes used for plotting indexed by the subfigure name
    """
    assert axs is not None, "No axes objects provided."

    hop_filter_data = inter_state.sel(frame=hops_mask)
    axs['nde'].set_ylabel(r'$\Delta E$ / eV')
    axs['nde'].minorticks_on()
    axs['nde'].set_xlabel(r"$\|\mathrm{NAC}_{i,j}\|_2$")
    axs['ntd'].tick_params(axis="x", labelbottom=False)

    if 'dip_trans' in hop_filter_data:
        axs['ntd'].set_ylabel(r"$\|\mathbf{\mu}_{i,j}\|_2$")
        axs['ntd'].minorticks_on()
        if 'dip_trans_norm' not in hop_filter_data:
            hop_filter_data['dip_trans_norm'] = keep_norming(hop_filter_data.dip_trans)

    if 'nacs' in hop_filter_data:
        if 'nacs_norm' not in hop_filter_data:
            hop_filter_data['nacs_norm'] = keep_norming(hop_filter_data.nacs)
    # print(hop_filter_data)

    def plot(label, yname, nacs_data):
        ax = axs[label]
        axx, axy = create_marginals(ax)
        bins = 100

        for i, (sc, data) in enumerate(nacs_data.groupby('statecomb')):
            if not state_selection.has_state_combination(sc):
                continue

            ydata = data[yname].squeeze()
            xdata = data['nacs_norm'].squeeze()
            # xmax = trunc_max(xdata)
            if yname == 'energy_interstate':
                ydata = convert_energy(ydata, to=energy.eV)

            ymax = calc_truncation_maximum(ydata)
            color = state_selection.get_state_combination_color(sc)
            # color = data['_color'].item()
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

    if 'energy_interstate' in inter_state:
        plot('nde', 'energy_interstate', hop_filter_data)

    if 'dip_trans' in inter_state:
        plot('ntd', 'dip_trans_norm', hop_filter_data)
    else:
        centertext(r"No $\mathbf{\mu}_{ij}$ data", axs['ntd'])
        # axs['ntd'].tick_params(axis='y', labelleft=False)
        axs['ntd'].get_yaxis().set_visible(False)
        axs['ntd'].get_xaxis().set_visible(False)

    return axs
