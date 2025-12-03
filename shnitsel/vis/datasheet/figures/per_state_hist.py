from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
import numpy as np

from shnitsel.core.typedefs import PerState

from .common import figaxs_defaults, centertext, symbols
from .hist import truncate_from_above


@figaxs_defaults(mosaic=[['energy', 'forces', 'dip_perm']], scale_factors=(1, 1 / 5))
def plot_per_state_histograms(
    per_state: PerState,
    axs: dict[str, Axes] | None = None,
    fig: Figure | SubFigure | None = None,
) -> dict[str, Axes]:
    """Function to plot the per-state energy, forces and permanent dipole histogram plots.

    Args:
        per_state (PerState): A dataset with per-state observable data.
        axs (dict[str, Axes] | None, optional): The map of subplot-axes. Keys identify the subplots (`energy`, `forces`, `dip_perm`) and the values are the axes to plot the subplot to. Defaults to None.
        fig (Figure | SubFigure | None, optional): Figure to generated axes from. Defaults to None.

    Returns:
        dict[str, Axes]: The axes dictionary after plotting.
    """
    assert axs is not None, "Could not obtain axes for plotting the graphs."
    for quantity in ['energy', 'forces', 'dip_perm']:
        ax = axs[quantity]
        if quantity not in per_state:
            centertext("No %s data" % symbols.get(quantity, quantity), ax)
            continue

        for state, data in per_state.groupby('state'):
            c = data['_color'].item()
            counts, edges, _ = ax.hist(
                truncate_from_above(data[quantity].squeeze().values, bins=100),
                color=c,
                alpha=0.2,
                bins=100,
            )
            ax.plot((edges[1:] + edges[:-1]) / 2, counts, c=c, lw=0.5)
            idxmax = np.argmax(counts)
            ax.text(
                edges[[idxmax, idxmax + 1]].mean(),
                counts[idxmax],
                r"$S_%d$" % (state - 1),
                c=c,
            )

        long_name = per_state[quantity].attrs.get('long_name')
        units = per_state[quantity].attrs.get('units')
        axs[quantity].set_xlabel(rf'{long_name} / {units}')

    # for quantity in ['forces', 'dip_perm']:
    #     axs[quantity].
    axs['energy'].set_ylabel('# points')
    return axs
