import logging
from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
import numpy as np

from shnitsel.data.dataset_containers import PerState
from shnitsel.filtering.state_selection import StateSelection

from .common import figaxs_defaults, centertext
from .hist import truncate_from_above

symbols = dict(
    energy=r"$E_i$", forces=r"$|\mathbf{F}_i|$", dip_perm=r"$|\mathbf{\mu}_i|$"
)


@figaxs_defaults(mosaic=[['energy', 'forces', 'dip_perm']], scale_factors=(1, 1 / 5))
def plot_per_state_histograms(
    per_state: PerState,
    state_selection: StateSelection,
    shape: tuple[int, int] | None = None,
    axs: dict[str, Axes] | None = None,
    fig: Figure | SubFigure | None = None,
) -> dict[str, Axes]:
    """Function to plot the per-state energy, forces and permanent dipole histogram plots.

    Parameters
    ----------
    per_state : PerState
        A dataset with per-state observable data.
    state_selection (StateSelection
        State selection object to limit the states included in plotting and to provide state names.
    axs : dict[str, Axes] | None, optional
        The map of subplot-axes. Keys identify the subplots (`energy`, `forces`, `dip_perm`) and the values are the axes to plot the subplot to. Defaults to None.
    shape : tuple[int, int] | None, optional
        Optional argument to reshape the per-state plots into another arrangement of rows and columns. First argument number of rows, second the number of columns.
    fig : Figure | SubFigure | None, optional
        Figure to generated axes from. Defaults to None.

    Returns
    -------
    dict[str, Axes]
        The axes dictionary after plotting.
    """
    assert axs is not None, "Could not obtain axes for plotting the graphs."

    if shape is not None:
        rows = shape[0]
        cols = shape[1]

        if rows > 1 or cols < 3:
            if rows * cols < 3:
                logging.warning(
                    "Specified shape does not accommodate 3 plots for state histograms."
                )
            else:
                for k, ax in list(axs.items()):
                    ax.remove()
                    del axs[k]
                total = rows * cols
                assert fig is not None, "Figure is required if shape is passed."
                ax_list = fig.subplots(rows, cols)
                ax_list = ax_list.flatten()
                for label, ax in zip(['energy', 'forces', 'dip_perm'], ax_list):
                    axs[label] = ax

                for pos in range(0, total, cols):
                    ax_list[pos].set_ylabel('# points')
        else:
            axs['energy'].set_ylabel('# points')

    for quantity in ['energy', 'forces', 'dip_perm']:
        ax = axs[quantity]
        if (
            not per_state.has_variable(quantity)
            or not per_state.dataset[quantity].notnull().any()
        ):
            centertext("No %s data" % symbols.get(quantity, quantity), ax)
            continue

        for state, data in per_state.dataset.groupby('state'):
            if not state_selection.has_state(state):
                continue

            if not data[quantity].notnull().any():
                continue

            color = state_selection.get_state_color(state)
            state_label = state_selection.get_state_tex_label(state)
            counts, edges, _ = ax.hist(
                truncate_from_above(data[quantity].squeeze().values, bins=100),
                color=color,
                alpha=0.2,
                bins=30,
            )
            ax.plot((edges[1:] + edges[:-1]) / 2, counts, c=color, lw=0.5)
            idxmax = np.argmax(counts)
            ax.text(
                edges[[idxmax, idxmax + 1]].mean(),
                counts[idxmax],
                f"${state_label}$",
                c=color,
            )

        long_name = symbols[quantity]  # per_state[quantity].attrs.get('long_name')
        units = per_state.data_vars[quantity].attrs.get('units')
        axs[quantity].set_xlabel(rf'{long_name} / {units}')

    # for quantity in ['forces', 'dip_perm']:
    #     axs[quantity].
    return axs
