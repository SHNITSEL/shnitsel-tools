import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure, SubFigure
import numpy as np
from shnitsel.data.dataset_containers import PerState
import xarray as xr

from shnitsel.filtering.state_selection import StateSelection
from shnitsel.vis.datasheet.figures.common import figaxs_defaults
from shnitsel.vis.plot.common import figax
from shnitsel.units.conversion import convert_energy
from shnitsel.units.definitions import energy
import scipy.stats as st


@figaxs_defaults(
    mosaic=[
        ['pc1'],
        ['pc2'],
    ],
    # scale_factors=(1 / 3, 4 / 5),
    # height_ratios=([1] * 5) + ([0.1] * 2),
)
def plot_energy_bands(
    per_state: PerState,
    pca_data: xr.DataArray,
    state_selection: StateSelection,
    hops_mask: xr.DataArray | None = None,
    fig: Figure | SubFigure | None = None,
    axs: dict[str, Axes] | None = None,
    colorbar_label: str | None = None,
    cmap: str | Colormap | None = None,
    cnorm: str | Normalize | None = None,
    cscale=None,
    band_kws=None,
    hops_kws=None,
) -> dict[str, Axes]:
    """Helper function to plot energy band graphs as a function of the principal components in pca_data.

    Parameters
    ----------
    per_state : PerState
        _description_
    pca_data : xr.DataArray
        _description_
    state_selection : StateSelection
        _description_
    hops_mask : xr.DataArray, optional
        _description_. Defaults to None.
    fig : Figure | SubFigure | None, optional
        _description_. Defaults to None.
    ax : dict[str, Axes], optional
        _description_. Defaults to None.
    colorbar_label : str | None, optional
        _description_. Defaults to None.
    cmap : str | Colormap | None, optional)
        _description_. Defaults to None.
    cnorm : str | Normalize | None, optional
        _description_. Defaults to None.
    cscale : _type_, optional)
        _description_. Defaults to None.
    band_kws : _type_, optional)
        _description_. Defaults to None.
    hops_kws : _type_, optional
        _description_. Defaults to None.

    Returns
    -------
    dict[str, Axes]
        _description_
    """
    assert axs is not None, "Could not acquire axes"
    # fig, ax = figax(fig=fig, ax=ax)
    band_kws = band_kws or {}
    # band_kws = {'alpha': 0.5, 's': 0.2, **band_kws}
    band_kws = {'alpha': 0.2, **band_kws}

    state_names = []
    state_colors = []

    for state_id in per_state.coords['state'].values:
        state_color = state_selection.get_state_color(state_id)

        state_names.append(f"${state_selection.get_state_tex_label(state_id)}$")
        state_colors.append(state_color)

    for pc_id, pc_ax in zip([0, 1], [axs["pc1"], axs["pc2"]]):
        for state_id, per_state_data in per_state.dataset.groupby('state'):
            state_color = state_selection.get_state_color(state_id)

            state_energy = per_state_data.energy
            state_energy = convert_energy(state_energy, to=energy.eV)

            state_pc1 = pca_data.isel(PC=pc_id)

            pc_centers = []
            pc_bounds = []

            for pc_key, pc_data in (
                state_energy.to_dataset()
                .assign({"pc": state_pc1})
                .groupby_bins('pc', bins=50)
            ):
                pc_mean = pc_data.energy.mean()
                num_points = len(pc_data.energy.values)
                if num_points > 3:
                    std_error = st.sem(pc_data.energy.values)
                    res_interval = st.t.interval(
                        0.95,
                        num_points - 1,
                        loc=pc_mean,
                        scale=std_error,
                    )
                    pc_min, pc_max = res_interval
                else:
                    pc_min = pc_data.energy.min()
                    pc_max = pc_data.energy.max()

                pc_centers.append(pc_data.pc.mean())
                pc_bounds.append([pc_min, pc_mean, pc_max])

            pc_bounds = np.array(pc_bounds)

            # Alternative scatter plot option
            # pc_ax.scatter(state_pc1, state_energy, c=state_color, **band_kws)
            pc_ax.fill_between(
                pc_centers,
                pc_bounds[:, 0],
                pc_bounds[:, 2],
                color=state_color,
                **band_kws,
            )

            pc_ax.plot(pc_centers, pc_bounds[:, 1], c=state_color)
            pc_ax.set_xlabel(f'PC{pc_id + 1}')
            pc_ax.set_ylabel(r'$E_{i}/eV$')

    num_states = len(state_names)
    cm = LinearSegmentedColormap.from_list(
        'state_color_map', state_colors, N=num_states
    )

    colorbar_axes = plt.colorbar(
        ScalarMappable(cmap=cm, norm=plt.Normalize(0, num_states - 1)),
        ticks=0.5 + np.arange(num_states) * (num_states - 1) / num_states,
        ax=[axs["pc1"], axs["pc2"]],
        label='State',
    )
    colorbar_axes.set_ticklabels(state_names)

    # if hops is not None:
    #     hops_kws = dict(s=0.5, c='limegreen') | (hops_kws or {})
    #     ax.scatter(hops.isel(PC=0), hops.isel(PC=1), **hops_kws)

    # TODO facilitate custom colorbar
    # fig.colorbar(cscale, ax=ax, label=colorbar_label, pad=0.02)

    # Alternative layout solution
    # d = make_axes_locatable(ax)
    # cax = d.append_axes("right", size="5%", pad="2%")
    # fig.colorbar(pc, cax=cax, label='dihedral')

    return axs
