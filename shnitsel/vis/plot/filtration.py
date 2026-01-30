from typing import Any, Sequence, overload
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from shnitsel.clean.common import (
    cum_max_quantiles,
    true_upto,
    _filter_mask_from_criterion_mask,
)
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.multi_indices import ensure_unstacked
from shnitsel.data.tree.node import TreeNode
from shnitsel.vis.support.multi_plot import MultiPlot

shnitsel_blue = (44 / 255, 62 / 255, 80 / 255)  # '#2c3e50'
shnitsel_yellow = '#C4A000'
shnitsel_magenta = '#7E5273'
text_color = '#fff'
text_backgroundcolor = (0, 0, 0, 0.2)


@overload
def check_thresholds(
    ds_or_da: TreeNode[Any, xr.Dataset | xr.DataArray | ShnitselDataset],
    quantiles: Sequence[float] | None = None,
) -> MultiPlot: ...
@overload
def check_thresholds(
    ds_or_da: xr.Dataset | xr.DataArray | ShnitselDataset,
    quantiles: Sequence[float] | None = None,
) -> Axes: ...


def check_thresholds(
    ds_or_da: xr.Dataset
    | xr.DataArray
    | ShnitselDataset
    | TreeNode[Any, xr.Dataset | xr.DataArray | ShnitselDataset],
    quantiles: Sequence[float] | None = None,
) -> Axes | MultiPlot:
    """Display graphs illustrating
        1. how many trajectories meet each criterion throughout, and
        2. quantiles of cumulative maxima over time for each criterion, indicating at what times a given
           proportion has failed the criterion

    Parameters
    ----------
    ds_or_da : xr.Dataset | xr.DataArray | ShnitselDataset | TreeNode[Any, xr.Dataset | xr.DataArray | ShnitselDataset]
        Data to plot. Can be flat or hierarchical format.
    quantiles : Sequence[float] | None, optional
        Quantiles to display and mark on the right-hand graph, by default None.

    Returns
    -------
        The matplotlib ``Axes`` object of the plots
    """
    if isinstance(ds_or_da, TreeNode):
        # TODO: FIXME: We need to accumulate the data across flat groups first.
        return MultiPlot(ds_or_da.map_data(check_thresholds, quantiles=quantiles))
    else:
        if isinstance(ds_or_da, xr.DataArray):
            filtranda = ds_or_da.copy()
        elif "filtranda" in ds_or_da:
            filtranda = ds_or_da['filtranda'].copy()
        elif "energy_filtranda" in ds_or_da and "length_filtranda" not in ds_or_da:
            filtranda = ds_or_da['energy_filtranda'].copy()
        elif "length_filtranda" in ds_or_da and "energy_filtranda" not in ds_or_da:
            filtranda = ds_or_da['length_filtranda'].copy()
        elif "length_filtranda" in ds_or_da and "energy_filtranda" in ds_or_da:
            filtranda = xr.concat(
                [
                    ds_or_da['energy_filtranda'].copy(),
                    ds_or_da['length_filtranda'].copy(),
                ],
                dim='criterion',
            )
        else:
            raise ValueError(
                "Dataset provided to `check_thresholds()` has no filtranda data set."
            )

        if 'frame' in filtranda.dims:
            filtranda = filtranda.assign_coords(
                {'is_frame': ('frame', np.ones(filtranda.sizes['frame']))}
            )
            # Assuming filtranda is a stacked Dataset/DataArray, unstack it
            if hasattr(filtranda, 'drop_dims'):
                # DataArrays don't have this method
                filtranda = filtranda.drop_dims(['trajectory'], errors='ignore')
            filtranda = filtranda.unstack('frame').rename({'atrajectory': 'trajectory'})
            filtranda['is_frame'] = filtranda['is_frame'].fillna(0).astype(bool)

        calculated_quantile_positions = cum_max_quantiles(
            filtranda, quantiles=quantiles
        )

        if 'thresholds' in filtranda.coords:
            # TODO: This is too complicated. Why calculate quantiles first and and then calculate true_upto?
            # Extract the true_upto per filtranda and then get the quantiles from the set of `true_upto`.
            good_throughout = (
                (filtranda < filtranda['thresholds']) | (~filtranda['is_frame'])
            ).all('time')
            filtranda['proportion'] = (
                good_throughout.sum('trajectory') / good_throughout.sizes['trajectory']
            )
            calculated_quantile_positions['intercept'] = true_upto(
                calculated_quantile_positions < filtranda['thresholds'], 'time'
            )

        fig, axs = plt.subplots(
            calculated_quantile_positions.sizes['criterion'],
            2,
            sharex='col',
            sharey='row',
            layout='constrained',
            width_ratios=[1, 2],
        )
        fig.set_size_inches(6, 2 * calculated_quantile_positions.sizes['criterion'])
        for (title, data), ax in zip(
            calculated_quantile_positions.groupby('criterion'), axs[:, 1]
        ):
            if 'thresholds' in data.coords:
                threshold = data.coords['thresholds'].item()
                ax.axhline(threshold, c=shnitsel_yellow)
            else:
                threshold = None

            for qval, qdata in data.groupby('quantile'):
                qdata = qdata.squeeze(['criterion', 'quantile'])

                ax.fill_between(
                    qdata.coords['time'], qdata, fc=(0, 0, 0, 0.2), ec=(0, 0, 0, 0)
                )
                ax.text(
                    qdata['time'][-1], qdata[-1], f"{qval * 100} %", va='center', c='k'
                )

                if threshold is not None:
                    t_icept = qdata['intercept'].item()
                    ax.vlines(t_icept, 0, threshold, color=shnitsel_yellow, ls=':')
                    ax.text(
                        t_icept,
                        threshold,
                        f"{qval * 100} % <{t_icept}",
                        ha='right',
                        va='center',
                        c=text_color,
                        backgroundcolor=text_backgroundcolor,
                        rotation='vertical',
                        fontsize=6,
                    )

        for (title, data), ax in zip(filtranda.groupby('criterion'), axs[:, 0]):
            data = data.squeeze('criterion')
            ax.set_ylabel(title)
            # NOTE (thevro): groupby('trajectory').max() behaves strangely
            # for unstacked format, so use groupby.map() instead
            max_value_per_traj = data.groupby('trajectory').map(lambda x: x.max())
            ax.hist(
                max_value_per_traj,
                density=True,
                cumulative=True,
                orientation='horizontal',
                color=shnitsel_blue,
            )
            if 'thresholds' in data.coords:
                threshold = data.coords['thresholds'].item()
                ax.axhline(threshold, c=shnitsel_yellow)
                # ax.text(
                #     0.5,
                #     threshold,
                #     str(threshold),
                #     ha='center',
                #     va='bottom',
                #     c=text_color,
                #     backgroundcolor=text_backgroundcolor,
                # )
                ax.text(
                    0.5,
                    threshold,
                    f"{threshold}\n{data.coords['proportion'].item() * 100:.0f} %",
                    ha='center',
                    va='center',
                    c=text_color,
                    backgroundcolor=text_backgroundcolor,
                )

        axs[-1, 0].set_xlabel('cumulative density\nof per-traj maxima')
        axs[-1, 1].set_xlabel('time / fs')
        return axs


def validity_populations(ds_or_da, intersections: bool = True) -> Axes:
    """Display two plots showing
    1. how many trajectories meet criteria (or combinations thereof) up to a given time
    2. how many frames would remain if the ensemble were transected at a given time
    (see :py:func:`shnitsel.clean.transect`)

    Parameters
    ----------
    ds_or_da
        Data to plot
    intersections, optional
        whether to plot intersections of criteria (how many trajectories still meet criterion 1 AND criterion 2)
        or to consider criteria independently

    Returns
    -------
        The matplotlib ``Axes`` object of the plots
    """
    # First make sure we have a DataArray rather than a Dataset
    if hasattr(ds_or_da, 'data_vars'):
        filtranda = ds_or_da['filtranda'].copy()
    else:
        filtranda = ds_or_da.copy()

    # Then make sure the DataArray is in the unstacked aka. layered format

    # mask, _ = ensure_unstacked(mask)
    # TODO: Use this local reimplementation of `ensure_unstacked()` to rewrite ensure_unstacked;
    # see also in `check_thresholds()`
    if 'frame' in filtranda.dims:
        filtranda = filtranda.assign_coords(
            {'is_frame': ('frame', np.ones(filtranda.sizes['frame']))}
        )
        # Assuming filtranda is a stacked Dataset/DataArray, unstack it
        if hasattr(filtranda, 'drop_dims'):
            # DataArrays don't have this method
            filtranda = filtranda.drop_dims(['trajectory'], errors='ignore')
        filtranda = filtranda.unstack('frame').rename({'atrajectory': 'trajectory'})
        filtranda['is_frame'] = filtranda['is_frame'].fillna(0).astype(bool)

    # Drop thresholds to stop interference with total_population below
    mask = (filtranda < filtranda['thresholds']).drop_vars('thresholds')

    mask = (
        mask.to_dataset('criterion')
        .assign({'total_population': mask.coords['is_frame']})
        .to_dataarray('criterion')
    )

    # For populations, we need a cumulative mask
    mask = mask.cumprod('time')

    # FIXME: Once ensure_unstacked is fixed or replaced with as_layered, we should
    # use 'trajectory' in the following rather than 'atrajectory'
    counts = mask.sum('trajectory')
    means = counts.mean('time')
    if intersections:
        counts = (
            mask.sortby(means, ascending=False).cumprod('criterion').sum('trajectory')
        )
    else:
        counts = counts.sortby(means, ascending=False)
    fig, axs = plt.subplots(2, 1)
    fig.set_size_inches(6, 8)
    for criterion in counts.coords['criterion'].data:
        data = counts.sel(criterion=criterion)
        axs[0].plot(data.coords['time'], data, label=criterion)
        axs[1].plot(data.coords['time'], data * data.coords['time'], label=criterion)
    if intersections:
        order = counts.coords['criterion'].data
        labels = [order[0]] + ['AND ' + x for x in order[1:]]
        axs[0].legend(labels)
    else:
        axs[0].legend()
    axs[0].set_ylabel('# trajectories')
    axs[1].set_ylabel('# frames if transected at time')
    axs[1].set_xlabel('time / fs')
    return axs
