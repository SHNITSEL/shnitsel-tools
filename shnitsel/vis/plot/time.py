"""General timeplots -- use on anything with a time coordinate"""
# import shnitsel as st
from logging import warning
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from shnitsel._contracts import needs
from shnitsel.vis.plot.common import figax
from shnitsel.analyze.stats import time_grouped_confidence_interval

def _set_axes(data, ax=None):
    _, ax = figax(ax=ax)

    ylabel = data.attrs.get('long_name', data.name or '')
    if (yunits := data.attrs.get('units')):
        ylabel += f' / {yunits}'
    ax.set_ylabel(ylabel)
    xlabel = 'time'
    if (xunits := data['time'].attrs.get('units')):
        xlabel += f' / {xunits}'
    ax.set_xlabel(xlabel)
    return ax

def plot_single(data, ax=None, time_coord='time'):
    """Plot some property of a single trajectory over time

    Parameters
    ----------
    data
        Data to plot
    ax
        A matplotlib ``Axes`` onto which to plot;
        if not provided, one will be created.

    Returns
    -------
        The ``Axes`` object used

    Raises
    ------
    ValueError
        If ``data`` has both 'state' and 'statecomb' dimensions.
    """
    if 'state' in data.dims and 'statecomb' in data.dims:
        raise ValueError("data shouldn't have both `state` and `statecomb` dimensions")
    if 'state' in data.dims:
        groupby='state'
        coord_name = 'state_names'
    elif 'statecomb' in data.dims:
        groupby='statecomb'
        coord_name = 'statecomb_names'
    else:
        data = data.expand_dims('_dummy').assign_coords(_dummy=[''])
        groupby = '_dummy'
        coord_name = '_dummy'

    _, ax = figax(ax=ax)
    for _, sdata in data.groupby(groupby):
        sdata = sdata.squeeze(groupby)
        # c = sdata['_color'].item()
        line2d = ax.plot(sdata[time_coord], sdata, lw=0.5)  # , c=c)
        ax.text(
            sdata[time_coord][-1],
            sdata[-1],
            sdata.coords[coord_name].item(),
            va='center',
            c=line2d[0].get_color(),
        )
    return _set_axes(data, ax)

def plot_ci(data, ax=None, time_coord='time'):
    """Plot some property of trajectories over time, aggregated as means and confidence intervals

    Parameters
    ----------
    data
        Data to plot
    ax
        A matplotlib ``Axes`` onto which to plot;
        if not provided, one will be created.

    Returns
    -------
        The ``Axes`` object used

    Raises
    ------
    ValueError
        If ``data`` has both 'state' and 'statecomb' dimensions.
    """
    _, ax = figax(ax=ax)

    if 'state' in data.dims and 'statecomb' in data.dims:
        raise ValueError("data shouldn't have both `state` and `statecomb` dimensions")
    elif 'state' in data.dims:
        dim='state'
        coord_name = 'state_names'
    elif 'statecomb' in data.dims:
        dim='statecomb'
        coord_name = 'statecomb_names'
    else:
        # TODO FIXME The expand_dims and squeeze steps shouldn't be necessary
        ci = time_grouped_confidence_interval(data.expand_dims('state')).squeeze('state')
        ax.fill_between(time_coord, 'upper', 'lower', data=ci, alpha=0.3)
        line2d = ax.plot(time_coord, 'mean', data=ci, lw=0.8)
        return _set_axes(data, ax)
        

    ci = time_grouped_confidence_interval(data)
    for _, sdata in ci.groupby(dim):
        sdata = sdata.squeeze(dim)
        ax.fill_between(time_coord, 'upper', 'lower', data=sdata, alpha=0.3)
        line2d = ax.plot(time_coord, 'mean', data=sdata, lw=0.8)
        ax.text(
            sdata[time_coord][-1],
            sdata['mean'][-1],
            sdata.coords[coord_name].item(),
            va='center',
            c=line2d[0].get_color(),
        )
    return _set_axes(data, ax)

def plot_many(data, ax=None, time_coord='time'):
    """Plot some property of trajectories over time as thin lines;
    state or statecomb is indicated by colour

    Parameters
    ----------
    data
        Data to plot
    ax
        A matplotlib ``Axes`` onto which to plot;
        if not provided, one will be created.

    Returns
    -------
        The ``Axes`` object used

    Raises
    ------
    ValueError
        If ``data`` has both 'state' and 'statecomb' dimensions.
    """
    _, ax = figax(ax=ax)

    if 'state' in data.dims and 'statecomb' in data.dims:
        raise ValueError("data shouldn't have both `state` and `statecomb` dimensions")
    elif 'state' in data.dims:
        dim='state'
        groupby = data.groupby('state')
        coord_name = 'state_names'
    elif 'statecomb' in data.dims:
        dim = 'statecomb'
        groupby = data.groupby('statecomb')
        coord_name = 'statecomb_names'
    else:
        dim = 'tmp'
        groupby = [(None, data)]
        for _, traj in data.groupby('atrajectory'):
            ax.plot(traj[time_coord], traj, lw=0.5, c='k')
        return _set_axes(data, ax)
    
    colors = iter(plt.get_cmap('tab10').colors)
    for _, sdata in groupby:
        sdata = sdata.squeeze(dim)
        label = sdata.coords[coord_name].item()
        c = next(colors)
        for _, traj in sdata.groupby('atrajectory'):
            ax.plot(traj[time_coord], traj, lw=0.5, label=label, c=c)
    # TODO: legend
    return _set_axes(data, ax)

def plot_shaded(data, ax, time_coord='time'):
    """Plot some property of trajectories over time, aggregated by using colour
    to show how many overlap at a given point

    Parameters
    ----------
    data
        Data to plot
    ax
        A matplotlib ``Axes`` onto which to plot;
        if not provided, one will be created.

    Returns
    -------
        The ``Axes`` object used.

    Raises
    ------
    ImportError
        If the ``datashader`` library is not installed.
    """
    try:
        import datashader as ds
    except ImportError as err:
        raise ImportError('plot_shaded requires the optional datashader dependency') from err
    try:
        import colorcet

        cmap = colorcet.bjy
    except ImportError:
        warning("colorcet package not installed; falling back on viridis cmap")
        cmap = plt.get_cmap('viridis')

    _, ax = figax(ax=ax)

    x = []
    y = []
    for _, traj in data.groupby('atrajectory'):
        x.append(traj.coords[time_coord].values)
        y.append(traj.values)
    df = pd.DataFrame({
        'x': pd.array(x, dtype='Ragged[float64]'),
        'y': pd.array(y, dtype='Ragged[float64]'),
    })
    cvs = ds.Canvas(plot_height=2000, plot_width=2000)
    agg = cvs.line(df, x='x', y='y', agg=ds.count(), line_width=5, axis=1)
    img = ds.tf.shade(agg, how='log', cmap=cmap)
    arr = np.array(img.to_pil())
    x0, x1 = agg.coords['x'].values[[0,-1]]
    y0, y1 = agg.coords['y'].values[[0, -1]]
    ax.imshow(arr, extent=[x0, x1, y0, y1], aspect='auto')
    return _set_axes(data, ax)

@needs(coords={"time"})
def timeplot(
    data: xr.DataArray,
    ax: plt.Axes | None = None,
    trajs: Literal['ci', 'shade', 'conv', None] = None,
    sep: bool = False,
    time_coord='time',
):
    """Plot some property of one or many trajectories over time,
    possibly aggregating over trajectories, and distinguishing
    different states/statecombs if applicable.

    Parameters
    ----------
    data
        Data to plot
    ax
        A matplotlib ``Axes`` onto which to plot;
        if not provided, one will be created.
    trajs
        How to aggregate trajectories, if at all:
            - ``None`` (default): do not aggregate
            - 'ci': aggregate by confidence interval
            - 'shade': use colour to represent overlap
                density using the datashader library
                (to produce what is sometimes called
                a hairplot)
    sep
        Whether to plot different states/statecombs
        separately; this will be done regardless when
        using ``trajs='shade'``.


    Returns
    -------
        The ``Axes`` object used

    Raises
    ------
    ValueError

        - If ``data`` has both 'state' and 'statecomb' dimensions.
        - If ``ax`` is passed when multiple ``Axes`` will be required
            as states/statecombs are to be plotted separately.
        - If ``trajs`` is set to a value other than ``None``, 'ci' or 'shade'
    NotImplementedError
        If ``trajs='conv'`` is used
    """
    if {'state', 'statecomb'}.issubset(data.dims):
        raise ValueError(
            "`data` should not have both 'state' and 'statecomb' dimensions"
        )
    state_dim = (
        'state'
        if 'state' in data.dims
        else 'statecomb'
        if 'statecomb' in data.dims
        else ''
    )

    if trajs in {'shade', 'conv'} and state_dim:
        sep = True

    if sep:
        if ax is not None:
            raise ValueError("Plotting multiple plots, so `ax` arg can not be used")
        nplots = data.sizes[state_dim]
        fig, axs = plt.subplots(1, nplots, layout='constrained', sharex=True)
        fig.set_size_inches(4 * nplots, 1.1414 * 4)
        res = []
        coord_name = state_dim + '_names'
        for (_, sdata), ax in zip(data.groupby(state_dim), axs):
            ax.set_title(sdata.coords[coord_name].item())
            sdata = sdata.squeeze(state_dim)
            res.append(
                timeplot(sdata, ax=ax, trajs=trajs, sep=False, time_coord=time_coord)
            )
        return res

    if 'trajid' not in data.coords and 'atrajectory' not in data.coords:
        assert trajs is None
        return plot_single(data, ax, time_coord=time_coord)
    if trajs == 'ci':
        return plot_ci(data, ax, time_coord=time_coord)
    elif trajs == 'shade':
        return plot_shaded(data, ax, time_coord=time_coord)
    elif trajs == 'conv':
        raise NotImplementedError(
            "Convolutions are not yet implemented here, "
            "please use xr_broaden_gauss manually"
        )
    elif trajs is None:
        return plot_many(data, ax, time_coord=time_coord)
    else:
        raise ValueError(
            f"`trajs` should be one of 'ci', 'shade' or None, rather than {trajs}"
        )
