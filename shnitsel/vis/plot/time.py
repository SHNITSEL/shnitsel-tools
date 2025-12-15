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

def set_axes(data, ax=None):
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

def plot_single(data, ax=None):
    if 'state' in data.dims and 'statecomb' in data.dims:
        raise ValueError("data shouldn't have both `state` and `statecomb` dimensions")
    if 'state' in data.dims:
        groupby='state'
        coord_name = 'state_names'
    if 'statecomb' in data.dims:
        groupby='statecomb'
        coord_name = 'statecomb_names'

    _, ax = figax(ax=ax)
    for _, sdata in data.groupby(groupby):
        sdata = sdata.squeeze(groupby)
        # c = sdata['_color'].item()
        line2d = ax.plot(sdata['time'], sdata, lw=0.5)#, c=c)
        ax.text(
            sdata['time'][-1],
            sdata[-1],
            sdata.coords[coord_name].item(),
            va='center',
            c=line2d[0].get_color())
    return set_axes(data, ax)

def plot_ci(data, ax=None):
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
        ax.fill_between('time', 'upper', 'lower', data=ci, alpha=0.3)
        line2d = ax.plot('time', 'mean', data=ci, lw=0.8)
        return set_axes(data, ax)
        

    ci = time_grouped_confidence_interval(data)
    for _, sdata in ci.groupby(dim):
        sdata = sdata.squeeze(dim)
        ax.fill_between('time', 'upper', 'lower', data=sdata, alpha=0.3)
        line2d = ax.plot('time', 'mean', data=sdata, lw=0.8)
        ax.text(
            sdata['time'][-1],
            sdata['mean'][-1],
            sdata.coords[coord_name].item(),
            va='center',
            c=line2d[0].get_color()
        )
    return set_axes(data, ax)

def plot_many(data, ax=None):
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
        for _, traj in data.groupby('trajid'):
            ax.plot(traj['time'], traj, lw=0.5, c='k')
        return set_axes(data, ax)
    
    colors = iter(plt.get_cmap('tab10').colors)
    for _, sdata in groupby:
        sdata = sdata.squeeze(dim)
        label = sdata.coords[coord_name].item()
        c = next(colors)
        for _, traj in sdata.groupby('trajid'):
            ax.plot(traj['time'], traj, lw=0.5, label=label, c=c)
    # TODO: legend
    # TODO: option/default of separate plot per state(comb)
    return set_axes(data, ax)

def plot_shaded(data, ax):
    # TODO: automatically make separate plot per state(comb)
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
    for _, traj in data.groupby('trajid'):
        x.append(traj.time.values)
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
    return set_axes(data, ax)

@needs(coords={"time"})
def timeplot(
    data: xr.DataArray,
    ax: plt.Axes | None = None,
    trajs: Literal['ci', 'shade', 'conv', None] = None
):
    if 'trajid' not in data.coords:
        assert trajs is None
        return plot_single(data, ax)
    if trajs == 'ci':
        return plot_ci(data, ax)
    elif trajs == 'shade':
        return plot_shaded(data, ax)
    elif trajs == 'conv':
        raise NotImplementedError(
            "Convolutions are not yet implemented here, "
            "please use xr_broaden_gauss manually"
        )
    elif trajs is None:
        return plot_many(data, ax)
    else:
        raise ValueError(
            f"`trajs` should be one of 'ci', 'shade' or None, rather than {trajs}"
        )
