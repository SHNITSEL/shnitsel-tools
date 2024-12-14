import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from . import (
  postprocess,
  xrhelpers)

def pca_line_plot(
  pca_res,
  hue,
  shadow=None,
  hue_label=None,
  hue_cmap=None,
  hue_qualitative=False,
  shadow_labels: dict|None =None,
  shadow_cmap=None,
  snips=None
):
    path_effects = mpl.patheffects
    LineCollection = mpl.collections.LineCollection
    BoundaryNorm = mpl.colors.BoundaryNorm
    ListedColormap = mpl.colors.ListedColormap
    mpatches = mpl.patches


    # DataArray.min() produces a DataArray, which range() doesn't like
    try:
        shadow=shadow.values 
    except AttributeError: pass # if shadow not a DataArray...

    if hue_label is None:
        hue_label = '{} / {}'.format(
            hue.attrs.get('unitdim'),
            hue.attrs.get('unit'))
    if hue_cmap is None:
        hue_cmap = 'binary_r'
    if shadow_labels is None:
        # NB. the keys follow SHARC's 1-based active state numbering convention.
        shadow_labels = {1: 'S$_0$', 2: 'S$_1$', 3: 'S$_2$'}
    if shadow_cmap is None:
        shadow_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            None, plt.cm.Pastel1(range(0,3)), 2)

    x = pca_res[:, 0]
    y = pca_res[:, 1]

    # Create a set of line segments so that we can color them individually
    # This creates the points as an N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    if snips is not None:
        segments = np.delete(segments, snips, axis=0)

    fig, axs = plt.subplots(1, 1)

    if shadow is not None:
        # Background shading:
        # Create a continuous norm to map from data points to colors
        # normed = plt.Normalize(shadow.min(), shadow.max())
        normed = plt.Normalize(0, 3)
        
        shadelc = LineCollection(
            segments, cmap=shadow_cmap, capstyle='round', norm=normed,)
            # path_effects=[path_effects.Stroke(joinstyle="round")])
        
        # Set the values used for colormapping
        shadelc.set_array(shadow)
        shadelc.set_linewidth(10)
        # lc.set_alpha(0.3)
        shading = axs.add_collection(shadelc)

        # something for plt.legend to get its colours from
        proxy_artists = [
            mpatches.Patch(
                color=shadow_cmap(i), label=shadow_labels[i])
            for i in range(shadow.min(), shadow.max() + 1)]

    # Main line:
    normed = plt.Normalize(hue.min(), hue.max())
    lc = LineCollection(
        segments,
        cmap=hue_cmap,
        norm=normed)
    lc.set_array(hue)
    lc.set_linewidth(2)
    line = axs.add_collection(lc)

    xpad = (x.max() - x.min())/50
    axs.set_xlim(x.min() - xpad, x.max() + xpad)
    ypad = (y.max() - y.min())/50
    axs.set_ylim(y.min() - ypad, y.max() + ypad)

    axs.set_xlabel('First principal component')
    axs.set_ylabel('Second principal component')

    if hue_qualitative:
        axs.legend(handles=line, loc='best')
    else:
        fig.colorbar(line).set_label(hue_label)

    # fig.colorbar(shading, ticks=range(0,3)).set_label('Active state')
    
    if shadow is None:
        legend = None
    else:
        legend = axs.legend(handles=proxy_artists, loc='best')
    # plt.close() # don't display duplicates, Jupyter!

    return fig, axs, legend

def pca_scatter_plot():
    ...

def timeplot(da, *, hue=None, time_axis='ts', delta_t=None, **kwargs):
    import seaborn as sns
    
    assert time_axis in {'ts', 'time'}
    xlabel = None

    default_delta_t = 0.5

    if delta_t is None:
        delta_t = da.attrs.get('delta_t', default_delta_t)

    df = (da
      .transpose('frame', ...)
      .to_pandas()
      .reset_index(['trajid', time_axis])
      .melt(id_vars=['trajid', time_axis])
    )

    xlabel = '$t$ / fs' if time_axis == 'time' else time_axis
    
    ax = sns.lineplot(df, x=time_axis, y='value', hue=hue, **kwargs)
    ax.set(xlabel=xlabel)
    return ax

def timeplot_interstate(da, delta_t=None, renamer=None, **kwargs):
    default_delta_t = 0.5

    if delta_t is None:
        delta_t = da.attrs.get('delta_t', default_delta_t)

    if 'statecomb' not in da.dims:
        da = postprocess.subtract_combinations(da, 'state', labels=True)
    da = postprocess.keep_norming(da)
    
    if 'statecomb' not in da.indexes:
        da = da.set_xindex(['from', 'to'])

    da = xrhelpers.flatten_midx(da, 'statecomb', renamer=renamer)

    # Does this really fall under this function's remit?
    da = postprocess.ts_to_time(da)
    
    return timeplot(da, hue='statecomb', time_axis='time', **kwargs)