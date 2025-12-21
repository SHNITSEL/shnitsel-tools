import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    """Draw multicoloured lines, and optionally
    larger "shadow" lines behind them

    Parameters
    ----------
    pca_res
        Data to plot the main line
    hue
        Data according to which to select colour
    shadow
        Data according to which to select the colours
        of the shadow
    hue_label
        Labels to describe different colours of the main line
    hue_cmap
        A colormap for the main line
    hue_qualitative
        Whether the colours of the main line are qualitative,
        by default False
    shadow_labels
        Labels to describe different colours of the main line
    shadow_cmap
        A colormap for the main line
    snips
        Indices indicating which points should *not* be connected;
        for instance, the end of one trajectory should probably
        not be connected to the next.

    Returns
    -------
    fig
        The matplotlib ``Figure`` object used for plotting
    ax
        The matplotlib ``Axes`` object used for plotting
    legend
        The object returned by ``matplotlib.pyplot.legend``
        for the main line
    """
    LineCollection = mpl.collections.LineCollection
    mpatches = mpl.patches


    # DataArray.min() produces a DataArray, which range() doesn't like
    try:
        shadow=shadow.values 
    except AttributeError: pass # if shadow not a DataArray...

    if hue_label is None:
        hue_label = '{} / {}'.format(
            hue.attrs.get('unitdim'),
            hue.attrs.get('units'),
        )
    if hue_cmap is None:
        hue_cmap = 'binary_r'
    if shadow_labels is None:
        # NB. the keys follow SHARC's 1-based active state numbering convention.
        shadow_labels = {1: 'S$_0$', 2: 'S$_1$', 3: 'S$_2$'}
    if shadow_cmap is None:
        shadow_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            None, plt.cm.Pastel1(range(0, 3)), 3
        )

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

    fig, ax = plt.subplots(1, 1)

    if shadow is not None:
        # Background shading:
        # Create a continuous norm to map from data points to colors
        # normed = plt.Normalize(shadow.min(), shadow.max())
        normed = plt.Normalize(0, 3)
        
        shadelc = LineCollection(
            segments, cmap=shadow_cmap, capstyle='round', norm=normed,)
            # path_effects=[mpl.path_effects.Stroke(joinstyle="round")])
        
        # Set the values used for colormapping
        shadelc.set_array(shadow)
        shadelc.set_linewidth(10)
        # lc.set_alpha(0.3)
        ax.add_collection(shadelc)

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
    line = ax.add_collection(lc)

    xpad = (x.max() - x.min())/50
    ax.set_xlim(x.min() - xpad, x.max() + xpad)
    ypad = (y.max() - y.min())/50
    ax.set_ylim(y.min() - ypad, y.max() + ypad)

    ax.set_xlabel('First principal component')
    ax.set_ylabel('Second principal component')

    if hue_qualitative:
        ax.legend(handles=line, loc='best')
    else:
        fig.colorbar(line).set_label(hue_label)

    # fig.colorbar(shading, ticks=range(0,3)).set_label('Active state')
    
    if shadow is None:
        legend = None
    else:
        legend = ax.legend(handles=proxy_artists, loc='best')
    # plt.close() # don't display duplicates, Jupyter!

    return fig, ax, legend
