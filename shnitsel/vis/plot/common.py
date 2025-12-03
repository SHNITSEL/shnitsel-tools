import io

import PIL
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure


def outlabel(ax, letter):
    fixedtrans = mpl.transforms.ScaledTranslation(
        -20 / 72, +7 / 72, ax.figure.dpi_scale_trans
    )
    transform = ax.transAxes + fixedtrans
    return ax.text(
        0.0,
        1.0,
        letter,
        transform=transform,
        va='bottom',
        fontweight='bold',
        bbox=dict(facecolor='0.9', edgecolor='none', pad=3.0),
    )


def inlabel(ax, letter):
    return ax.annotate(
        letter,
        xy=(1, 1),
        xycoords='axes fraction',
        xytext=(-1, -0.5),
        textcoords='offset fontsize',
        va='top',
        fontweight='bold',
        bbox=dict(facecolor='0.9', edgecolor='none', pad=3.0),
    )


def figax(
    fig: Figure | SubFigure | None = None,
    ax: Axes | None = None,
) -> tuple[Figure | SubFigure, Axes]:
    """
    Create figure and axes-object if an axes-object is not supplied.
    """
    if fig is None and ax is None:
        fig, ax = plt.subplots(1, 1)
    elif fig is None:
        assert ax is not None
        fig = ax.figure
    elif ax is None:
        ax = fig.subplots(1, 1)

    assert isinstance(fig, Figure) or isinstance(fig, SubFigure)
    assert isinstance(ax, Axes)
    return fig, ax


def extrude(x, y, xmin, xmax, ymin, ymax):
    # for extrusion, flip negative rays into quadrant 1
    if x < 0:
        xlim = -xmin # positive
        xsgn = -1
    else:
        xlim = xmax
        xsgn = 1
    if y < 0:
        ylim = -ymin # positive
        ysgn = -1
    else:
        ylim = ymax
        ysgn = 1
    # now extrude
    x2 = abs(ylim*x/y)  # try extruding till we reach the top
    if x2 <= xlim: # have we dropped off the right?
        y2 = ylim  # if not, go with this
    else:          # but if we would have dropped off the right
        x2 = xlim  # just go as far right as possible instead
        y2 = abs(xlim*y/x)
    return x2*xsgn, y2*ysgn


def mpl_imshow_png(ax, png, **imshow_kws):
    buffer = io.BytesIO()
    buffer.write(png)
    buffer.seek(0)
    img_array = np.array(PIL.Image.open(buffer))
    ax.axis('off')
    return ax.imshow(img_array, rasterized=True, **imshow_kws)
