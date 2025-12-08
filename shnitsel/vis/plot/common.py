import io

import PIL
import matplotlib as mpl
from matplotlib.image import AxesImage
from matplotlib.text import Annotation, Text
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure


def outlabel(ax: Axes, label: str) -> Text:
    """Helper function to put a label besides the graph in `ax` (outside the axes).

    Args:
        ax (Axes): The axes object to put the label next to.
        label (str): The label text to add.

    Returns:
        Text: The created mpl.Text object.
    """
    fixedtrans = mpl.transforms.ScaledTranslation(
        -20 / 72, +7 / 72, ax.figure.dpi_scale_trans
    )
    transform = ax.transAxes + fixedtrans
    return ax.text(
        0.0,
        1.0,
        label,
        transform=transform,
        va='bottom',
        fontweight='bold',
        bbox=dict(facecolor='0.9', edgecolor='none', pad=3.0),
    )


def inlabel(ax: Axes, label: str) -> Annotation:
    """Helper function to add a text label inside of the axes to `ax`.

    Args:
        ax (Axes): The axes object to put the label into.
        label (str): The label text to add.

    Returns:
        Annotation: The resulting mpl.Annotation object representing the inserted label.
    """
    return ax.annotate(
        label,
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

    Args:
        fig (Figure | SubFigure | None, optional): The optional figure to use as a basis for `ax` if the latter is not provided. Defaults to None.
        ax (Axes | None, optional): The axes object provided. Will be used to populate `fig` if provided.. Defaults to None.

    Returns:
        tuple[Figure | SubFigure, Axes]: A complete combination of figure and axes.
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


def extrude(
    x: float, y: float, xmin: float, xmax: float, ymin: float, ymax: float
) -> tuple[float, float]:
    """Calculate the endpoint of extrusion of the point (x,y) from point (0,0) until it intersects either x or y boundary.

    Args:
        x (float): x coordinate of the base point
        y (float): y coorindate of the base point
        xmin (float): Lower x limit
        xmax (float): Upper x limit
        ymin (float): Lower y limit
        ymax (float): Upper y limit

    Returns:
        tuple[float, float]: The position at the end of the extrusion, where the origin-ray through (x,y) intersects the boundary of the axes.
    """
    # TODO: Document
    # for extrusion, flip negative rays into quadrant 1
    if x < 0:
        xlim = -xmin  # positive
        xsgn = -1
    else:
        xlim = xmax
        xsgn = 1
    if y < 0:
        ylim = -ymin  # positive
        ysgn = -1
    else:
        ylim = ymax
        ysgn = 1
    # now extrude
    x2 = abs(ylim * x / y)  # try extruding till we reach the top
    if x2 <= xlim:  # have we dropped off the right?
        y2 = ylim  # if not, go with this
    else:  # but if we would have dropped off the right
        x2 = xlim  # just go as far right as possible instead
        y2 = abs(xlim * y / x)
    return x2 * xsgn, y2 * ysgn


def mpl_imshow_png(ax: Axes, png: bytes, **imshow_kws) -> AxesImage:
    """Helper function to display an image from a bytestream input, e.g. an encoded png in axes.

    Removes axes labels from `ax`.

    Args:
        ax (Axes): The axes to plot the encoded image to.
        png (bytes): The bytestream of the encoded image.

    Returns:
        AxesImage: The image plotted to the axes `ax`.
    """
    buffer = io.BytesIO()
    buffer.write(png)
    buffer.seek(0)
    img_array = np.array(PIL.Image.open(buffer))
    ax.axis('off')
    return ax.imshow(img_array, rasterized=True, **imshow_kws)
