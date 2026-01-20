import io
from xml.etree.ElementTree import Element

import PIL
import matplotlib as mpl
from matplotlib.image import AxesImage
from matplotlib.text import Annotation, Text
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
import logging
from svgpathtools import Document, Path, Line, CubicBezier, QuadraticBezier, Arc


def outlabel(ax: Axes, label: str) -> Text:
    """Adds a label just outside the top left corner of a plot (outside the axes).

    Parameters
    ----------
    ax : Axes
        The ``Axes`` object to annotate
    label : str
        The label to be added

    Returns
    -------
        The :py:class:`matplotlib.text.Text` instance created
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

    Parameters
    ----------
    ax : Axes
        The ``Axes`` object to annotate
    label : str
        The label to be added

    Returns
    -------
        The :py:class:`matplotlib.text.Annotation` instance created representing the inserted label.
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

    Parameters
    ----------
    fig : Figure | SubFigure | None, optional
        The optional figure to use as a basis for `ax` if the latter is not provided. Defaults to None.
    ax : Axes | None, optional
        The axes object provided. Will be used to populate `fig` if provided.. Defaults to None.

    Returns
    -------
    tuple[Figure | SubFigure, Axes]
        A complete combination of figure and axes.
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

    Parameters
    ----------
    x, y : float
        Coordinates of the vector to extrapolate
    xmin, xmax, ymin, ymax : float
        Bounds of the rectangle to the edge of which the
        ray should be extended

    Returns
    -------
    tuple[float, float]
        The position at the end of the extrusion, where the origin-ray through (x,y) intersects the boundary of the axes.
    """
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

    Parameters
    ----------
    ax : Axes
        The ``Axes`` object into which to plot
    png : bytes
        The bytestream data of the image to plot

    Returns
    -------
        ``AxesImage``, as returned by ``ax.imshow``
    """
    buffer = io.BytesIO()
    buffer.write(png)
    buffer.seek(0)
    img_array = np.array(PIL.Image.open(buffer))
    ax.axis('off')
    return ax.imshow(img_array, rasterized=True, **imshow_kws)


def mpl_svg_into_axes(ax: Axes, svg_string: str, chord_length: float = 1e-2) -> Axes:
    """Helper function to plot an SVG image represented by its str representation into
    provided `Axes`.

    Used to plot SVG graphics into a set of axes instead of a pixelated PNG.

    Parameters
    ----------
    ax : Axes
        The ``Axes`` object into which to plot
    svg_string : str
        The bytestream data of the image to plot
    chord_length,optional : float
        Length of cords of bezier curves to be drawn. Defaults to 1e-2.

    Returns
    -------
        ``Axes``, after the plotting of
    """
    doc = Document.from_svg_string(svg_string=svg_string)
    doc_paths = doc.paths()
    zorder_pos = np.linspace(1.0, 2.0, num=len(doc_paths))
    path: Path
    for path, zorder_pos_path in zip(doc_paths, zorder_pos):
        # print(path.__dict__)

        is_filled = False
        fill_color = None
        fill_opacity = 1.0

        is_stroked = False
        stroke_color = None
        stroke_width = None
        stroke_opacity = 1.0
        stroke_linecap = None
        join_style = 'miter'

        tree_elem: Element = path._tree_element

        if tree_elem is None:
            tree_elem: Element = path.element  # type: ignore

        if tree_elem is not None:
            if tree_elem.attrib is not None:
                if 'fill' in tree_elem.attrib:
                    fill_color = tree_elem.attrib['fill']
                    if fill_color == 'none':
                        is_filled = False
                    else:
                        is_filled = True

                if 'stroke' in tree_elem.attrib:
                    stroke_color = tree_elem.attrib['stroke']
                    if stroke_color == 'none':
                        is_stroked = False
                    else:
                        is_stroked = True

                if 'stroke-width' in tree_elem.attrib:
                    is_stroked = True
                    stroke_width = tree_elem.attrib['stroke-width']

                if 'stroke-linecap' in tree_elem.attrib:
                    stroke_linecap = tree_elem.attrib['stroke-linecap']

                if 'opacity' in tree_elem.attrib:
                    fill_opacity = float(tree_elem.attrib['opacity'])
                    stroke_opacity = fill_opacity

                if 'style' in tree_elem.attrib:
                    # Read CSS style info:
                    # print(tree_elem.attrib['style'])
                    single_css_attribs = list(tree_elem.attrib['style'].split(';'))
                    css_attrib_pairs = [
                        res
                        for part in single_css_attribs
                        if len(res := part.split(':')) == 2
                    ]
                    css_attribs = {k.strip(): v.strip() for (k, v) in css_attrib_pairs}
                    # print(css_attribs)

                    if 'fill' in css_attribs:
                        fill_color = css_attribs['fill']
                        if fill_color == 'none':
                            is_filled = False
                        else:
                            is_filled = True
                    if 'opacity' in css_attribs:
                        fill_opacity = float(css_attribs['opacity'])
                        stroke_opacity = fill_opacity
                    if 'stroke' in css_attribs:
                        stroke_color = css_attribs['stroke']
                        if stroke_color == 'none':
                            is_stroked = False
                        else:
                            is_stroked = True

                    if 'stroke-linecap' in css_attribs:
                        stroke_linecap = css_attribs['stroke-linecap']

                    if 'stroke-width' in css_attribs:
                        stroke_width = css_attribs['stroke-width']
                    if 'stroke-opacity' in css_attribs:
                        stroke_opacity = float(css_attribs['stroke-opacity'])
                    if 'stroke-linejoin' in css_attribs:
                        join_style = css_attribs['stroke-linejoin']

        full_path = []
        for segment in path:
            # TODO: COLOR!
            # color = path.
            pts = None
            if isinstance(segment, Line):
                start = segment.start
                end = segment.end
                pts = [start, end]
            elif isinstance(segment, QuadraticBezier):
                segment_length: float = segment.length()  # type: ignore
                num_lines = int(np.ceil(segment_length / chord_length))
                pts = list(segment.points(np.linspace(0, 1, num_lines + 1)))
                pass
            elif isinstance(segment, CubicBezier):
                segment_length: float = segment.length()  # type: ignore
                num_lines = int(np.ceil(segment_length / chord_length))
                pts = list(segment.points(np.linspace(0, 1, num_lines + 1)))
                pass
            elif isinstance(segment, Arc):
                arc_length: float = segment.length()  # type: ignore
                num_lines = int(np.ceil(arc_length / chord_length))
                pts = [segment.point(t) for t in np.linspace(0, 1, num_lines + 1)]
            else:
                logging.info(f"Ignoring part of svg path of type {type(segment)}.")
            if pts is not None:
                full_path += pts

        if len(full_path) > 0:
            full_path = np.array(full_path)

            if is_filled:
                if fill_opacity is None:
                    fill_opacity = 1.0

                poly = ax.fill(
                    np.real(full_path),
                    np.imag(full_path),
                    color=fill_color or 'k',
                    # alpha=fill_opacity,
                    zorder=zorder_pos_path,
                    capstyle=stroke_linecap if stroke_linecap is not None else 'butt',
                    edgecolor=(stroke_color or 'k') if is_stroked else None,
                )
                # print(poly, fill_color, fill_opacity)
            elif is_stroked:
                if stroke_opacity is None:
                    # logging.warning(f"lining {path}")
                    stroke_opacity = 1.0

                line = ax.plot(
                    np.real(full_path),
                    np.imag(full_path),
                    color=stroke_color or 'k',
                    # alpha=stroke_opacity,
                    solid_joinstyle=join_style or 'miter',
                    # capstyle=stroke_linecap
                    # if stroke_linecap is not None
                    # else 'butt',
                    lw=float(stroke_width.strip("ptxemc "))
                    if stroke_width is not None
                    else None,
                )

    # ax.axis('off')
    ax.invert_yaxis()
    return ax
