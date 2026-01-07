from typing import Literal

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from shnitsel.analyze.pca import pca_and_hops
from shnitsel.geo.geocalc import distance, angle, dihedral
from . import pca_biplot as pb
from .common import figax
from shnitsel.bridges import default_mol, set_atom_props


def fit_kdes(noodle, geo_prop, geo_filter):
    """Fit multiple subsets of data using kernel density esimation

    Parameters
    ----------
    noodle
        Data to fit
    geo_prop
        Data by which to filter ``noodle`` into different subsets to be
        fitted independently
    geo_filter
        A list of ``(minimum, maximum)`` pairs, to which ``geo_prop`` will
        be compared; one KDE will be fitted for the data falling in the
        range defined by each pair.

    Returns
    -------
        A list of kernels; one for each pair in ``geo_filter``.

    Raises
    ------
    ValueError
        If any of the ``geo_filter`` ranges is such that no points from
        ``geo_prop`` fall within it
    """
    kernels = []
    for p1, p2 in geo_filter:
        mask = (p1 < geo_prop) & (geo_prop < p2)
        subset = noodle.sel(frame=mask).T
        if subset.size == 0:
            raise ValueError(f"No points in range {p1} < x < {p2}")
        kernels.append(stats.gaussian_kde(subset))
    return kernels


def eval_kdes(kernels: list, xx, yy):
    """Evaluate each of a list of kernels over the same meshgrid

    Parameters
    ----------
    kernels
        A list of kernel density estimate objects (from scikit-learn)
    xx
        The meshgrid values for x
    yy
        The meshgrid values for y (shape must much ``xx``)

    Returns
    -------
        A list of arrays of the same shape as ``xx`` and ``yy``;
        each representing the height each fitted kernel assigns
        to the corresponding x, y coordinate from ``xx`` and ``yy``.
    """
    xys = np.c_[xx.ravel(), yy.ravel()].T
    Zs = []
    for k in kernels:
        Z = k.evaluate(xys)
        Z = Z.reshape(xx.shape) / Z.max()
        Zs.append(Z)
    return Zs


def get_xx_yy(noodle, fineness=500, extension=0.1):
    """Produce a meshgrid over the range of the data

    Parameters
    ----------
    noodle
        The data from which to derive a range
    fineness
        The fineness of the grid, by default 500
    extension
        How far to extend beyond the range,
        as a fraction of the mean-extremum difference, by default 0.1

    Returns
    -------
        Two numpy NDArrays, as returned by np.meshgrid
    """
    means = noodle.mean(dim='frame').values
    mins = noodle.min(dim='frame').values
    mins -= (means - mins) * extension
    maxs = noodle.max(dim='frame').values
    maxs += (maxs - means) * extension
    ls = np.linspace(mins, maxs, num=fineness).T
    xx, yy = np.meshgrid(ls[0], ls[1])
    return xx, yy


def fit_and_eval_kdes(noodle, geo_prop, geo_filter, fineness=500, extension=0.1):
    """Get kernel density estimates and a range-appropriate meshgrid for data

    Parameters
    ----------
    noodle
        Data to fit
    geo_prop
        Data by which to filter ``noodle`` into different subsets to be
        fitted independently
    geo_filter
        A  list of ``(minimum, maximum)`` pairs, to which ``geo_prop`` will
        be compared; one KDE will be fitted for the data falling in the
        range defined by each pair.
    fineness, optional
        The finess of the meshgrid, by default 500
    extension, optional
        How far to extend the meshgrid beyond the range of the data,
        as a fraction of the mean-extremum difference, by default 0.1

    Returns
    -------
    xx, yy
        Meshgrid values
    Zs
        The KDE-interpolated height corresponding to x and y values in
        xx and yy respectively
    """
    noodle = noodle.transpose('frame', 'PC')  # required order for the following 3 lines

    xx, yy = get_xx_yy(noodle, fineness=fineness, extension=extension)
    kernels = fit_kdes(noodle, geo_prop, geo_filter)
    return xx, yy, eval_kdes(kernels, xx, yy)


def plot_kdes(xx, yy, Zs, colors=None, levels=None, fill=True, fig=None, ax=None):
    """Plot contours of kernel density estimates

    Parameters
    ----------
    xx
        An array of x values
    yy
        An array of y values (must have the same shape as ``xx``)
    Zs
        A list of arrays of z values (each array must have the same
        shape as ``xx`` and ``yy``)
    colors
        A set of colours accepted by matplotlib (e.g. a colormap)
    levels
        Determines the number and positions of the contour lines / regions.
        (Passed to ``matplotlib.pyplot.contour``)
    fill
        Whether to fill in the outlined contours
        (i.e. whether to use ``matplotlib.pyplot.contour`` or
        ``matplotlib.pyplot.contourf``).
    fig
        A matplotlib ``Figure`` object into which to draw
        (if not provided, a new one will be created)
    ax
        A matplotlib ``Axes`` object into which to draw
        (if not provided, a new one will be created)
    """
    fig, ax = figax(fig=fig, ax=ax)
    if colors is None:
        if len(Zs) == 2:
            colors = ['purple', 'green']
        else:
            colors = plt.get_cmap('tab10').colors

    for Z, c in zip(Zs, colors):
        if fill:
            ax.contourf(xx, yy, Z, levels=levels, colors=c, alpha=0.1)
        ax.contour(xx, yy, Z, levels=levels, colors=c, linewidths=0.5)


def biplot_kde(
    frames,
    at1: int = 0,
    at2: int = 1,
    at3: int | None = None,
    at4: int | None = None,
    geo_filter: list[tuple[float, float]] | None = None,
    levels: int | list[float] | None = None,
    scatter_color: Literal['time', 'geo'] = 'time',
    fill: bool = True,
    nbins=4,
    mean=False
):
    """\
    Generates a biplot that visualizes PCA projections and kernel density estimates (KDE) 
    of a property (distance, angle, dihedral angle) describing the geometry of specified
    atoms. The property is chosen based on the number of atoms specified:
    
    * 2 atoms => distance
    * 3 atoms => angle
    * 4 atoms => dihedral angle

    Parameters
    ----------
    frames
        A dataset containing trajectory frames with atomic coordinates.
    at1, at2, at3, at4
        Indices of the first, second, third and fourth atoms for geometric property calculation.
    geo_filter
        A list of tuples representing ranges. A KDE is plotted for each range, indicating the distribution of
        points for which the value of the geometry feature falls in that range.
    levels
        Contour levels for the KDE plot. Defaults to [0.08, 1]. This parameter is passed to
        matplotlib.axes.Axes.contour.
    scatter_color
        Must be one of 'time' or 'geo'. If 'time', the scatter-points will be colored based on the time coordinate;
        if 'geo', the scatter-points will be colored based on the relevant geometry feature (see above).
    fill
        Whether to plot filled contours (``fill=True``, uses ``ax.contourf``)
        or just contour lines (``fill=False``, uses ``ax.contour``).
    nbins
        number of bins to be visualized, must be an integer between 1 and 4

    Returns
    -------
    kde_dat
        The computed KDE data for the atom-atom distance distribution.

    Notes
    -----
    * Computes a geometric property of the specified atoms across all frames.
    * Uses kernel density estimation (KDE) to analyze the distance distributions.
    * Performs PCA on trajectory pairwise distances and visualizes clustering of structural changes.
    * Produces a figure with PCA projection, cluster analysis, and KDE plots.
    """
    if scatter_color not in {'time', 'geo'}:
        raise ValueError("`scatter_color` must be 'time' or 'geo'")

    if levels is None:
        levels = [0.08, 1]

    match at1, at2, at3, at4:
        case at1, at2, None, None:
            # compute distance between atoms at1 and at2
            geo_prop = distance(frames['atXYZ'], at1, at2)
            if not geo_filter:
                geo_filter = [(0, 3), (5, 100)]
        case at1, at2, at3, None:
            # compute angle between vectors at1 - at2 and at2 - at3
            assert at3 is not None  # to satisfy the typechecker
            geo_prop = angle(frames['atXYZ'], at1, at2, at3, deg=True)
            if not geo_filter:
                geo_filter = [(0, 80), (110, 180)]
        case at1, at2, at3, at4:
            # compute dihedral defined as angle between normals to planes (at1, at2, at3) and (at2, at3, at4)
            assert at3 is not None
            assert at4 is not None
            geo_prop = dihedral(frames['atXYZ'], at1, at2, at3, at4, deg=True)
            if not geo_filter:
                geo_filter = [(0, 80), (110, 180)]

    # prepare layout
    fig, oaxs = plt.subplots(1, 2, layout='constrained', width_ratios=[3, 2])
    fig.set_size_inches(8.27, 11.69 / 3)  # a third of a page, spanning both columns
    gs = oaxs[0].get_subplotspec().get_gridspec()
    for ax in oaxs:
        ax.remove()
    pcasf = fig.add_subfigure(gs[0])
    pcaax = pcasf.subplots(1, 1)
    structsf = fig.add_subfigure(gs[1])
    structaxs = structsf.subplot_mosaic('ab\ncd')

    # prepare data
    noodle, hops = pca_and_hops(frames, mean=mean)
    kde_data = fit_and_eval_kdes(noodle, geo_prop, geo_filter, fineness=100)
    d = pb.pick_clusters(frames, nbins=nbins, mean=mean)
    loadings, clusters, picks = d['loadings'], d['clusters'], d['picks']
    mol = default_mol(frames)
    mol = set_atom_props(mol, atomLabel=True, atomNote=[''] * mol.GetNumAtoms())

    if scatter_color == 'time':
        noodleplot_c = None
        noodleplot_cmap = None
    elif scatter_color == 'geo':
        noodleplot_c = geo_prop
        noodleplot_cmap = 'PRGn'
    else:
        assert False

    pb.plot_noodleplot(
        noodle,
        hops,
        c=noodleplot_c,
        cmap=noodleplot_cmap,
        ax=pcaax,
        noodle_kws=dict(alpha=1, marker='.'),
        hops_kws=dict(c='r', s=0.2),
    )

    # in case more clusters were found than we have room for:
    picks = picks[:4]

    pb.plot_clusters_grid(
        loadings,
        [clusters[i] for i in picks],
        ax=pcaax,
        axs=structaxs,
        mol=mol,
        labels=list('abcd'),
    )
    xx, yy, Zs = kde_data
    plot_kdes(xx, yy, Zs, levels=levels, fill=fill, ax=pcaax)

    return kde_data


def plot_cdf_for_kde(z, level, ax=None):
    """Plot the cumulative density for a KDE, to show what
    proportion of points are contained by contours at a
    given density ``level``

    Parameters
    ----------
    z
        The values from the kernel evaluated over the input
        space
    level
        The cumulative density corresponding to this level
        will be marked on the graph
    ax
        A :py:class:`matplotlib.axes.Axes` object into which
        to plot. (If not provided, one will be created.)

    Returns
    -------
    y
        The proportion of points contained by contours placed
        at density ``level``
    """
    fig, ax = figax(ax=ax)
    bins, edges, _ = ax.hist(
        z,
        bins=1000,
        range=(0, 1.1 * level),
        cumulative=True,
        density=True,
        histtype='step',
    )
    y = bins[abs(edges - level).argmin()]
    ax.plot([0, level], [y, y], c='r')
    ax.plot([level, level], [0, y], c='r')
    return y
