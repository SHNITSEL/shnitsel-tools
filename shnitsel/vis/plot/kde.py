from typing import Iterable, Literal, Sequence

from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.figure import Figure, SubFigure
import numpy as np
from numpy.typing import ArrayLike
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

from shnitsel.analyze.pca import PCAResult, pca_and_hops
from shnitsel.data.dataset_containers import Frames, Trajectory, wrap_dataset
from shnitsel.filtering.structure_selection import StructureSelection
from shnitsel.geo.geocalc_.angles import angle
from shnitsel.geo.geocalc_.dihedrals import dihedral
from shnitsel.geo.geocalc_.distances import distance

# from shnitsel.geo.geocalc import distance, angle, dihedral
from . import pca_biplot as pb
from .common import figax
from shnitsel.bridges import construct_default_mol, set_atom_props
import xarray as xr


def _fit_kdes(
    pca_data: xr.DataArray,
    geo_property: xr.DataArray,
    geo_kde_ranges: Sequence[tuple[float, float]],
) -> Sequence[stats.gaussian_kde]:
    """\
    Fit a set of KDEs to the `pca_data`, after it has been split into subsets based on the values of
    `geo_property`. 
    
    The parameter `geo_kde_ranges` specifies the subsets of the values of `geo_property`
    that should be filtered into the same subset. 
    Returns one KDE for each such subset.


    Parameters
    ----------
    pca_data : xr.DataArray
        The pca data for which KDEs should be fitted on the various ranges.
    geo_property : xr.DataArray
        The geometric property that the data should be clustered/filtered by.
    geo_kde_ranges : Sequence[tuple[float, float]]
        The sequence of (distinct) ranges of values of the geometric property
        that the `pca_data` should be divided by.

    Returns
    ----------
    Sequence[stats.gaussian_kde]
        The sequence of fitted KDEs (kernels) for each range of `geo_kde_ranges`.
    Raises
    ------
    ValueError
        If any of the ``geo_filter`` ranges is such that no points from
        ``geo_prop`` fall within it
    """
    kernels = []
    for p1, p2 in geo_kde_ranges:
        mask = (p1 < geo_property) & (geo_property < p2)
        subset = pca_data.sel(frame=mask).T  # Swap leading frame dimension to the end?
        if subset.size == 0:
            raise ValueError(f"No points in range {p1} < x < {p2}")
        kernels.append(stats.gaussian_kde(subset))
    return kernels


def _eval_kdes(
    kernels: Sequence[stats.gaussian_kde], xx: np.ndarray, yy: np.ndarray
) -> Sequence[np.ndarray]:
    """Evaluate all fitted gaussian kernel density estimators on a mesh-grid
    and return the results.


    Parameters
    ----------
        kernels : Sequence[stats.gaussian_kde]
            The transformed pca data to get the supporting mesh grid for.
        xx : np.ndarray
            The x coordinates of the mesh grid.
        yy : np.ndarray
            The y coordinates of the mesh grid.

    Returns
    ----------
        Sequence[np.ndarray]
            The sequence of evaluated approximate probability densities
            at the positions described by `xx` and `yy` for each and every
            individual KDE provided in `kernels`.
    """
    xys = np.c_[xx.ravel(), yy.ravel()].T
    Zs = []
    for k in kernels:
        Z = k.evaluate(xys)
        Z = Z.reshape(xx.shape) / Z.max()
        Zs.append(Z)
    return Zs


def _get_xx_yy(
    pca_data: xr.DataArray, num_steps: int = 500, extension: float = 0.1
) -> tuple[np.ndarray, np.ndarray]:
    """Get appropriately over-sized mesh-grids for x and y coordinates
    with an excess overhang of `extension` relative to the min/max-to-mean distance
    and `num_steps` intermediate steps between the upper and lower bound.

    Statistical properties will be derived from `pca_data`.


    Parameters
    ----------
    pca_data: xr.DataArray
        The transformed pca data to get the supporting mesh grid for.
    num_steps, optional : int
        Number of intermediate steps to generate in the grid. Defaults to 500.
    extension, optional : float
        Excess overhang beyond minima and maxima in x and y direction
        relative to their distance from the mean. Defaults to 0.1.

    Returns
    ----------
    tuple[np.ndarray, np.ndarray]
        First the numpy array holding x positions of a meshgrid
        Then the array holding y positions of a meshgrid.
    """
    means: np.ndarray = pca_data.mean(dim='frame').values
    mins: np.ndarray = pca_data.min(dim='frame').values
    mins -= (means - mins) * extension
    maxs: np.ndarray = pca_data.max(dim='frame').values
    maxs += (maxs - means) * extension
    ls = np.linspace(mins, maxs, num=num_steps).T
    xx, yy = np.meshgrid(ls[0], ls[1])
    return xx, yy


def _fit_and_eval_kdes(
    pca_data: PCAResult,
    geo_property: xr.DataArray,
    geo_kde_ranges: Sequence[tuple[float, float]],
    num_steps: int = 500,
    extension: float = 0.1,
) -> tuple[np.ndarray, np.ndarray, Sequence[np.ndarray]]:
    """Fit KDEs for each range of the `geo_kde_ranges` and filter by the value of `geo_property`
    being within the respective range.
    Then return a mesh grid and the evaluation of these kernel estimators on that mash grid.


    Parameters
    ----------
    pca_data: xr.DataArray
        The transformed pca data to get the supporting mesh grid for and extract
        the KDEs from.
    geo_property : xr.DataArray
        The geometric property that the data should be clustered/filtered by.
    geo_kde_ranges : Sequence[tuple[float, float]]
        The sequence of (distinct) ranges of values of the geometric property
        that the `pca_data` should be divided by.
    num_steps, optional : int
        Number of intermediate steps to generate in the grid. Defaults to 500.
    extension, optional : float
        Excess overhang beyond minima and maxima in x and y direction
        relative to their distance from the mean. Defaults to 0.1.

    Returns
    ----------
        tuple[np.ndarray, np.ndarray, Sequence[np.ndarray]]
        First the numpy array holding x positions of a meshgrid.
        Then the array holding y positions of a meshgrid.
        Last the Sequence of KDE evaluations on the meshgrid for each filter range.
    """
    # TODO: FIXME: We should be able to deal with a missing `frame` dimension.
    pca_data_da = pca_data.projected_inputs.transpose(
        'frame', 'PC'
    )  # required order for the following 3 lines

    xx, yy = _get_xx_yy(pca_data_da, num_steps=num_steps, extension=extension)
    kernels = _fit_kdes(pca_data_da, geo_property, geo_kde_ranges)
    return xx, yy, _eval_kdes(kernels, xx, yy)


def _plot_kdes(
    xx: np.ndarray,
    yy: np.ndarray,
    Zs: Sequence[np.ndarray],
    colors: Iterable | None = None,
    contour_levels: int | list[float] | None = None,
    contour_fill: bool = True,
    fig: Figure | SubFigure | None = None,
    ax: Axes | None = None,
):
    """Plot contours of kernel density estimates

    Parameters
    ----------
    xx : np.ndarray
        An array of x values
    yy : np.ndarray
        An array of y values (must have the same shape as ``xx``)
    Zs : Sequence[np.ndarray]
        A list of arrays of z values (each array must have the same
        shape as ``xx`` and ``yy``)
    colors : Iterable, optional
        A set of colours accepted by matplotlib (e.g. a colormap) of at least the same length as Zs
    contour_levels : int | list[float], optional
        Determines the number and positions of the contour lines / regions.
        (Passed to ``matplotlib.pyplot.contour``)
    contour_fill : bool, optional
        Whether to fill in the outlined contours
        (i.e. whether to use ``matplotlib.pyplot.contour`` or
        ``matplotlib.pyplot.contourf``).
    fig : Figure | SubFigure, optional
        A matplotlib ``Figure`` object into which to draw
        (if not provided, a new one will be created)
    ax : Axes, optional
        A matplotlib ``Axes`` object into which to draw
        (if not provided, a new one will be created)
    """
    fig, ax = figax(fig=fig, ax=ax)
    if colors is None:
        if len(Zs) == 2:
            colors = ['purple', 'green']
        else:
            colors = plt.get_cmap('inferno')(range(len(Zs)))

    for Z, c in zip(Zs, colors):
        if contour_fill:
            ax.contourf(xx, yy, Z, levels=contour_levels, colors=c, alpha=0.1)
        ax.contour(xx, yy, Z, levels=contour_levels, colors=c, linewidths=0.5)


def biplot_kde(
    frames: xr.Dataset | Frames | Trajectory,
    at1: int = 0,
    at2: int = 1,
    at3: int | None = None,
    at4: int | None = None,
    feature_selection: StructureSelection | None = None,
    geo_kde_ranges: Sequence[tuple[float, float]] | None = None,
    scatter_color_property: Literal['time', 'geo'] = 'time',
    geo_cmap: str | None = 'PRGn',  # any valid cmap type
    time_cmap: str | None = 'cividis',  # any valid cmap type
    contour_levels: int | list[float] | None = None,
    contour_colors: list[str] | None = None,
    contour_fill: bool = True,
    num_bins: Literal[1, 2, 3, 4] = 4,
    fig: mpl.figure.Figure | None = None,
    center_mean: bool = False,
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
    at1, at2, at3, at4: int
        Indices of the first, second, third and fourth atoms for geometric property calculation.
    feature_selection: StructureSelection, optional
        An optional selection of features/structure to use for the PCA analysis.
    geo_kde_ranges
        A Sequence of tuples representing ranges. A KDE is plotted for each range, indicating the distribution of
        points for which the value of the geometry feature falls in that range.
        Default values are chosen depending on the type of feature that should be analyzed. 
    contour_levels
        Contour levels for the KDE plot. Either the number of contour levels as an int or the list of floating 
        point values at which the contour lines should be drawn. Defaults to [0.08, 1]. 
        This parameter is passed to matplotlib.axes.Axes.contour.
    scatter_color_property
        Must be one of 'time' or 'geo'. If 'time', the scatter-points will be colored based on the time coordinate;
        if 'geo', the scatter-points will be colored based on the relevant geometry feature (see above).
    geo_cmap
        The Colormap to use for the noodleplot, if ``scatter_color='geo'``; this also determines contour
        colors unless ``contour_colors`` is set.
    time_cmap
        The Colormap to use for the noodleplot, if ``scatter_color='time'``.
    contour_fill
        Whether to plot filled contours (``contour_fill=True``, uses ``ax.contourf``)
        or just contour lines (``contour_fill=False``, uses ``ax.contour``).
    contour_colors
        An iterable (not a Colormap) of colours (in a format matplotlib will accept) to use for the contours.
        By default, the ``geo_cmap`` will be used; this defaults to 'PRGn'.
    num_bins
        number of bins to be visualized, must be an integer between 1 and 4
    fig
        matplotlib.figure.Figure object into which the plot will be drawn;
        if not provided, one will be created using ``plt.figure(layout='constrained')``
    center_mean
        Flag whether PCA data should be mean-centered before analysis. Defaults to False.

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
    assert at1 is not None and at2 is not None, (
        "Indices of first two atoms to `biplot_kde()` must not be None."
    )

    if scatter_color_property not in {'time', 'geo'}:
        raise ValueError("`scatter_color` must be 'time' or 'geo'")

    if contour_levels is None:
        contour_levels = [0.08, 1]

    wrapped_ds = wrap_dataset(frames, Frames | Trajectory)

    match at1, at2, at3, at4:
        case at1, at2, None, None:
            # compute distance between atoms at1 and at2
            geo_prop = distance(wrapped_ds.positions, at1, at2)
            if not geo_kde_ranges:
                geo_kde_ranges = [(0, 3), (5, 100)]
        case at1, at2, at3, None:
            # compute angle between vectors at1 - at2 and at2 - at3
            assert at3 is not None  # to satisfy the typechecker
            geo_prop = angle(wrapped_ds.positions, at1, at2, at3, deg=True)
            if not geo_kde_ranges:
                geo_kde_ranges = [(0, 80), (110, 180)]
        case at1, at2, at3, at4:
            # compute dihedral defined as angle between normals to planes (at1, at2, at3) and (at2, at3, at4)
            assert at3 is not None
            assert at4 is not None
            geo_prop = dihedral(wrapped_ds.positions, at1, at2, at3, at4, deg=True)
            if not geo_kde_ranges:
                geo_kde_ranges = [(0, 80), (110, 180)]

    # prepare layout
    if fig is None:
        fig = plt.figure(layout='constrained')

    oaxs = fig.subplots(1, 2, width_ratios=[3, 2])

    fig.set_size_inches(8.27, 11.69 / 3)  # a third of a page, spanning both columns
    gs = oaxs[0].get_subplotspec().get_gridspec()
    for ax in oaxs:
        ax.remove()
    pcasf = fig.add_subfigure(gs[0])
    pcaax = pcasf.subplots(1, 1)
    structsf = fig.add_subfigure(gs[1])
    structaxs = structsf.subplot_mosaic('ab\ncd')

    # prepare data
    pca_data, hops_mask = pca_and_hops(
        wrapped_ds, feature_selection=feature_selection, center_mean=center_mean
    )
    kde_data = _fit_and_eval_kdes(pca_data, geo_prop, geo_kde_ranges, num_steps=100)
    d = pb.pick_clusters(pca_data, num_bins=num_bins, center_mean=center_mean)
    loadings, clusters, picks = d['loadings'], d['clusters'], d['picks']
    mol = construct_default_mol(wrapped_ds)
    mol = set_atom_props(mol, atomLabel=True, atomNote=[''] * mol.GetNumAtoms())

    if scatter_color_property == 'time':
        noodleplot_c = None
        noodleplot_cmap = None
    elif scatter_color_property == 'geo':
        noodleplot_c = geo_prop
        noodleplot_cmap = 'PRGn'
    else:
        assert False

    # noodleplot_c, noodleplot_cmap = {
    #     'time': (None, time_cmap),
    #     'geo': (geo_prop, geo_cmap),
    # }[scatter_color_property]

    # noodle_cnorm = mpl.colors.Normalize(noodleplot_c.min(), noodle_c.max())
    # noodle_cscale = mpl.cm.ScalarMappable(norm=noodle_cnorm, cmap=noodle_cmap)

    pb.plot_noodleplot(
        pca_data.projected_inputs,
        hops_mask,
        c=noodleplot_c,
        cmap=noodleplot_cmap,
        # cnorm=noodle_cnorm,
        ax=pcaax,
        noodle_kws=dict(alpha=1, marker='.'),
        hops_kws=dict(c='r', s=0.2),
    )

    # in case more clusters were found than we have room for:
    picks = picks[:4]

    print(pca_data.explain_loadings())

    pb.plot_clusters_grid(
        loadings,
        [clusters[i] for i in picks],
        ax=pcaax,
        axs=structaxs,
        mol=mol,
        labels=list('abcd'),
    )

    if contour_colors is None:
        contour_colors = plt.get_cmap(noodleplot_cmap)(
            np.linspace(0, 1, len(contour_levels))
        )

    xx, yy, Zs = kde_data
    _plot_kdes(
        xx,
        yy,
        Zs,
        colors=contour_colors,
        contour_levels=contour_levels,
        contour_fill=contour_fill,
        ax=pcaax,
    )

    return kde_data


def plot_cdf_for_kde(
    z: np.ndarray, contour_level: float, ax: Axes | None = None
) -> float:
    """Plot the cumulative density for a KDE, to show what
    proportion of points are contained by contours at a
    given density ``level``

    Parameters
    ----------
    z : np.ndarray
        The values from the kernel evaluated over the input
        space
    contour_level : float
        The cumulative density corresponding to this level
        will be marked on the graph
    ax : Axes, optional
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
        range=(0, 1.1 * contour_level),
        cumulative=True,
        density=True,
        histtype='step',
    )
    y = float(bins[abs(edges - contour_level).argmin()])
    ax.plot([0, contour_level], [y, y], c='r')
    ax.plot([contour_level, contour_level], [0, y], c='r')
    return y
