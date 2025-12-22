from typing import Iterable, Literal, Sequence

from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.figure import Figure, SubFigure
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from shnitsel.analyze.pca import pca_and_hops
from shnitsel.core.typedefs import Frames
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
            The sequence of fitted KDEs for each range of `geo_kde_ranges`.
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
    pca_data: xr.DataArray,
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
    pca_data = pca_data.transpose(
        'frame', 'PC'
    )  # required order for the following 3 lines

    xx, yy = _get_xx_yy(pca_data, num_steps=num_steps, extension=extension)
    kernels = _fit_kdes(pca_data, geo_property, geo_kde_ranges)
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
    frames: Frames,
    at1: int = 0,
    at2: int = 1,
    at3: int | None = None,
    at4: int | None = None,
    geo_kde_ranges: Sequence[tuple[float, float]] | None = None,
    scatter_color_property: Literal['time', 'geo'] = 'time',
    contour_levels: int | list[float] | None = None,
    contour_fill: bool = True,
    num_bins: Literal[1, 2, 3, 4] = 4,
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
    contour_fill
        Whether to plot filled contours (``contour_fill=True``, uses ``ax.contourf``)
        or just contour lines (``contour_fill=False``, uses ``ax.contour``).
    num_bins
        number of bins to be visualized, must be an integer between 1 and 4
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
    assert (
        at1 is not None and at2 is not None
    ), "Indices of first two atoms to `biplot_kde()` must not be None."

    if scatter_color_property not in {'time', 'geo'}:
        raise ValueError("`scatter_color` must be 'time' or 'geo'")

    if contour_levels is None:
        contour_levels = [0.08, 1]

    match at1, at2, at3, at4:
        case at1, at2, None, None:
            # compute distance between atoms at1 and at2
            geo_prop = distance(frames['atXYZ'], at1, at2)
            if not geo_kde_ranges:
                geo_kde_ranges = [(0, 3), (5, 100)]
        case at1, at2, at3, None:
            # compute angle between vectors at1 - at2 and at2 - at3
            assert at3 is not None  # to satisfy the typechecker
            geo_prop = angle(frames['atXYZ'], at1, at2, at3, deg=True)
            if not geo_kde_ranges:
                geo_kde_ranges = [(0, 80), (110, 180)]
        case at1, at2, at3, at4:
            # compute dihedral defined as angle between normals to planes (at1, at2, at3) and (at2, at3, at4)
            assert at3 is not None
            assert at4 is not None
            geo_prop = dihedral(frames['atXYZ'], at1, at2, at3, at4, deg=True)
            if not geo_kde_ranges:
                geo_kde_ranges = [(0, 80), (110, 180)]

    # prepare layout
    # TODO: this should probably be standardized somewhere.
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
    pca_data, hops = pca_and_hops(frames, center_mean=center_mean)
    kde_data = _fit_and_eval_kdes(pca_data, geo_prop, geo_kde_ranges, num_steps=100)
    d = pb.pick_clusters(frames, num_bins=num_bins, center_mean=center_mean)
    loadings, clusters, picks = d['loadings'], d['clusters'], d['picks']
    mol = construct_default_mol(frames)
    mol = set_atom_props(mol, atomLabel=True, atomNote=[''] * mol.GetNumAtoms())

    if scatter_color_property == 'time':
        noodleplot_c = None
        noodleplot_cmap = None
    elif scatter_color_property == 'geo':
        noodleplot_c = geo_prop
        noodleplot_cmap = 'PRGn'
    else:
        assert False

    pb.plot_noodleplot(
        pca_data,
        hops,
        c=noodleplot_c,
        cmap=noodleplot_cmap,
        ax=pcaax,
        noodle_kws=dict(alpha=1, marker='.'),
        hops_kws=dict(c='r', s=0.2),
    )

    # in case more clusters were found than we have room for:
    picks = picks[:4]

    pb.plot_clusters3(
        loadings,
        [clusters[i] for i in picks],
        ax=pcaax,
        axs=structaxs,
        mol=mol,
        labels=list('abcd'),
    )
    xx, yy, Zs = kde_data
    _plot_kdes(
        xx, yy, Zs, contour_levels=contour_levels, contour_fill=contour_fill, ax=pcaax
    )

    return kde_data


def plot_cdf_for_kde(z: np.ndarray, contour_level: float, ax: Axes | None = None):
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
