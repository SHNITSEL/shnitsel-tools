from math import ceil
from typing import Any, Callable, Iterable, TYPE_CHECKING

from matplotlib.colors import Normalize
from matplotlib.figure import Figure, SubFigure
import numpy as np
from numpy.typing import NDArray
import xarray as xr
import matplotlib as mpl

from matplotlib.axes import Axes
from matplotlib.pyplot import subplot_mosaic

from scipy import stats
from sklearn.cluster import KMeans

from shnitsel.analyze.generic import get_standardized_pairwise_dists
from shnitsel.analyze.pca import pca

from .common import figax, extrude, mpl_imshow_png
from ...rd import highlight_pairs

if TYPE_CHECKING:
    from rdkit.Chem import Mol


def plot_noodleplot(
    noodle: NDArray | xr.DataArray,
    hops=None,
    fig: Figure | SubFigure | None = None,
    ax: Axes | None = None,
    c: NDArray | xr.DataArray = None,
    colorbar_label: str | None = None,
    cmap: str | None = None,
    cnorm: str | Normalize | None = None,
    noodle_kws: dict[str, Any] = None,
    hops_kws: dict[str, Any] = None,
) -> Axes:
    """Plot a matplotlib scatter-plot with specialized defaults

    Parameters
    ----------
    noodle
        Data to plot as a quanititatively colour-coded scatter-plot
    hops
        The coordinates of projections representing hopping points, by default None
    fig
        The :py:class:`matplotlib.pyplot.figure.Figure` object onto which to plot
        (If not provided, one will be created.)
    ax
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        (If not provided, one will be created.)
    c
        Data according to which to colour the main scatter-plot
    colorbar_label
        Label to describe the colorbar
    cmap, optional
        A :py:class:`matplotlib.colors.Colormap` used to colour the main scatter-plot
    cnorm, optional
        A :py:class:`matplotlib.colors.Normalize` object applied to
        the values in ``c`` before passing to ``cmap``.
        If not provided, linear normalization is used.
    noodle_kws, optional
        Keyword arguments for the main scatter-plot
    hops_kws, optional
        Keyword arguments for the hopping-point scatter-plot

    Returns
    -------
        The :py:class:matplotlib.axes.Axes` instance used
    """
    fig, ax = figax(fig=fig, ax=ax)
    if c is None:
        c = noodle['time']
        c_is_time = True
    else:
        c_is_time = False

    if colorbar_label is not None:
        pass
    elif hasattr(c, 'attrs') and 'long_name' in c.attrs:
        colorbar_label = c.attrs['long_name']
        if 'units' in c.attrs:
            colorbar_label = f"{colorbar_label} / {c.attrs['units']}"
    elif hasattr(c, 'name'):
        colorbar_label = c.name
        if 'units' in c.attrs:
            colorbar_label = f"{colorbar_label} / {c.attrs['units']}"
    elif c_is_time:
        colorbar_label = '$t$ / fs'

    cmap = cmap or mpl.colormaps['cividis_r']
    if isinstance(cmap, str):
        cmap = mpl.colormaps[cmap]
    cnorm = cnorm or mpl.colors.Normalize(c.min(), c.max())  # type: ignore

    # TODO: remove groupby? Needed only for line-plot or for legend
    # for trajid, traj in noodle.groupby('trajid'):
    #     ctraj = c.sel(trajid=trajid)
    noodle_kws = noodle_kws or {}
    noodle_kws = {'alpha': 0.5, 's': 0.2, **noodle_kws}
    sc = ax.scatter(
        noodle.isel(PC=0), noodle.isel(PC=1), c=c, cmap=cmap, norm=cnorm, **noodle_kws
    )

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    if hops is not None:
        hops_kws = dict(s=0.5, c='limegreen') | (hops_kws or {})
        ax.scatter(hops.isel(PC=0), hops.isel(PC=1), **hops_kws)

    fig.colorbar(sc, ax=ax, label=colorbar_label, pad=0.02)

    # Alternative layout solution
    # d = make_axes_locatable(ax)
    # cax = d.append_axes("right", size="5%", pad="2%")
    # fig.colorbar(pc, cax=cax, label='dihedral')

    assert isinstance(ax, Axes)
    return ax


# TODO: implement plotting of noodleplot using multi-coloured lines


def get_loadings(frames: xr.Dataset, mean=False) -> xr.DataArray:
    """Get the loadings for the PCA of pairwise distances
    for the positional data in ``frames``.

    Parameters
    ----------
    frames
        A Dataset with an 'atXYZ' data_var, which should have
        'atom' and 'direction' dimensions.

    Returns
    -------
        A DataArray of loadings with dimensions
        'PC' (principal component) and 'atomcomb' (atom combination,
        one for each pair of atoms).
    """
    atXYZ = frames['atXYZ']
    descr = get_standardized_pairwise_dists(atXYZ, mean=mean)
    _, pca_obj = pca(descr, 'atomcomb', return_pca_object=True)

    return xr.DataArray(
        data=pca_obj[-1].components_,
        dims=['PC', 'atomcomb'],
        coords=dict(atomcomb=descr.atomcomb),
        attrs={'natoms': frames.sizes['atom']},
    )


def plot_loadings(ax: Axes, loadings: xr.DataArray):
    """Plot all loadings as arrows.

    Parameters
    ----------
    ax
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        the loadings.
    loadings
        A DataArray of PCA loadings including an 'atomcomb' dimension;
        as produced by :py:func:`shnitsel.vis.plot.pca_biplot.get_loadings`.
    """
    for _, pcs in loadings.groupby('atomcomb'):
        assert len(pcs) == 2
        pc1, pc2 = pcs.item(0), pcs.item(1)
        ax.arrow(0, 0, pc1, pc2)
        a1, a2 = int(pcs['atomcomb_from']), int(pcs['atomcomb_to'])
        ax.text(pc1, pc2, f"{a1},{a2}")


def cluster_general(decider: Callable[[int, int], bool], n: int) -> list[list[int]]:
    """Cluster indices iteratively according to a provided function.

    Parameters
    ----------
    decider
        A function to decide whether two points can potentially share
        a cluster.
    n
        The number of indices to cluster.

    Returns
    -------
        A list of clusters, where each cluster is represented as a
        list of indices.
    """
    clustered = np.full((n,), False)
    clusters = []
    # for each item, if it has not been clustered,
    # put those later items which have not yet been clustered
    # in a cluster with it
    for i in range(n):
        if clustered[i]:
            continue

        cluster = [i]
        for j in range(i + 1, n):
            if clustered[j]:
                continue
            if decider(i, j):
                cluster.append(j)
                clustered[j] = True

        clusters.append(cluster)

    return clusters


def cluster_loadings(loadings: xr.DataArray, cutoff: float = 0.05) -> list[list[int]]:
    """Cluster loadings iteratively based on proximity on the
    principal component manifold

    Parameters
    ----------
    loadings
        A DataArray of loadings
    cutoff, optional
        An upper bound on the possible distances between a point
        in a cluster and other points, within which they will still
        be assigned to the smae cluster, by default 0.05

    Returns
    -------
        A list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``.
    """

    def dist(i, j, l):
        pc1, pc2 = l.isel(atomcomb=j).values - l.isel(atomcomb=i).values
        return (pc1**2 + pc2**2) ** 0.5

    def decider(i, j):
        nonlocal loadings, cutoff, dist
        return dist(i, j, loadings) <= cutoff

    n = loadings.sizes['atomcomb']
    return cluster_general(decider, n)


def plot_clusters(
    loadings: xr.DataArray,
    clusters: list[list[int]],
    ax: Axes | None = None,
    labels: list[str] | None = None,
):
    """Plot clusters of PCA loadings

    Parameters
    ----------
    loadings
        A DataArray of PCA loadings including an 'atomcomb' dimension;
        as produced by :py:func:`shnitsel.vis.plot.pca_biplot.get_loadings`.
    clusters
        A list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``; as produced
        by :py:func:`shnitsel.vis.plot.pca_biplot.get_clusters`.
    ax
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        (If not provided, one will be created.)
    labels
        Labels for the loadings; if not provided, loadings will be labelled
        according to indices of the atoms to which they relate.
    """
    fig, ax = figax(ax=ax)
    for i, cluster in enumerate(clusters):
        acs = loadings.isel(atomcomb=cluster)
        x, y = acs.mean(dim='atomcomb')
        s = (
            labels[i]
            if labels is not None
            else ' '.join([f'({a1},{a2})' for a1, a2 in acs.atomcomb.values])
        )
        ax.arrow(0, 0, x, y)
        ax.text(x, y, s)


def _get_clusters_coords(loadings, atomcomb_clusters):
    return np.array(
        [
            loadings.isel(atomcomb=c).mean(dim='atomcomb').values
            for c in atomcomb_clusters
        ]
    )


def _separate_angles(points: NDArray, min_angle: float = 10) -> dict[int, float]:
    """Group points based on their polar angles, and work out scale factors
    by which to place labels along the ray from origin to point when annotating
    points, intending to avoid overlaps between labels.

    Parameters
    ----------
    points
        An array of shape (npoints, 2)
    min_angle, optional
        The minimal difference in argument (angle from positive x-axis, in degrees),
        of two points, below which they will be considered part of the
        same cluster; by default 10

    Returns
    -------
        A dictionary mapping from indices (corresponding to ``points``)
        to scalefactors used to extrude the label away from the loading.
    """
    angles = [float(np.degrees(np.arctan2(x, y))) for x, y in points]

    def decider(i, j):
        nonlocal angles
        return (
            abs(angles[i] - angles[j]) <= min_angle
        )  # degrees. Edge case: -179 and 179

    angle_clusters = cluster_general(decider, len(angles))
    scalefactors = {}

    def calc_dist(point):
        x, y = point
        return (x**2 + y**2) ** 0.5

    for angle_cluster in angle_clusters:
        if len(angle_cluster) < 2:
            continue
        dists = np.array(
            [(idx, calc_dist(points[idx])) for idx in angle_cluster],
            dtype=[('idx', int), ('dist', float)],
        )
        dists.sort(order='dist')
        factor: float = 1
        for idx, dist in dists[::-1]:
            scalefactors[idx] = factor  # less extrusion for the smaller radius
            factor *= 0.8
    return scalefactors


def _filter_cluster_coords(coords, n):
    radii = [(x**2 + y**2) ** 0.5 for x, y in coords]
    angles = [np.degrees(np.arctan2(x, y)) for x, y in coords]
    res = set(np.argsort(radii)[-(n - 2) :])
    avg = np.mean(angles)
    splay = [abs(avg - angle) for angle in angles]
    return res.union(np.argsort(splay)[-2:])


def plot_clusters_insets(
    ax, loadings, clusters, mol, min_angle=10, inset_scale=1, show_at_most=None
):
    """Plot selected clusters of the loadings of a pairwise distance PCA,
    and interpretations of those loadings, as highlighted molecular structures inset upon
    the loadings plot.

    Parameters
    ----------
    ax
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        the loadings
    loadings
        A DataArray of PCA loadings including an 'atomcomb' dimension;
        as produced by :py:func:`shnitsel.vis.plot.pca_biplot.get_loadings`.
    clusters
        A list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``; as produced
        by :py:func:`shnitsel.vis.plot.pca_biplot.get_clusters`.
    mol
        An RDKit ``Mol`` object to be used for structure display.
    min_angle
        Where multiple clusters of loadings lie in similar directions from
        the origin, they will be grouped together and only their member with the
        greatest radius will be annotated with a highlighted structure.
    inset_scale
        A factor by which to scale the size of the inset highlighted structures.
    show_at_most
        Maximal number of clusters to show; if the number of clusters is greater than
        this value, the clusters with smallest radius will be excluded so that only this
        many remain.
    """
    points = _get_clusters_coords(loadings, clusters)
    if show_at_most is not None:
        indices = _filter_cluster_coords(points, show_at_most)
    else:
        indices = range(len(clusters))
    scalefactors = _separate_angles(points, min_angle)

    for i, cluster in enumerate(clusters):
        acs = loadings.isel(atomcomb=cluster)
        x, y = acs.mean(dim='atomcomb')
        arrow_color = 'k' if i in indices else (0, 0, 0, 0.5)
        ax.arrow(
            0, 0, x, y, head_width=0.01, length_includes_head=True, color=arrow_color
        )

        scale = scalefactors.get(i, 1)

        x2, y2 = extrude(x, y, *ax.get_xlim(), *ax.get_ylim())
        x2 *= 0.8 * scale
        y2 *= 0.8 * scale

        if i not in indices:
            continue

        ax.plot([x, x2], [y, y2], '--', c='darkgray', lw=0.5)

        ymin, ymax = ax.get_ylim()
        inset_size = inset_scale * np.array([7, 10]) * (ymax - ymin) / 65
        iax = ax.inset_axes([x2, y2, *inset_size], transform=ax.transData)
        iax.set_anchor('SW')  # keep bottom-left corner of image at arrow tip!

        png = highlight_pairs(mol, acs.atomcomb.values)
        mpl_imshow_png(iax, png)


# Compatability with old notebooks:
plot_clusters2 = plot_clusters_insets


def _get_axs(clusters, labels):
    naxs = min(len(clusters), len(labels))
    ncols = ceil(naxs**0.5)
    nblanks = naxs % ncols
    flat = labels[:naxs] + [None] * nblanks
    mosaic = np.array(flat).reshape(-1, ncols)
    _, axs = subplot_mosaic(mosaic)
    return axs


def plot_clusters_grid(
    loadings: xr.DataArray,
    clusters: list[list[int]],
    ax: Axes | None = None,
    labels: list[str] | None = None,
    axs: dict[Axes] | None = None,
    mol: 'Mol | None' = None,
):
    """Plot selected clusters of the loadings of a pairwise distance PCA,
    and interpretations of those loadings:

        - On the left, a large plot of selected clusters of loadings indicated as arrows
        - On the right, a grid of structures corresponding to
        structures of loadings; the pairs involved in the cluster
        are represented by colour-coding the atoms of the structures.

    Parameters
    ----------
    loadings
        A DataArray of PCA loadings including an 'atomcomb' dimension;
        as produced by :py:func:`shnitsel.vis.plot.pca_biplot.get_loadings`.
    clusters
        A list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``; as produced
        by :py:func:`shnitsel.vis.plot.pca_biplot.get_clusters`.
    ax
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        the loadings
        (If not provided, one will be created.)
    labels
        Labels for the loadings; if not provided, loadings will be labelled
        according to indices of the atoms to which they relate.
    axs
        A dictionary mapping from plot labels to :py:class:`matplotlib.pyplot.axes.Axes`
        objects
        (If not provided, one will be created.)
    mol
        An RDKit ``Mol`` object to be used for structure display
    """
    fig, ax = figax(ax=ax)
    if labels is None:
        labels = list('abcdefghijklmnopqrstuvwxyz')

    if axs is None:
        axs = _get_axs(clusters, labels)

    for mol_ax in axs.values():
        mol_ax.axis('off')

    for i, cluster in enumerate(clusters):
        acs = loadings.isel(atomcomb=cluster)
        x, y = acs.mean(dim='atomcomb')
        s = labels[i]
        ax.arrow(0, 0, x, y, head_width=0.01, length_includes_head=True)

        x2, y2 = extrude(x, y, *ax.get_xlim(), *ax.get_ylim())

        ax.plot([x, x2], [y, y2], '--', c='k', lw=0.5)
        ax.text(x2, y2, s)

        if axs is not None and mol is not None:
            png = highlight_pairs(mol, acs.atomcomb.values)
            mpl_imshow_png(axs[s], png)
            axs[s].set_title(s)

# Compatability with old notebooks:
plot_clusters3 = plot_clusters_grid

def circbins(
    angles, nbins=4, center=0
) -> tuple[list[Iterable[int]], list[tuple[float, float]]]:
    """Bin angular data by clustering unit-circle projections

    Parameters
    ----------
    angles
        Angles in degrees
    nbins, optional
        Number of bins to return, by default 4
    center, optional
        No longer used, by default 0

    Returns
    -------
    bins
        Indices of angles belonging to each bin
    edges
        Tuple giving a pair of boundary angles for each bin;
        the order of the bins corresponds to the order used in ``bins``
    """

    def proj(x):
        "project angles in degrees onto unit circle"
        x = x * np.pi / 180
        return np.c_[np.cos(x), np.sin(x)]

    kmeans = KMeans(n_clusters=nbins)

    labels = kmeans.fit_predict(proj(angles))

    space = np.linspace(0, 360, num=10)
    sample = kmeans.predict(proj(space))
    mask = np.diff(sample) != 0
    mask = np.concat([mask, [sample[-1] == sample[0]]])
    edgeps = space[mask]
    edges = list(zip(edgeps, np.roll(edgeps, -1)))
    label_at_edge = sample[mask][:-1]
    bins = [np.flatnonzero(labels == label) for label in label_at_edge]
    return bins, edges


def plot_bin_edges(
    angles: NDArray,
    radii: NDArray,
    bins: list[Iterable[int]],
    edges: list[tuple[float, float]],
    picks: list[int],
    ax: Axes,
    labels: list[str],
):
    """Illustrate how angles have been binned.

    Parameters
    ----------
    angles
        A 1D array of angles in degrees.
    radii
        A 1D array of radii, with order corresponding
        to ``angles``.
    bins
        Lists of bins, each bin represented as a list
        of indices.
    edges
        A pair of edges (angles in degrees) for each bin
        in ``bins``.
    picks
        A list of indices indicating which cluster has
        been chosen from each bin.
    ax
        An matplotlib ``Axes`` object onto which to plot;
        this should be set up with polar projection.
    labels
        One label for each entry in ``picks``.
    """
    rangles = np.radians(angles)

    for e in np.radians(edges):
        ax.plot([e, e], [0, 0.4], c='gray', ls='--', lw='1')

    for a, r, s in zip(rangles[picks], radii[picks], labels[: len(picks)]):
        ax.text(a, r, s, ha='left', va='bottom', fontsize=6)

    for b, c in zip(bins, list('rgbm')):
        # ax.plot(x, y)
        # colors = ['r' if x else 'b' for x in mask]
        ax.scatter(rangles[b], radii[b], c='gray', s=5)

    ax.scatter(rangles[picks], radii[picks], c='k', s=5)

    ax.set_rlabel_position(200)


def pick_clusters(frames, nbins, mean=False):
    """Calculate pairwise-distance PCA, cluster the loadings
    and pick a representative subset of the clusters.

    Parameters
    ----------
    frames
        An :py:class:`xarray.Dataset` with an 'atXYZ' variable
        having an 'atom' dimension
    nbins
        The number of bins to use when binning clusters of
        loadings according to the angle they make to the x-axis
        on the projection manifold

    Returns
    -------
        A dictionary with the following key-value pairs:

            - loadings: the loadings of the PCA
            - clusters: a list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``; as produced
        by :py:func:`shnitsel.vis.plot.pca_biplot.get_clusters`.
            - picks: the cluster chosen from each bin of clusters
            - angles: the angular argument (rotation from the positive x-axis) of each
            cluster center
            - center: the circular mean of the angle of all picked clusters
            - radii: The distance of each cluster from the origin
            - bins: Indices of angles belonging to each bin
            - edges: Tuple giving a pair of boundary angles for each bin;
        the order of the bins corresponds to the order used in ``bins``
    """
    loadings = get_loadings(frames, mean)
    clusters = cluster_loadings(loadings)
    points = _get_clusters_coords(loadings, clusters)

    angles = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    radii = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
    center = stats.circmean(angles, high=180, low=-180)

    picks, bins, edges = _binning_with_min_entries(
        nbins=nbins, angles=angles, radii=radii, return_bins_edges=True
    )
    # bins, edges = circbins(angles, nbins=4, center=center)
    # picks = [b[np.argmax(radii[b])] for b in bins]

    return dict(
        loadings=loadings,
        clusters=clusters,
        picks=picks,
        angles=angles,
        center=center,
        radii=radii,
        bins=bins,
        edges=edges,
    )


def _binning_with_min_entries(
    nbins: int,
    angles: NDArray,
    radii: NDArray,
    min_entries: int = 4,
    max_attempts: int = 10,
    return_bins_edges: bool = False,
):
    attempts = 0
    bins, edges = circbins(angles=angles, nbins=nbins)

    # Repeat binning until all bins have at least 'min_entries' or exceed max_attempts
    while any(arr.size == 0 for arr in bins) and attempts < max_attempts:
        print(
            f"Less than {min_entries} directions found, procedure repeated with another binning."
        )
        nbins += 1  # Increase the number of bins
        bins, edges = circbins(angles, nbins)
        attempts += 1

    # If max attempts were reached without satisfying condition
    if attempts >= max_attempts:
        print(f"Max attempts ({max_attempts}) reached. Returning current bins.")

    picks = [b[np.argmax(radii[b])] for b in bins]

    if return_bins_edges:
        return picks, bins, edges
    else:
        return picks
