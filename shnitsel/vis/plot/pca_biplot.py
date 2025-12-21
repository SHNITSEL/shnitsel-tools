from typing import Any, Callable, Iterable

from matplotlib.colors import Normalize
from matplotlib.figure import Figure, SubFigure
import numpy as np
import xarray as xr
import matplotlib as mpl

from matplotlib.axes import Axes

from scipy import stats
from sklearn.cluster import KMeans

from shnitsel.analyze.generic import get_standardized_pairwise_dists
from shnitsel.analyze.pca import pca

from .common import figax, extrude, mpl_imshow_png
from ...rd import highlight_pairs


def plot_noodleplot(
    noodle: np.NDArray | xr.DataArray,
    hops=None,
    fig: Figure | SubFigure | None = None,
    ax: Axes | None = None,
    c: np.NDArray | xr.DataArray = None,
    colorbar_label: str | None = None,
    cmap: str | None = None,
    cnorm: str | Normalize | None = None,
    cscale: mpl.cm.ScalarMappable = None,
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
    cscale, optional
        A :py:class:`matplotlib.cm.ScalarMappable`
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
    elif hasattr(c, 'name'):
        colorbar_label = c.name
    elif c_is_time:
        colorbar_label = '$t$ / fs'

    cmap = cmap or mpl.colormaps['cividis_r']
    if isinstance(cmap, str):
        cmap = mpl.colormaps[cmap]
    cnorm = cnorm or mpl.colors.Normalize(c.min(), c.max())  # type: ignore
    cscale = cscale or mpl.cm.ScalarMappable(  # type: ignore
        norm=cnorm, cmap=cmap
    )

    # TODO: remove groupby? Needed only for line-plot or for legend
    # for trajid, traj in noodle.groupby('trajid'):
    #     ctraj = c.sel(trajid=trajid)
    noodle_kws = noodle_kws or {}
    noodle_kws = {'alpha': 0.5, 's': 0.2, **noodle_kws}
    ax.scatter(noodle.isel(PC=0), noodle.isel(PC=1), c=cmap(cnorm(c)), **noodle_kws)

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    if hops is not None:
        hops_kws = dict(s=0.5, c='limegreen') | (hops_kws or {})
        ax.scatter(hops.isel(PC=0), hops.isel(PC=1), **hops_kws)

    # TODO facilitate custom colorbar
    fig.colorbar(cscale, ax=ax, label=colorbar_label, pad=0.02)

    # Alternative layout solution
    # d = make_axes_locatable(ax)
    # cax = d.append_axes("right", size="5%", pad="2%")
    # fig.colorbar(pc, cax=cax, label='dihedral')

    assert isinstance(ax, Axes)
    return ax


# TODO: finish later!
# def plot_noodleplot_lines(
#     noodle,  # hops,
#     ax=None,
#     cmap=None,
#     cnorm=None,
#     cscale=None,
# ):
#     fig, ax = figax(ax=ax)
#     points = noodle.values
#     # One traj per line
#     for trajid, traj in noodle.groupby('trajid'):
#         segments = np.concatenate([points[:-1], points[1:]], axis=1)
#         lc = mpl.collections.LineCollection([[0, 0]])

#     segments
#     lc
#     return ax


def get_loadings(frames, mean=False):
    atXYZ = frames['atXYZ']
    descr = get_standardized_pairwise_dists(atXYZ, mean=mean)
    _, pca_obj = pca(descr, 'atomcomb', return_pca_object=True)

    return xr.DataArray(
        data=pca_obj[-1].components_,
        dims=['PC', 'atomcomb'],
        coords=dict(atomcomb=descr.atomcomb),
        attrs={'natoms': frames.sizes['atom']},
    )


def plot_loadings(ax, loadings):
    for _, pcs in loadings.groupby('atomcomb'):
        assert len(pcs) == 2
        pc1, pc2 = pcs.item(0), pcs.item(1)
        ax.arrow(0, 0, pc1, pc2)
        a1, a2 = int(pcs['from']), int(pcs['to'])
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


def get_clusters_coords(loadings, atomcomb_clusters):
    return np.array(
        [
            loadings.isel(atomcomb=c).mean(dim='atomcomb').values
            for c in atomcomb_clusters
        ]
    )


def separate_angles(points, min_angle=10):
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


def filter_cluster_coords(coords, n):
    radii = [(x**2 + y**2) ** 0.5 for x, y in coords]
    angles = [np.degrees(np.arctan2(x, y)) for x, y in coords]
    res = set(np.argsort(radii)[-(n - 2) :])
    avg = np.mean(angles)
    splay = [abs(avg - angle) for angle in angles]
    return res.union(np.argsort(splay)[-2:])


def plot_clusters2(
    ax, loadings, clusters, mol, min_angle=10, inset_scale=1, show_at_most=None
):
    points = get_clusters_coords(loadings, clusters)
    if show_at_most is not None:
        indices = filter_cluster_coords(points, show_at_most)
        # clusters = [c for i, c in enumerate(clusters) if i in indices]
        # points = [p for i, p in enumerate(points) if i in indices]
    else:
        indices = range(len(clusters))
    scalefactors = separate_angles(points, min_angle)

    # else:
    # indices = range(len(clusters))
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
        # display(Image(png))  # DEBUG
        mpl_imshow_png(iax, png)


def plot_clusters3(loadings, clusters, ax=None, labels=None, axs=None, mol=None):
    fig, ax = figax(ax=ax)
    if labels is None:
        labels = list('abcdefghijklmnopqrstuvwxyz')

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


##################
# New stuff for angle binning!


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


def plot_bin_edges(angles, radii, bins, edges, picks, ax, labels):
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
            - picks:
            - angles:
            - center:
            - radii:
            - bins:
            - edges:
    """
    loadings = get_loadings(frames, mean)
    clusters = cluster_loadings(loadings)
    points = get_clusters_coords(loadings, clusters)

    angles = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    radii = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
    center = stats.circmean(angles, high=180, low=-180)

    picks, bins, edges = binning_with_min_entries(
        nbins=nbins, angles=angles, center=center, radii=radii, return_bins_edges=True
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


def binning_with_min_entries(
    nbins,
    angles,
    center,
    radii,
    min_entries=4,
    max_attempts=10,
    return_bins_edges=False,
):
    attempts = 0
    bins, edges = circbins(angles=angles, nbins=nbins, center=center)

    # Repeat binning until all bins have at least 'min_entries' or exceed max_attempts
    while any(arr.size == 0 for arr in bins) and attempts < max_attempts:
        print(
            f"Less than {min_entries} directions found, procedure repeated with another binning."
        )
        nbins += 1  # Increase the number of bins
        bins, edges = circbins(angles, nbins, center=center)
        attempts += 1

    # If max attempts were reached without satisfying condition
    if attempts >= max_attempts:
        print(f"Max attempts ({max_attempts}) reached. Returning current bins.")

    picks = [b[np.argmax(radii[b])] for b in bins]

    if return_bins_edges:
        return picks, bins, edges
    else:
        return picks
