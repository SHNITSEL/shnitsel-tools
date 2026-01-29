import logging
from math import ceil
from typing import Any, Callable, Iterable, Sequence, TYPE_CHECKING, overload

from matplotlib.colors import Normalize, Colormap
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
from shnitsel.analyze.pca import PCAResult, pca
from shnitsel.data.dataset_containers import wrap_dataset, Trajectory, Frames
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.support_functions import tree_merge, tree_zip

from .common import figax, extrude, mpl_imshow_png
from ...rd import highlight_pairs

# if TYPE_CHECKING:
from rdkit.Chem import Mol


def plot_noodleplot(
    noodle: xr.DataArray | TreeNode[Any, xr.DataArray],
    hops_mask: xr.DataArray | TreeNode[Any, xr.DataArray] | None = None,
    fig: Figure | SubFigure | None = None,
    ax: Axes | None = None,
    c: NDArray | xr.DataArray | TreeNode[Any, xr.DataArray] | None = None,
    colorbar_label: str | None = None,
    cmap: str | Colormap | None = None,
    cnorm: Normalize | None = None,
    cscale=None,
    noodle_kws: dict | None = None,
    hops_kws: dict | None = None,
    rasterized: bool = True,
) -> Axes:
    """Create a `noodle` plot, i.e. a line or scatter plot of PCA-decomposed data.

    Parameters
    ----------
    noodle : xr.DataArray | TreeNode[Any, xr.DataArray]
        PCA decomposed data.
    hops_mask : xr.DataArray | TreeNode[Any, xr.DataArray], optional
        DataArray holding hopping-point information of the trajectories. Defaults to None.
    fig : Figure | SubFigure | None, optional
        Figure to plot the graph into. Defaults to None.
    ax : Axes, optional
        The axes to plot into. Will be generated from `fig` if not provided. Defaults to None.
    c : xr.DataArray | TreeNode[Any, xr.DataArray], optional
        The data to use for assigning the color to each individual data point. Defaults to None.
    colorbar_label : str | None, optional
        Label to plot next to the colorbar. If not provided will wither be taken from the `long_name` attribute or `name` attribute of the data or defaults to `t/fs`.
    cmap : str | Colormap | None, optional
        Colormap for plotting the datapoints. Defaults to None.
    cnorm : Normalize | None, optional
        Normalization method to map data to the colormap. Defaults to None.
    cscale : _type_, optional)
        The colorbar scale mapping that is used for creating the colorbar gradient. Defaults to None.
    noodle_kws : dict, optional
        Keywords arguments for the noodle/PCA plot. Defaults to None.
    hops_kws : dict, optional
        Keyword arguments for plotting the hopping points. Defaults to None.
    rasterized : bool, optional
        Flag to control whether the plot will be rasterized. Defaults to True.

    Returns
    -------
    Axes
        The `:py:class:matplotlib.axes.Axes` after plotting to them
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
        colorbar_label = str(c.attrs['long_name'])
        if 'units' in c.attrs:
            colorbar_label = f"${colorbar_label}$ / {c.attrs['units']}"
    elif hasattr(c, 'name') and c.name is not None:
        colorbar_label = str(c.name)
        if 'units' in c.attrs:
            colorbar_label = f"{colorbar_label} / {c.attrs['units']}"
    elif c_is_time:
        colorbar_label = '$t$ / fs'

    cmap = cmap or mpl.colormaps['cividis_r']
    if isinstance(cmap, str):
        cmap = mpl.colormaps[cmap]

    # TODO: remove groupby? Needed only for line-plot or for legend
    # for trajid, traj in noodle.groupby('trajid'):
    #     ctraj = c.sel(trajid=trajid)
    noodle_kws = noodle_kws or {}
    noodle_kws = {'alpha': 0.5, 's': 0.2, **noodle_kws}
    if isinstance(noodle, TreeNode):
        noodle_scatter = noodle.as_stacked
    else:
        noodle_scatter = noodle

    if isinstance(c, TreeNode):
        color_scatter = c.as_stacked
    else:
        color_scatter = c

    cnorm = cnorm or mpl.colors.Normalize(color_scatter.min(), color_scatter.max())  # type: ignore

    assert isinstance(noodle_scatter, xr.DataArray)
    assert isinstance(color_scatter, xr.DataArray)

    sc = ax.scatter(
        noodle_scatter.isel(PC=0).values,
        noodle_scatter.isel(PC=1).values,
        c=color_scatter,
        cmap=cmap,
        norm=cnorm,
        rasterized=rasterized,
        **noodle_kws,
    )

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    if hops_mask is not None:
        hops_kws = dict(s=0.5, c='limegreen') | (hops_kws or {})
        if isinstance(noodle, TreeNode):
            assert isinstance(hops_mask, TreeNode), (
                "If pca projections are tree data, the hops mask must be a tree of equal shape."
            )

            def apply_mask(x: tuple[xr.DataArray, xr.DataArray]) -> xr.DataArray | None:
                if any(x[1].values):
                    return (x[0])[(x[1])]
                return None

            merged_tree = tree_zip(noodle, hops_mask, res_data_type=tuple)

            hops_noodle = merged_tree.map_data(
                apply_mask, keep_empty_branches=False
            ).as_stacked
        else:
            assert not isinstance(hops_mask, TreeNode), (
                "If pca projections are no tree, the hops mask must be not be a tree either."
            )

            hops_noodle = noodle[hops_mask]

        ax.scatter(
            hops_noodle.isel(PC=0).values,
            hops_noodle.isel(PC=1).values,
            rasterized=rasterized,
            **hops_kws,
        )

    fig.colorbar(sc, ax=ax, label=colorbar_label, pad=0.02)

    # Alternative layout solution
    # d = make_axes_locatable(ax)
    # cax = d.append_axes("right", size="5%", pad="2%")
    # fig.colorbar(pc, cax=cax, label='dihedral')

    assert isinstance(ax, Axes)
    return ax


# TODO: implement plotting of noodleplot using multi-coloured lines


def get_loadings(
    frames_or_pca_result: xr.Dataset | Frames | Trajectory | PCAResult,
    center_mean: bool = False,
) -> xr.DataArray:
    """Get the loadings for the PCA of pairwise distances
    for the positional data in ``frames``.

    Parameters
    ----------
    frames : xr.Dataset | Frames | Trajectory
        A Dataset with an 'atXYZ' data_var, which should have
        'atom' and 'direction' dimensions.
    center_mean :  bool, optional
        Whether  centering of the mean should be should be applied, by default `False`

    Returns
    -------
    xr.DataArray
        A DataArray of loadings with dimensions
        'PC' (principal component) and 'descriptor' (atom combination,
        one for each pair of atoms).
    """
    pca_res: PCAResult
    attrs = {}
    if not isinstance(frames_or_pca_result, PCAResult):
        wrapped_ds = wrap_dataset(frames_or_pca_result, Frames | Trajectory)

        descr = get_standardized_pairwise_dists(
            frames_or_pca_result, center_mean=center_mean
        )
        pca_res = pca(descr, 'descriptor')
        attrs = {'natoms': wrapped_ds.sizes['atom']}
    else:
        pca_res = frames_or_pca_result

    return xr.DataArray(
        # data=pca_obj[-1].components_,
        data=pca_res.principal_components,
        dims=['PC', 'descriptor'],
        coords=dict(
            descriptor=pca_res.inputs.descriptor,
            descriptor_tex=pca_res.inputs.descriptor_tex,
            feature_indices=pca_res.inputs.feature_indices,
        ),
        attrs=attrs,
    )


def plot_loadings(ax: Axes, loadings: xr.DataArray):
    """Plot all loadings as arrows.

    Parameters
    ----------
    ax : Axes
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        the loadings.
    loadings: xr.DataArray
        A DataArray of PCA loadings including an 'descriptor' dimension;
        as produced by :py:func:`shnitsel.vis.plot.pca_biplot.get_loadings`.
    """
    # TODO: FIXME: This needs to be reconciled ith pairwise distances using StructureSelection now.
    raise NotImplementedError("Descriptor decomposition not implemented")
    for _, pcs in loadings.groupby('descriptor'):
        assert len(pcs) == 2
        pc1, pc2 = pcs.item(0), pcs.item(1)
        ax.arrow(0, 0, pc1, pc2)
        a1, a2 = int(pcs['atomcomb_from']), int(pcs['atomcomb_to'])
        ax.text(pc1, pc2, f"{a1},{a2}")


def cluster_general(decider: Callable[[int, int], bool], n: int) -> list[list[int]]:
    """Cluster indices iteratively according to a provided function.

    Parameters
    ----------
    decider : Callable[[int, int], bool]
        A function to decide whether two points can potentially share
        a cluster.
    n : int
        The number of indices to cluster.

    Returns
    -------
    list[list[int]]
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
    loadings : xr.DataArray
        A DataArray of loadings
    cutoff : float, optional
        An upper bound on the possible distances between a point
        in a cluster and other points, within which they will still
        be assigned to the smae cluster, by default 0.05

    Returns
    -------
    list[list[int]]
        A list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``.
    """

    def dist(i, j, l):
        pc1, pc2 = l.isel(descriptor=j).values - l.isel(descriptor=i).values
        return (pc1**2 + pc2**2) ** 0.5

    def decider(i, j):
        nonlocal loadings, cutoff, dist
        return dist(i, j, loadings) <= cutoff

    n = loadings.sizes['descriptor']
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
        A DataArray of PCA loadings including an 'descriptor' dimension;
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
        acs = loadings.isel(descriptor=cluster)
        x, y = acs.mean(dim='descriptor')
        s = (
            labels[i]
            if labels is not None
            else ' '.join([f'({a1},{a2})' for a1, a2 in acs.descriptor.values])
        )
        ax.arrow(0, 0, x, y)
        ax.text(x, y, s)


def _get_clusters_coords(loadings, descriptor_clusters):
    return np.array(
        [
            loadings.isel(descriptor=c).mean(dim='descriptor').values
            for c in descriptor_clusters
        ]
    )


def _separate_angles(points: NDArray, min_angle: float = 10) -> dict[int, float]:
    """Group points based on their polar angles, and work out scale factors
    by which to place labels along the ray from origin to point when annotating
    points, intending to avoid overlaps between labels.

    Parameters
    ----------
    points : NDArray
        An array of shape (npoints, 2)
    min_angle : float, optional
        The minimal difference in argument (angle from positive x-axis, in degrees),
        of two points, below which they will be considered part of the
        same cluster; by default 10

    Returns
    -------
    dict[int, float]
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
    ax: Axes,
    loadings: xr.DataArray,
    clusters: list[list[int]],
    mol: Mol,
    min_angle: float = 10,
    inset_scale: float = 1,
    show_at_most: int | None = None,
):
    """Plot selected clusters of the loadings of a pairwise distance PCA,
    and interpretations of those loadings, as highlighted molecular structures inset upon
    the loadings plot.

    Parameters
    ----------
    ax : Axes
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        the loadings
    loadings : xr.DataArray
        A DataArray of PCA loadings including an 'descriptor' dimension;
        as produced by :py:func:`shnitsel.vis.plot.pca_biplot.get_loadings`.
    clusters : list[list[int]]
        A list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``; as produced
        by :py:func:`shnitsel.vis.plot.pca_biplot.get_clusters`.
    mol : Mol
        An RDKit ``Mol`` object to be used for structure display.
    min_angle : float, optional
        Where multiple clusters of loadings lie in similar directions from
        the origin, they will be grouped together and only their member with the
        greatest radius will be annotated with a highlighted structure.
        This is the angle in degrees for the grouping behavior, by default 10.
    inset_scale : float, optional
        A factor by which to scale the size of the inset highlighted structures.
    show_at_most : int, optional
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
        acs = loadings.isel(descriptor=cluster)
        x, y = acs.mean(dim='descriptor')
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

        png = highlight_pairs(mol, acs.feature_indices.values)
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
    axs: dict[str, Axes] | None = None,
    mol: Mol | None = None,
):
    """Plot selected clusters of the loadings of a pairwise distance PCA,
    and interpretations of those loadings:

        - On the left, a large plot of selected clusters of loadings indicated as arrows
        - On the right, a grid of structures corresponding to
        structures of loadings; the pairs involved in the cluster
        are represented by colour-coding the atoms of the structures.

    Parameters
    ----------
    loadings : xr.DataArray
        A DataArray of PCA loadings including an 'descriptor' dimension;
        as produced by :py:func:`shnitsel.vis.plot.pca_biplot.get_loadings`.
    clusters : list[list[int]]
        A list of clusters, where each cluster is represented as a
        list of indices corresponding to ``loadings``; as produced
        by :py:func:`shnitsel.vis.plot.pca_biplot.get_clusters`.
    ax : Axes, optional
        The :py:class:`matplotlib.pyplot.axes.Axes` object onto which to plot
        the loadings
        (If not provided, one will be created.)
    labels : list[str], optional
        Labels for the loadings; if not provided, loadings will be labelled
        according to indices of the atoms to which they relate.
    axs : dict[str, Axes], optional
        A dictionary mapping from plot labels to :py:class:`matplotlib.pyplot.axes.Axes`
        objects
        (If not provided, one will be created.)
    mol : Mol, optional
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
        acs = loadings.isel(descriptor=cluster)
        x, y = acs.mean(dim='descriptor')
        s = labels[i]
        ax.arrow(0, 0, x, y, head_width=0.01, length_includes_head=True)

        x2, y2 = extrude(x, y, *ax.get_xlim(), *ax.get_ylim())

        ax.plot([x, x2], [y, y2], '--', c='k', lw=0.5)
        ax.text(x2, y2, s)

        if axs is not None and mol is not None:
            png = highlight_pairs(mol, acs.feature_indices.values)
            mpl_imshow_png(axs[s], png)
            axs[s].set_title(s)


# Compatability with old notebooks:
plot_clusters3 = plot_clusters_grid


def circbins(
    angles: np.ndarray,
    num_bins: int = 4,
) -> tuple[Sequence[np.ndarray], list[tuple[float, float]]]:
    """Bin angular data by clustering unit-circle projections

    Parameters
    ----------
    angles : np.ndarray
        Angles in degrees
    num_bins : int, optional
        Number of bins to return, by default 4

    Returns
    -------
    bins : Sequence[np.ndarray]
        Indices of angles belonging to each bin as an np.ndarray
    edges : list[tuple[float, float]]
        Tuple giving a pair of boundary angles for each bin;
        the order of the bins corresponds to the order used in ``bins``
    """

    def proj(x):
        "project angles in degrees onto unit circle"
        x = x * np.pi / 180
        return np.c_[np.cos(x), np.sin(x)]

    kmeans = KMeans(n_clusters=num_bins)
    labels = kmeans.fit_predict(proj(angles.astype(float)).astype(float))
    space = np.linspace(0.0, 360.0, num=10).astype(float)
    sample = kmeans.predict(proj(space).astype(float))
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
    angles : NDArray
        A 1D array of angles in degrees.
    radii : NDArray
        A 1D array of radii, with order corresponding
        to ``angles``.
    bins : list[Iterable[int]]
        Lists of bins, each bin represented as a list
        of indices.
    edges : list[tuple[float, float]]
        A pair of edges (angles in degrees) for each bin
        in ``bins``.
    picks : list[int]
        A list of indices indicating which cluster has
        been chosen from each bin.
    ax : Axes
        An matplotlib ``Axes`` object onto which to plot;
        this should be set up with polar projection.
    labels : list[str]
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


@overload
def pick_clusters(
    frames: Frames | xr.Dataset | PCAResult,
    num_bins: int,
    center_mean: bool = False,
) -> dict: ...


@overload
def pick_clusters(
    frames: TreeNode[Any, PCAResult],
    num_bins: int,
    center_mean: bool = False,
) -> TreeNode[Any, dict]: ...


# TODO: FIXME: This function does too much at once. It should either allow for general PCA configuration or allow for PCA results to be passed in. Only allowing pwdist pca is quite restrictive
def pick_clusters(
    frames: Frames | xr.Dataset | PCAResult | TreeNode[Any, PCAResult],
    num_bins: int,
    center_mean: bool = False,
) -> dict | TreeNode[Any, dict]:
    """Calculate pairwise-distance PCA, cluster the loadings
    and pick a representative subset of the clusters.

    Parameters
    ----------
    frames : Frames | xr.Dataset | PCAResult
        An :py:class:`xarray.Dataset` with an 'atXYZ' variable
        having an 'atom' dimension to calculate a pwdist PCA on or the result of a previously executed PCA.
    num_bins : int
        The number of bins to use when binning clusters of
        loadings according to the angle they make to the x-axis
        on the projection manifold
    center_mean : bool, optional
        Flag to apply mean centering before the analysis, by default Faule

    Returns
    -------
    dict
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
    TreeNode[Any, dict]
        If provided with a tree as input, this is returned per input leaf as a tree again
    """
    if isinstance(frames, TreeNode):
        return frames.map_data(pick_clusters)

    if isinstance(frames, PCAResult):
        loadings = frames.loadings
    else:
        wrapped_ds = wrap_dataset(frames, Frames | Trajectory)
        loadings = get_loadings(wrapped_ds, center_mean)
    clusters = cluster_loadings(loadings)
    points = _get_clusters_coords(loadings, clusters)

    angles = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    radii = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
    # center = stats.circmean(angles, high=180, low=-180)

    picks, bins, edges = _binning_with_min_entries(
        num_bins=num_bins, angles=angles, radii=radii, return_bins_edges=True
    )
    # bins, edges = circbins(angles, nbins=4, center=center)
    # picks = [b[np.argmax(radii[b])] for b in bins]

    return dict(
        loadings=loadings,
        clusters=clusters,
        picks=picks,
        angles=angles,
        radii=radii,
        bins=bins,
        edges=edges,
    )


def _binning_with_min_entries(
    num_bins: int,
    angles: NDArray,
    radii: NDArray,
    min_entries: int = 4,
    max_attempts: int = 10,
    return_bins_edges: bool = False,
) -> Sequence[int] | tuple[Sequence[int], Sequence[NDArray], list[tuple[float, float]]]:
    attempts = 0
    bins: Sequence[NDArray]
    edges: list[tuple[float, float]]
    bins, edges = circbins(angles=angles, num_bins=num_bins)

    # Repeat binning until all bins have at least 'min_entries' or exceed max_attempts
    while any(arr.size == 0 for arr in bins) and attempts < max_attempts:
        logging.info(
            f"Less than {min_entries} directions found, procedure repeated with another binning."
        )
        num_bins += 1  # Increase the number of bins
        bins, edges = circbins(angles, num_bins)
        attempts += 1

    # If max attempts were reached without satisfying condition
    if attempts >= max_attempts:
        logging.warning(
            f"Max attempts ({max_attempts}) reached. Returning current bins."
        )

    picks: Sequence[int] = [bin[np.argmax(radii[bin])] for bin in bins]

    if return_bins_edges:
        return picks, bins, edges
    else:
        return picks
