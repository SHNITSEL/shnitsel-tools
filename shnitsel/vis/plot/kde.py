import logging
from typing import Any, Iterable, Literal, Sequence

from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.figure import Figure, SubFigure
import numpy as np
from numpy.typing import ArrayLike
import rdkit
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

from shnitsel.analyze.hops import hops_mask_from_active_state
from shnitsel.analyze.pca import PCAResult, pca, pca_and_hops
from shnitsel.data.dataset_containers import Frames, Trajectory, wrap_dataset
from shnitsel.data.dataset_containers.multi_layered import MultiSeriesLayered
from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.state_selection import StateSelection, StateSelectionDescriptor
from shnitsel.filtering.structure_selection import (
    AngleDescriptor,
    BondDescriptor,
    DihedralDescriptor,
    PyramidsDescriptor,
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc_ import pyramids
from shnitsel.geo.geocalc_.angles import angle
from shnitsel.geo.geocalc_.dihedrals import dihedral
from shnitsel.geo.geocalc_.distances import distance

from rdkit.Chem import Mol

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
            logging.warning(f"No points in range {p1} < x < {p2} for KDE fit")
            # raise ValueError(f"No points in range {p1} < x < {p2}")
            kernels.append(None)
        else:
            try:
                kernels.append(stats.gaussian_kde(subset))
            except Exception as e:
                logging.warning(f"{e}")
                kernels.append(None)
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
        if k is None:
            Zs.append(None)
        else:
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
    geo_property: xr.DataArray | TreeNode[Any, xr.DataArray],
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
    pca_data_da = pca_data.projected_inputs

    # Convert data to flat formats for operations
    if isinstance(pca_data_da, TreeNode):
        pca_data_da = pca_data_da.as_stacked
    assert isinstance(pca_data_da, xr.DataArray)

    if isinstance(geo_property, TreeNode):
        geo_property = geo_property.as_stacked
    assert isinstance(geo_property, xr.DataArray)
    pca_data_da = pca_data_da.transpose(
        'frame', 'time', 'PC', missing_dims='ignore'
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
        if Z is None:
            continue

        if contour_fill:
            ax.contourf(xx, yy, Z, levels=contour_levels, colors=c, alpha=0.1)
        ax.contour(xx, yy, Z, levels=contour_levels, colors=c, linewidths=0.5)


def biplot_kde(
    frames: xr.Dataset | ShnitselDataset | TreeNode[Any, ShnitselDataset | xr.Dataset],
    *ids: int,
    pca_data: TreeNode[Any, PCAResult] | PCAResult | None = None,
    state_selection: StateSelection | StateSelectionDescriptor | None = None,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    mol: rdkit.Chem.Mol | None = None,
    geo_kde_ranges: Sequence[tuple[float, float]] | None = None,
    scatter_color_property: Literal['time', 'geo'] = 'time',
    geo_feature: BondDescriptor
    | AngleDescriptor
    | DihedralDescriptor
    | PyramidsDescriptor
    | None = None,
    geo_cmap: str | None = 'PRGn',  # any valid cmap type
    time_cmap: str | None = 'cividis',  # any valid cmap type
    contour_levels: int | list[float] | None = None,
    contour_colors: list[str] | None = None,
    contour_fill: bool = True,
    num_bins: Literal[1, 2, 3, 4] = 4,
    fig: Figure | None = None,
    center_mean: bool = False,
) -> Figure | Sequence[Figure]:
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
        This needs to correspond to the data that was the input to `pca_data` if that parameter is provided.
    *ids: int
        Indices for atoms to be used in `geo_feature` if `geo_feature` is not set. 
        Note that pyramidalization angles cannot reliably be provided in this format.
    pca_data : PCAResult, optional
        A PCA result to use for the analysis. If not provided, will perform PCA analysis based on `structure_selection` or a
        generic pairwise distance PCA on `frames`.
        Accordingly, if provided, the parameter `frames` needs to correspond to the input provided to obtain the value in `
    structure_selection: StructureSelection | StructureSelectionDescriptor, optional
        An optional selection of features/structure to use for the PCA analysis.
    geo_kde_ranges : Sequence[tuple[float, float]], optional
        A Sequence of tuples representing ranges. A KDE is plotted for each range, indicating the distribution of
        points for which the value of the geometry feature falls in that range.
        Default values are chosen depending on the type of feature that should be analyzed. 
    contour_levels :  int | list[float], optional
        Contour levels for the KDE plot. Either the number of contour levels as an int or the list of floating 
        point values at which the contour lines should be drawn. Defaults to [0.08, 1]. 
        This parameter is passed to matplotlib.axes.Axes.contour.
    scatter_color_property : {'time', 'geo'}, default='time'
        Must be one of 'time' or 'geo'. If 'time', the scatter-points will be colored based on the time coordinate;
        if 'geo', the scatter-points will be colored based on the relevant geometry feature (see above).
    geo_cmap : str, default = 'PRGn'
        The Colormap to use for the noodleplot, if ``scatter_color='geo'``; this also determines contour
        colors unless ``contour_colors`` is set.
    time_cmap : str, default = 'cividis'
        The Colormap to use for the noodleplot, if ``scatter_color='time'``.
    contour_fill : bool, default = True
        Whether to plot filled contours (``contour_fill=True``, uses ``ax.contourf``)
        or just contour lines (``contour_fill=False``, uses ``ax.contour``).
    contour_colors : list[str], optional
        An iterable (not a Colormap) of colours (in a format matplotlib will accept) to use for the contours.
        By default, the ``geo_cmap`` will be used; this defaults to 'PRGn'.
    num_bins : {1, 2, 3, 4}, default = 4
        number of bins to be visualized, must be an integer between 1 and 4
    fig : mpl.figure.Figure, optional
        matplotlib.figure.Figure object into which the plot will be drawn;
        if not provided, one will be created using ``plt.figure(layout='constrained')``
    center_mean : bool, default = False
        Flag whether PCA data should be mean-centered before analysis. Defaults to False.

    Returns
    -------
    Figure
        The single figure of the PCA result, if the PCA result was not provided as a tree or on-the go PCA did not yield a tree result.
    Sequence[Figure]
        The sequence of all figures, one for each individual PCA result if the provided or obtained PCA result was a tree structure.

    Notes
    -----
    * Computes a geometric property of the specified atoms across all frames.
    * Uses kernel density estimation (KDE) to analyze the distance distributions.
    * Performs PCA on trajectory pairwise distances and visualizes clustering of structural changes.
    * Produces a figure with PCA projection, cluster analysis, and KDE plots.
    """

    if pca_data is None:
        # prepare data
        pca_data = pca(
            frames, structure_selection=structure_selection, center_mean=center_mean
        )

    if pca_data is not None and isinstance(pca_data, TreeNode):

        def single_pca_map(x: TreeNode[Any, PCAResult]) -> TreeNode[Any, Figure] | None:
            assert x.is_leaf
            if not x.has_data:
                return None

            pca_path = x.path
            pca_res = x.data

            frame_input_data = pca_res.inputs

            # TODO: FIXME: Make sourcing of input frames more robust. Maybe keep actual inputs on result?

            try:
                input_path = x._parent.path if x._parent is not None else "."
                if input_path.startswith("/"):
                    input_path = "." + input_path
                frame_input_data = frames[input_path]
            except:
                pass

            fig: Figure = biplot_kde(
                frame_input_data,
                *ids,
                pca_data=pca_res,
                state_selection=state_selection,
                structure_selection=structure_selection,
                mol=mol,
                geo_kde_ranges=geo_kde_ranges,
                scatter_color_property=scatter_color_property,
                geo_feature=geo_feature,
                geo_cmap=geo_cmap,  # any valid cmap type
                time_cmap=time_cmap,  # any valid cmap type
                contour_levels=contour_levels,
                contour_colors=contour_colors,
                contour_fill=contour_fill,
                num_bins=num_bins,
                center_mean=center_mean,
            )  # type: ignore # For single PCA, we get single result.
            assert isinstance(fig, Figure)
            fig.suptitle("PCA:" + pca_path)
            return x.construct_copy(data=fig)

        mapped_biplots = pca_data.map_filtered_nodes(
            lambda x: x.is_leaf, single_pca_map
        )
        assert mapped_biplots is not None, (
            "Failed to apply biplot to individual results in tree structure. Was the tree empty?"
        )
        return list(mapped_biplots.collect_data())

    try:
        hops_mask = hops_mask_from_active_state(
            frames, hop_type_selection=state_selection
        )
    except:
        logging.warning("Could not obtain `hops` mask from `frames` input.")
        hops_mask = None

    if scatter_color_property not in {'time', 'geo'}:
        raise ValueError("`scatter_color` must be 'time' or 'geo'")

    if contour_levels is None:
        contour_levels = [0.08, 1]

    if isinstance(frames, TreeNode):
        tree_mode = True
    else:
        tree_mode = False

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

    d = pb.pick_clusters(pca_data, num_bins=num_bins, center_mean=center_mean)
    loadings, clusters, picks = d['loadings'], d['clusters'], d['picks']

    res_mol: Mol
    if mol is None:
        if (
            structure_selection is None
            or not isinstance(structure_selection, StructureSelection)
            or structure_selection.mol is None
        ):
            # if tree_mode:
            #     res_mol = list(
            #         frames.map_data(lambda x: (x.__mol.item() if "__mol" in x else None) if isinstance(x, xr.DataArray) else wrap_dataset(x).mol).collect_data()
            #     )[0]
            if tree_mode:
                res_mol = list(frames.map_data(construct_default_mol).collect_data())[0]
            else:
                # print(f"{frames=}")
                # wrapped_da = wrap_dataset(frames)
                # print(f"{wrapped_ds=}")
                # res_mol = wrapped_da.mol
                res_mol = construct_default_mol(frames)

        else:
            res_mol = Mol(structure_selection.mol)
    else:
        res_mol = mol

    res_mol = set_atom_props(
        res_mol, atomLabel=True, atomNote=[''] * res_mol.GetNumAtoms()
    )

    if scatter_color_property == 'time':
        noodleplot_c = None
        noodleplot_cmap = time_cmap
        kde_data = None
    elif scatter_color_property == 'geo':
        if geo_feature is None:
            # Try and use additional positional parameters.
            geo_feature = tuple(ids)

        assert geo_feature is not None and len(geo_feature) >= 2, (
            "If the scatter property is set to `geo`, the `geo_feature` parameter of `biplot_kde()` must not be None."
        )

        wrapped_ds = wrap_dataset(frames)
        colorbar_label = None

        match geo_feature:
            case (atc, (at1, at2, at3)):
                # compute pyramidalization as described by the center atom `atc` and the neighbor atoms `at1, at2, at3`
                geo_prop = pyramids.pyramidalization_angle(
                    wrapped_ds.positions, atc, at1, at2, at3, deg=True
                )
                if not geo_kde_ranges:
                    geo_kde_ranges = [(-90, -10), (-10, 10), (10, 90)]
                colorbar_label = f"pyr({atc}, ({at1}, {at2}, {at3}))/°"
            case (at1, at2):
                # compute distance between atoms at1 and at2
                geo_prop = distance(wrapped_ds.positions, at1, at2)
                if not geo_kde_ranges:
                    geo_kde_ranges = [(0, 3), (5, 100)]
                colorbar_label = (
                    f'dist({at1}, {at2}) / {geo_prop.attrs.get("units", "Bohr")}'
                )
            case (at1, at2, at3):
                # compute angle between vectors at1 - at2 and at2 - at3
                assert at3 is not None  # to satisfy the typechecker
                geo_prop = angle(wrapped_ds.positions, at1, at2, at3, deg=True)
                if not geo_kde_ranges:
                    geo_kde_ranges = [(0, 80), (110, 180)]
                colorbar_label = (
                    f'angle({at1}, {at2}, {at3}) / {geo_prop.attrs.get("units", "°")}'
                )
            case (at1, at2, at3, at4):
                # compute dihedral defined as angle between normals to planes (at1, at2, at3) and (at2, at3, at4)
                assert at3 is not None
                assert at4 is not None
                geo_prop = dihedral(wrapped_ds.positions, at1, at2, at3, at4, deg=True)
                if not geo_kde_ranges:
                    geo_kde_ranges = [(0, 80), (110, 180)]
                colorbar_label = f'dih({at1}, {at2}, {at3}, {at4}) / {geo_prop.attrs.get("units", "°")}'
            case _:
                raise ValueError(
                    "The value provided to `biplot_kde()` as a `geo_feature` tuple does not constitute a Feature descriptor"
                )
        kde_data = _fit_and_eval_kdes(pca_data, geo_prop, geo_kde_ranges, num_steps=100)
        noodleplot_c = geo_prop
        noodleplot_cmap = geo_cmap
    else:
        assert False

    # noodleplot_c, noodleplot_cmap = {
    #     'time': (None, time_cmap),
    #     'geo': (geo_prop, geo_cmap),
    # }[scatter_color_property]

    # noodle_cnorm = mpl.colors.Normalize(noodleplot_c.min(), noodle_c.max())
    # noodle_cscale = mpl.cm.ScalarMappable(norm=noodle_cnorm, cmap=noodle_cmap)
    pca_noodles: TreeNode[Any, xr.DataArray] | xr.DataArray
    pca_noodles = pca_data.projected_inputs

    # TODO: FIXME: Noodle plot seems to have issues with rendering when passed data as a tree?
    pb.plot_noodleplot(
        pca_noodles,
        hops_mask,
        c=noodleplot_c,
        cmap=noodleplot_cmap,
        # cnorm=noodle_cnorm,
        colorbar_label=colorbar_label,
        ax=pcaax,
        noodle_kws=dict(alpha=1, marker='.'),
        hops_kws=dict(c='r', s=0.2),
    )

    # in case more clusters were found than we have room for:
    picks = picks[:4]

    # print(pca_data.explain_loadings())

    pb.plot_clusters_grid(
        loadings,
        [clusters[i] for i in picks],
        ax=pcaax,
        axs=structaxs,
        mol=res_mol,
        labels=list('abcd'),
    )

    if contour_colors is None:
        contour_colors = plt.get_cmap(noodleplot_cmap)(
            np.linspace(0, 1, len(contour_levels))
        )

    if kde_data:
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

    # TODO: FIXME: Should this really return the KDE data?
    # return kde_data
    return fig


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
