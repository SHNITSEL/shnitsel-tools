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
from shnitsel.core.typedefs import DimName
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
from .kde import (
    _fit_and_eval_kdes,
    plot_distribution_on_mesh,
    plot_distribution_heatmap,
)


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
        raise ValueError(
            f"Unsupported coloring option `{scatter_color_property}` only supported options are `geo` or `time`."
        )

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
        plot_distribution_on_mesh(
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
