import logging
from typing import Any, Iterable, Literal, Sequence, overload

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


def fit_kde(
    data: xr.DataArray, leading_dim: DimName | None = None
) -> stats.gaussian_kde | None:
    """\
    Fit a set of KDEs to the `data`.
    
    The parameter `geo_kde_ranges` specifies the subsets of the values of `geo_property`
    that should be filtered into the same subset. 
    Returns one KDE for each such subset.


    Parameters
    ----------
    data : xr.DataArray
        The data for which KDEs should be fitted on the various ranges. 
        Usually PCA data or data obtained from another kind of clustering or dimensionality reduction.
    leading_dim: DimName, optional
        Name of the leading dimension which enumerates the different samples. 
        Needs to be shifted to the back of the dataset. 
        If not provided, will be guessed from 'time' or 'frame'.

    Returns
    ----------
    stats.gaussian_kde
        The fitted KDE (kernel density estimator) for the entire dataset.
    None 
        If data was empty.
    """

    if leading_dim is None:
        if 'frame' in data.dims:
            leading_dim = 'frame'
        elif 'time' in data.dims:
            leading_dim = 'time'
        else:
            raise ValueError(
                "Could not guess leading dimension name from `data` input."
            )
    if data.size == 0:
        logging.warning(f"No points for KDE to fit")
        return None
    else:
        base = data.transpose(..., leading_dim) # Swap leading dimension to the end

        try:
            return stats.gaussian_kde(base)
        except Exception as e:
            logging.warning(f"Failed to fit KDE: {e}")
            return None


def filter_by_property_range(
    data: xr.DataArray,
    property_values: xr.DataArray,
    property_ranges: Sequence[tuple[float, float]],
    leading_dim: DimName | None = None,
) -> Sequence[xr.DataArray]:
    """\
    Filter data in a data array by the values of a property. 
    For each range in `property_ranges`, extract the set of data points in that range and
    returns a Sequence of `xarray.DataArray` instances holding the `data` points 
    with the property in that range.


    Parameters
    ----------
    data : xr.DataArray
        The data which should be filtered for the various ranges. 
        Usually PCA data or data obtained from another kind of clustering or dimensionality reduction.
    property_values : xr.DataArray
        The (geometric) property that the data should be clustered/filtered by.
    property_ranges : Sequence[tuple[float, float]]
        The sequence of (distinct) ranges of values of the (geometric) property
        that the `pca_data` should be divided by.
    leading_dim: DimName, optional
        Name of the leading dimension which enumerates the different samples. 
        Needs to be shifted to the back of the dataset. 
        If not provided, will be guessed from 'time' or 'frame'.

    Returns
    ----------
    Sequence[xr.DataArray]
        The sequence of subsets of data with `property_values` in the respective `property_range` intervals.
    """
    if leading_dim is None:
        if 'frame' in data.dims:
            leading_dim = 'frame'
        elif 'time' in data.dims:
            leading_dim = 'time'
        else:
            raise ValueError(
                "Could not guess leading dimension name from `data` input."
            )

    filtered_data = []
    for p1, p2 in property_ranges:
        mask = (p1 < property_values) & (property_values < p2)
        subset = data.sel({leading_dim: mask})
        if subset.size == 0:
            logging.warning(
                f"No data points in range {p1} < x < {p2} for property filter"
            )
        filtered_data.append(subset)
    return filtered_data


def fit_filtered_kdes(
    data: xr.DataArray,
    geo_property: xr.DataArray,
    geo_kde_ranges: Sequence[tuple[float, float]],
    leading_dim: DimName | None = None,
) -> Sequence[stats.gaussian_kde | None]:
    """\
    Fit a set of KDEs to the `data` (usuallyn PCA or other dimension reduced data), 
    after it has been split into subsets based on the values of
    `geo_property`. 
    
    The parameter `geo_kde_ranges` specifies the subsets of the values of `geo_property`
    that should be filtered into the same subset. 
    Returns one KDE for each such subset.


    Parameters
    ----------
    data : xr.DataArray
        The data for which KDEs should be fitted on the various ranges. 
        Usually PCA data or data obtained from another kind of clustering or dimensionality reduction.
    geo_property : xr.DataArray
        The geometric property that the data should be clustered/filtered by.
    geo_kde_ranges : Sequence[tuple[float, float]]
        The sequence of (distinct) ranges of values of the geometric property
        that the `pca_data` should be divided by.
    leading_dim: DimName, optional
        Name of the leading dimension which enumerates the different samples. 
        Needs to be shifted to the back of the dataset. 
        If not provided, will be guessed from 'time' or 'frame'.


    Returns
    ----------
    Sequence[stats.gaussian_kde | None]
        The sequence of fitted KDEs (kernels) for each range of `geo_kde_ranges`.
        If there were no points in the respective range, None is returned as the entry of that range.

    Raises
    ------
    ValueError
        If any of the ``geo_filter`` ranges is such that no points from
        ``geo_prop`` fall within it
    """
    if leading_dim is None:
        if 'frame' in data.dims:
            leading_dim = 'frame'
        elif 'time' in data.dims:
            leading_dim = 'time'
        else:
            raise ValueError(
                "Could not guess leading dimension name from `data` input."
            )

    filtered_data = filter_by_property_range(data, geo_property, geo_kde_ranges)

    kernels = [fit_kde(subset, leading_dim=leading_dim) for subset in filtered_data]
    return kernels


@overload
def _eval_kdes_on_mesh_grid(
    kernels: stats.gaussian_kde,
    xx: np.ndarray,
    yy: np.ndarray,
    normalize: bool = True,
) -> np.ndarray: ...


@overload
def _eval_kdes_on_mesh_grid(
    kernels: Sequence[stats.gaussian_kde | None],
    xx: np.ndarray,
    yy: np.ndarray,
    normalize: bool = True,
) -> Sequence[np.ndarray]: ...


def _eval_kdes_on_mesh_grid(
    kernels: stats.gaussian_kde | Sequence[stats.gaussian_kde | None],
    xx: np.ndarray,
    yy: np.ndarray,
    normalize: bool = True,
) -> np.ndarray | Sequence[np.ndarray | None]:
    """Evaluate all fitted gaussian kernel density estimators on a mesh-grid
    and return the results.


    Parameters
    ----------
        kernels : stats.gaussian_kde | Sequence[stats.gaussian_kde | None]
            The transformed data fitted with gaussian kernel estimators to
            evaluate on the mesh grid.
        xx : np.ndarray
            The x coordinates of the mesh grid.
        yy : np.ndarray
            The y coordinates of the mesh grid.
        normalize : bool, default=True
            Flag whether the evaluated kdes should be normalized by their maximum value.

    Returns
    ----------
        np.ndarray
            If a single estimator was provided, a single result will be returned
        Sequence[np.ndarray]
            The sequence of evaluated approximate probability densities
            at the positions described by `xx` and `yy` for each and every
            individual KDE provided in `kernels`.
    """
    if isinstance(kernels, Sequence):
        return [
            _eval_kdes_on_mesh_grid(kernels=kernel, xx=xx, yy=yy, normalize=normalize)
            if kernel is not None
            else None
            for kernel in kernels
        ]

    xys = np.c_[xx.ravel(), yy.ravel()].T
    if kernels is None:
        return None
    else:
        Z = kernels.evaluate(xys)
        Z = Z.reshape(xx.shape)
        if normalize:
            Z /= Z.max()
        return Z


def _get_oversized_meshgrid(
    data: xr.DataArray,
    num_steps: int = 500,
    extension: float = 0.1,
    leading_dim: DimName | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Get appropriately over-sized mesh-grids for x and y coordinates
    with an excess overhang of `extension` relative to the min/max-to-mean distance
    and `num_steps` intermediate steps between the upper and lower bound.

    Statistical properties will be derived from `data`.

    Note that `data` should be 2-dimensional except for `leading_dim`.


    Parameters
    ----------
    data: xr.DataArray
        The (pca) data to get the supporting mesh grid for.
    num_steps, optional : int, default=500
        Number of intermediate steps to generate in the grid. Defaults to 500.
    extension, optional : float, default=0.1
        Excess overhang beyond minima and maxima in x and y direction
        relative to their distance from the mean. Defaults to 0.1.
    leading_dim: DimName, optional
        Name of the leading dimension which enumerates the different samples.
        Needs to be shifted to the back of the dataset.
        If not provided, will be guessed from 'time' or 'frame'.

    Returns
    ----------
    tuple[np.ndarray, np.ndarray]
        First the numpy array holding x positions of a meshgrid
        Then the array holding y positions of a meshgrid.
    """
    if leading_dim is None:
        if 'frame' in data.dims:
            leading_dim = 'frame'
        elif 'time' in data.dims:
            leading_dim = 'time'
        else:
            raise ValueError(
                "Could not guess leading dimension name from `data` input."
            )

    means: np.ndarray = data.mean(dim=leading_dim).values
    mins: np.ndarray = data.min(dim=leading_dim).values
    mins -= (means - mins) * extension
    maxs: np.ndarray = data.max(dim=leading_dim).values
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

    xx, yy = _get_oversized_meshgrid(
        pca_data_da, num_steps=num_steps, extension=extension
    )
    kernels = fit_filtered_kdes(pca_data_da, geo_property, geo_kde_ranges)
    return xx, yy, _eval_kdes_on_mesh_grid(kernels, xx, yy)


def plot_distribution_on_mesh(
    xx: np.ndarray,
    yy: np.ndarray,
    Zs: np.ndarray | Sequence[np.ndarray],
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
    Zs : np.ndarray | Sequence[np.ndarray]
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
        if isinstance(Zs, np.ndarray) or len(Zs) == 1:
            colors = ['purple']
        elif len(Zs) == 2:
            colors = ['purple', 'green']
        else:
            colors = plt.get_cmap('inferno')(range(len(Zs)))

    if isinstance(Zs, Sequence):
        for Z, c in zip(Zs, colors):
            if Z is None:
                continue
            if contour_fill:
                ax.contourf(xx, yy, Z, levels=contour_levels, colors=c, alpha=0.1)
            ax.contour(xx, yy, Z, levels=contour_levels, colors=c, linewidths=0.5)
    else:
        if Zs is not None:
            if contour_fill:
                ax.contourf(
                    xx, yy, Zs, levels=contour_levels, colors=colors[0], alpha=0.1
                )
            ax.contour(
                xx, yy, Zs, levels=contour_levels, colors=colors[0], linewidths=0.5
            )
    return fig, ax


def plot_distribution_heatmap(
    xx: np.ndarray,
    yy: np.ndarray,
    Zs: np.ndarray,
    cmap: str | Colormap | None = None,
    fig: Figure | SubFigure | None = None,
    ax: Axes | None = None,
):
    """Plots kernel density estimate as a heatmap

    Parameters
    ----------
    xx : np.ndarray
        An array of x values
    yy : np.ndarray
        An array of y values (must have the same shape as ``xx``)
    Zs : np.ndarray
        The density information sampled on the meshgrid denoted by `xx`, `yy`.
    cmap : str | Colormap , optional
        The colormap to use to plot the heatmap
    fig : Figure | SubFigure, optional
        A matplotlib ``Figure`` object into which to draw
        (if not provided, a new one will be created)
    ax : Axes, optional
        A matplotlib ``Axes`` object into which to draw
        (if not provided, a new one will be created)
    """
    fig, ax = figax(fig=fig, ax=ax)
    if cmap is None:
        cmap = plt.get_cmap('inferno')

    ax.pcolormesh(xx, yy, Zs, cmap=cmap, rasterized=True)
    return fig, ax


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
