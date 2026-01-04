from typing import Hashable, overload

from shnitsel import _state
from shnitsel._contracts import needs
import xarray as xr

from shnitsel.analyze.generic import get_standardized_pairwise_dists
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.multi_indices import mdiff
from sklearn.decomposition import PCA as sk_PCA

from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline

from shnitsel.core.typedefs import AtXYZ
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.structure_selection import StructureSelection
from shnitsel.geo.geocalc import get_bats


# TODO: Make inputs type configurable.
class PCAResult:
    _pca_inputs: xr.DataArray
    _pca_pipeline: Pipeline
    _pca_dimension: Hashable
    _pca_components: xr.DataArray
    _pca_object: sk_PCA
    _pca_inputs_projected: xr.DataArray

    def __init__(
        self,
        pca_inputs: xr.DataArray,
        pca_dimension: Hashable,
        pca_pipeline: Pipeline,
        pca_object: sk_PCA,
        pca_projected_inputs: xr.DataArray,
    ):
        self._pca_inputs = pca_inputs
        self._pca_pipeline = pca_pipeline
        self._pca_dimension = pca_dimension
        self._pca_components = xr.DataArray(
            pca_object.components_,
            coords=[
                pca_projected_inputs.coords['PC'],
                pca_inputs.coords[pca_dimension],
            ],
        )
        self._pca_object = pca_object
        self._pca_inputs_projected = pca_projected_inputs
        # TODO: Get the projected inputs?

    @property
    def inputs(self) -> xr.DataArray:
        return self._pca_inputs

    @property
    def fitted_pca_object(self) -> sk_PCA:
        return self._pca_object

    @property
    def pca_mapped_dimension(self) -> Hashable:
        return self._pca_dimension

    @property
    def pca_pipeline(self) -> Pipeline:
        return self._pca_pipeline

    @property
    def principal_components(self) -> xr.DataArray:
        return self._pca_components

    @property
    def projected_inputs(self) -> xr.DataArray:
        return self._pca_inputs_projected

    def project_array(self, other_da: xr.DataArray) -> xr.DataArray:
        return xr.apply_ufunc(
            self._pca_pipeline.transform,
            other_da,
            input_core_dims=[[self._pca_dimension]],
            output_core_dims=[['PC']],
        )


# TODO: Make signature consistent with `pca()` and standardize extraction of hops mask
@needs(coords_or_vars={'atXYZ', 'astate'})
def pca_and_hops(
    frames: Frames | Trajectory | xr.Dataset,
    feature_selection: StructureSelection | None = None,
    center_mean: bool = False,
    n_components: int = 2,
) -> tuple[PCAResult, xr.DataArray]:
    """
    Get PCA points and info on which of them represent hops

    _extended_summary_

    Parameters
    ----------
    frames : xr.Dataset
        A Dataset containing 'atXYZ' and 'astate' variables
    feature_selection : StructureSelection, optional
        An optional selection of features to calculate and base the PCA fitting on.
        If not provided, will calculate a PCA for full pairwise distances.
    center_mean : bool
        Center mean data before pca if True, by default: False.
    n_components : int, optional
        The number of principal components to return, by default 2, by default 2

    Returns
    -------
    tuple[PCAResult, xr.DataArray]
        A tuple of the following two parts:
        - pca_res
            The object result of the call to `pca()` holding all results of the pca analysis (see documentation of `pca()`).
        - hopping_point_masks
            The mask of the hopping point events. Can be used to only extract the hopping point PCA results from the projected input result in pca_res.
    """
    if feature_selection is None:
        pca_res = pairwise_dists_pca(
            frames['atXYZ'], n_components=n_components, center_mean=center_mean
        )
    else:
        pca_res = pca_with_features(
            frames.atXYZ, feature_selection=feature_selection, n_components=n_components
        )

    hops_positions = mdiff(frames['astate']) != 0

    return pca_res, hops_positions


@overload
def pca(
    data: ShnitselDB[Trajectory | Frames],
    dim: None = None,
    feature_selection: StructureSelection | None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> ShnitselDB[PCAResult]:
    """Specialization for the pca being mapped over grouped data in a ShnitselDB tree structure"""
    ...


@overload
def pca(
    data: Trajectory
    | Frames
    | ShnitselDB[Trajectory | Frames]
    | xr.Dataset
    | xr.DataArray,
    dim: None = None,
    feature_selection: StructureSelection | None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> PCAResult | ShnitselDB[PCAResult]:
    """Perform a PCA decomposition on features derived from `data` using the structural features flagged in `feature_selection`.
    Will not directly run PCA on the input data.

    Yiealds

    Parameters
    ----------
    data : Trajectory | Frames | ShnitselDB[Trajectory | Frames] | xr.Dataset | xr.DataArray
        The data for which the features for the PCA should be calculated.
        Assumes that either the format provides `position` information, can be converted to a `Trajectory` or `Frames`
        instance or that it is an `atXYZ` data array holding positional information from which the features
        for the PCA can be calculated.
    dim : None
        Unused in this specialization. If `dim` is not None with `data` of type other than `xr.DataArray`, an exception will be Raised.
    feature_selection : StructureSelection, optional
        The flagged geometric features that should be calculated from `data`, by default None.
        If set to None, full pairwise distances between all positions in the `data` will be
    n_components : int, optional
        The number of principal components to be computed, by default 2
    center_mean : bool, optional
        Flag to center data before being passed to the PCA if set to `True`, by default `False`.

    Returns
    -------
    PCAResult
        The result of running the PCA analysis on the features selected in `feature_selection` or a full pairwise distance PCA
        extracted from `data`.
    """
    ...


@overload
def pca(
    data: xr.DataArray,
    dim: Hashable,
    feature_selection: None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> PCAResult:
    """Perform a PCA decomposition directly on the data in `data` along the dimension `dim`
    without first deriving features from `data` to run the PCA on.

    Parameters
    ----------
    data : xr.DataArray
        The data for which the PCA should be calculated
    dim : Hashable
        The dimension along which the PCA should be performed
    feature_selection : None, optional
        Unused if `dim` is set, by default None
    n_components : int, optional
        The number of principal components to be computed, by default 2
    center_mean : bool, optional
        Flag to center data before being passed to the PCA if set to `True`, by default `False`.

    Returns
    -------
    PCAResult
        The result of running the PCA analysis on the `data` array along the dimension `dim`.
    """
    ...


def pca(
    data: Trajectory
    | Frames
    | ShnitselDB[Trajectory | Frames]
    | xr.Dataset
    | xr.DataArray,
    dim: Hashable | None = None,
    feature_selection: StructureSelection | None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> PCAResult | ShnitselDB[PCAResult]:
    """
    Function to perform a PCA decomposition on the `data` of various origins and formats.

    Can accept either full trajectory data in types of `Frames`, `Trajectory` or `ShnitselDB`
    hierarchical formats or as a raw `xr.Dataset`.
    Alternatively, the dataarray


    Parameters
    ----------
    da : xr.DataArray
        A DataArray with at least a dimension with a name matching `dim`
    dim : Hashable
        The name of the array-dimension to reduce (i.e. the axis along which different
        features lie)
    n_components : int, optional
        The number of principal components to return, by default 2
    center_mean : bool, optional
        Flag to center data before being passed to the PCA if set to `True`, by default `False`.

    Returns
    -------
    PCAResult
        The full information obtained by the fitting of the result.
        Contains the inputs for the PCA result, the principal components,
        the mapped values for the inputs, the full pipeline to apply the PCA
        transformation again to other data.

        The mapped inputs are a DataArray with the same dimensions as ``da``, except for the dimension
        indicated by `dim`, which is replaced by a dimension ``PC`` of size ``n_components``.

        ``result.principal_components`` holds the fitted principal components.
        ``result.projected_inputs`` provides the PCA projection result when applied to the inputs.

    Examples:
    ---------
    >>> pca_results1 = pca(data1)
    >>> pca_results1.projected_inputs  # See the loadings
    >>> pca_results2 = pca_results1.project_array(data2)
    """
    if dim is not None:
        assert isinstance(data, xr.DataArray), (
            "If an analysis dimension `dim` is provided, the `data` parameter must be of type `xarray.DataArray`"
        )
        return pca_direct(data, dim=dim, n_components=n_components)
    else:
        # We need to calculate features first.

        if isinstance(data, ShnitselDB):
            raise NotImplementedError()
        else:
            feature_array: xr.DataArray
            if isinstance(data, Trajectory) or isinstance(data, Frames):
                # extract positional data
                data = data.positions
            elif isinstance(data, xr.Dataset):
                # Extract positional data from `atXYZ` variable.
                data = data.atXYZ

            # At this point, we should have an instance of xr.DataArray in `data` or someone provided us with a weird
            # unsupported input type
            if isinstance(data, xr.DataArray):
                # We got a `atXYZ` or positions array
                if feature_selection is None:
                    # Get array with standardized pairwise distance features.
                    # Need to rename to ensure the relevant dimension is called `descriptor` and not `atomcomb`
                    feature_array = get_standardized_pairwise_dists(
                        data, center_mean=center_mean
                    ).rename({'atomcomb': 'descriptor'})
                else:
                    feature_array = get_bats(
                        data, structure_selection=feature_selection
                    )

                return pca_direct(
                    feature_array, dim='descriptor', n_components=n_components
                )
            else:
                raise ValueError(
                    "Provided instance of `data` could not be used to extract positional data "
                    "(or the `atXYZ` variable specifically) required for feature calculation for the PCA."
                )


def pca_direct(data: xr.DataArray, dim: Hashable, n_components: int = 2) -> PCAResult:
    """Wrapper function to directly apply the PCA decomposition to the values in a dataarray.

    Contrary to the `pca()` function, the features for the pca are not derived from the first `data` parameter

    Parameters
    ----------
    data : xr.DataArray
        A DataArray with at least a dimension with a name matching `dim`
    dim : Hashable
        The name of the array-dimension to reduce (i.e. the axis along which different
        features lie)
    n_components : int, optional
        The number of principal components to return, by default 2

    Returns
    -------
    PCAResult
        The full information obtained by the fitting of the result.
        Contains the inputs for the PCA result, the principal components,
        the mapped values for the inputs, the full pipeline to apply the PCA
        transformation again to other data.

        The mapped inputs are a DataArray with the same dimensions as ``da``, except for the dimension
        indicated by `dim`, which is replaced by a dimension ``PC`` of size ``n_components``.

    Examples:
    ---------
    >>> pca_results1 = pca(data1, 'features')
    >>> pca_results1.projected_inputs  # See the loadings
    >>> pca_results2 = pca_results1.project_array(data2)
    """
    scaler = MinMaxScaler()
    pca_object = sk_PCA(n_components=n_components)

    pipeline = Pipeline([('scaler', scaler), ('pca', pca_object)])

    pca_res: xr.DataArray = xr.apply_ufunc(
        pipeline.fit_transform,
        data,
        input_core_dims=[[dim]],
        output_core_dims=[['PC']],
    )
    loadings = xr.DataArray(
        pipeline[-1].components_, coords=[pca_res.coords['PC'], data.coords[dim]]
    )

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        # TODO: Potentially remove? The new result type holds way more data/information

        def use_to_transform(other_da: xr.DataArray):
            return xr.apply_ufunc(
                pipeline.transform,
                other_da,
                input_core_dims=[[dim]],
                output_core_dims=[['PC']],
            )

        accessor_object = getattr(pca_res, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.loadings = loadings
        accessor_object.pca_object = pipeline
        accessor_object.use_to_transform = use_to_transform

    pca_result_wrapper = PCAResult(
        pca_inputs=data,
        pca_object=pipeline[-1],
        pca_dimension=dim,
        pca_projected_inputs=pca_res,
        pca_pipeline=pipeline,
    )

    return pca_result_wrapper


# Alternative names
principal_component_analysis = pca
PCA = pca
