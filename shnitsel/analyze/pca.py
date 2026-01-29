from typing import Any, Generic, Hashable, Mapping, TypeVar, overload

import numpy as np

from shnitsel import _state
from shnitsel._contracts import needs
import xarray as xr

from shnitsel.analyze.generic import get_standardized_pairwise_dists
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.multi_indices import mdiff
from sklearn.decomposition import PCA as sk_PCA
from .hops import hops_mask_from_active_state

from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline

from shnitsel.data.tree.data_group import DataGroup
from shnitsel.data.tree.data_leaf import DataLeaf
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.structure_selection import (
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc import get_bats
from shnitsel.analyze.generic import norm

OriginType = TypeVar('OriginType')
ResultType = TypeVar('ResultType')
DataType = TypeVar('DataType')


class PCAResult(
    Generic[OriginType, ResultType],
):
    """Class to hold the results of a PCA analysis.

    Also retains input data as well as corresponding results of the PCA decomposition.
    Input and output types are parametrized to allow for tree structures to be accurately represented.

    Provides accessors for all result meta data as well as the method `project_array(data_array)` to
    project another array of appropriate shape with dimension `pca_mapped_dimension` to the PCA
    principal components.

    Parameters
    ----------
    OriginType
        The type of the original intput data. Should either be xr.DataArray for simple types, meaning we were provided a feature array
        or a flat DataGroup with xr.DataArrays in its leaves for tree types.
    ResultType
        Matching structure to `OriginType` but with the projected PCA decomposed input data as data within it.
        Either an xr.DataArray or a DataGroup same as for `OriginType`.
    """

    _pca_inputs: OriginType
    _pca_pipeline: Pipeline
    _pca_dimension: Hashable
    _pca_components: xr.DataArray
    _pca_object: sk_PCA
    _pca_inputs_projected: ResultType

    def __init__(
        self,
        pca_inputs: OriginType,
        pca_dimension: Hashable,
        pca_pipeline: Pipeline,
        pca_object: sk_PCA,
        pca_projected_inputs: ResultType,
    ):
        self._pca_inputs = pca_inputs
        self._pca_pipeline = pca_pipeline
        self._pca_dimension = pca_dimension
        if isinstance(pca_inputs, xr.DataArray):
            assert isinstance(pca_projected_inputs, xr.DataArray), (
                "If inputs are provided as a single data array, the results must also be a single data array"
            )
            coord_initial = [
                pca_projected_inputs.coords['PC'],
                pca_inputs.coords[pca_dimension],
            ]
            self._pca_components = xr.DataArray(
                pca_object.components_,
                coords=coord_initial,
            ).assign_coords(
                PCAResult.get_extra_coords_for_loadings(pca_inputs, pca_dimension)
            )
        elif isinstance(pca_inputs, TreeNode):
            assert isinstance(pca_projected_inputs, TreeNode), (
                "If inputs are provided in tree shape, the projected inputs must also be provided in tree shape."
            )
            inputs_collected = list(pca_inputs.collect_data())
            outputs_collected = list(pca_projected_inputs.collect_data())
            assert all(isinstance(x, xr.DataArray) for x in inputs_collected), (
                "Tree-shaped inputs of PCA are not of data type xr.DataArray"
            )
            assert all(isinstance(x, xr.DataArray) for x in outputs_collected), (
                "Tree-shaped results of PCA are not of data type xr.DataArray"
            )
            coord_initial = [
                outputs_collected[0].coords['PC'],
                inputs_collected[0].coords[pca_dimension],
            ]

            self._pca_components = xr.DataArray(
                pca_object.components_,
                coords=coord_initial,
            ).assign_coords(
                PCAResult.get_extra_coords_for_loadings(
                    inputs_collected[0], pca_dimension
                )
            )
        self._pca_object = pca_object
        self._pca_inputs_projected = pca_projected_inputs
        # TODO: Get the projected inputs?

    @property
    def inputs(self) -> OriginType:
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
    def loadings(self) -> xr.DataArray:
        return self._pca_components

    @property
    def projected_inputs(self) -> ResultType:
        return self._pca_inputs_projected

    @property
    def results(self) -> ResultType:
        return self.projected_inputs

    def get_most_significant_loadings(
        self, top_n_per: int = 5, top_n_total: int = 5
    ) -> tuple[Mapping[Hashable, xr.DataArray], xr.DataArray]:
        """Function to retrieve the most significant loadings in the
        PCA result for each individual component and in total.

        You can configure the amount of

        Parameters
        ----------
        top_n_per : int, optional
            Number of top (most significant absolute loading) n loadings per component, by default 5
        top_n_total : int, optional
            Number of overall top (i.e. most significant by 2-norm of their loadings across all PC) n features across all components, by default 5

        Returns
        -------
        tuple[Mapping[Hashable, xr.DataArray], xr.DataArray]
            First the mapping of each PC to the array holding the data of all their most significant loadings.
            Second the overall most significant loadings across all components.
        """
        loadings = self.loadings

        per_pc_results = {}
        for pc in loadings.PC.values:
            component = loadings.sel(PC=pc)
            # print(component)
            # print(component.values)

            top_n_per_local = min(component.sizes[self.pca_mapped_dimension], top_n_per)

            abs_loading = np.abs(component)
            top_arg_indices = np.argpartition(abs_loading, -top_n_per_local)[
                -top_n_per_local:
            ]
            top_arg_coords = component.coords[self.pca_mapped_dimension].values[
                top_arg_indices
            ]

            # print(top_arg_indices)
            # print(top_arg_coords)

            per_pc_results[pc] = component.sel(
                {self.pca_mapped_dimension: top_arg_coords}
            )

        top_n_total_local = min(loadings.sizes[self.pca_mapped_dimension], top_n_total)
        total_abs_loadings = norm(loadings, dim='PC')

        top_arg_indices = np.argpartition(total_abs_loadings, -top_n_total_local)[
            -top_n_total_local:
        ]
        top_arg_coords = loadings.coords[self.pca_mapped_dimension].values[
            top_arg_indices
        ]

        # print(top_arg_indices)
        # print(top_arg_coords)
        total_pc_results = loadings.sel({self.pca_mapped_dimension: top_arg_coords})
        # print(component.feature_indices)
        return per_pc_results, total_pc_results

    def explain_loadings(self, top_n_per: int = 5, top_n_total: int = 5) -> str:
        """Generate a textual explanation of the top influential loadings in the PCA result.

        Tries to put the results of `get_most_significant_loadings()` into a textual form.

        Parameters
        ----------
        top_n_per : int, optional
            Number of top (most significant absolute loading) n loadings per component, by default 5
        top_n_total : int, optional
            Number of overall top (i.e. most significant by 2-norm of their loadings across all PC) n features across all components, by default 5

        Returns
        -------
        str
            A text describing the results of the principal components analysis.
        """
        per_top, total_top = self.get_most_significant_loadings(
            top_n_per=top_n_per, top_n_total=top_n_total
        )

        explanation: str = ""

        total_expl = f"Maximum contributing features overall:\n"
        for feature, indices, coeff in zip(
            total_top.descriptor.values,
            total_top.feature_indices.values,
            norm(total_top, dim='PC').values,
        ):
            total_expl += f" {feature} (weight: {coeff}) (Idxs: {indices}) \n"
        explanation += total_expl + "\n\n"

        for pc in per_top:
            loadings = per_top[pc]

            pc_expl = f"Maximum contributing features to component {pc} :\n"
            for feature, indices, coeff in zip(
                loadings.descriptor.values,
                loadings.feature_indices.values,
                loadings.values,
            ):
                pc_expl += f" {feature}  (weight: {coeff}) (Idxs: {indices}) \n"
            explanation += pc_expl + "\n"
        return explanation

    def project_array(self, other_da: xr.DataArray) -> xr.DataArray:
        return xr.apply_ufunc(
            self._pca_pipeline.transform,
            other_da,
            input_core_dims=[[self._pca_dimension]],
            output_core_dims=[['PC']],
        )

    @staticmethod
    def get_extra_coords_for_loadings(
        data: xr.DataArray, dim: Hashable
    ) -> Mapping[Hashable, xr.DataArray]:
        # coords = {'PC': pca_res.coords['PC']}
        # coords.update(
        #     {
        #         key: coord
        #         for key, coord in data.coords.items()
        #         if dim in coord.dims and key != dim
        #     }
        # )
        coords = {
            key: coord
            for key, coord in data.coords.items()
            if dim in coord.dims and key != dim
        }
        return coords


@overload
def pca_and_hops(
    frames: TreeNode[Any, ShnitselDataset | xr.Dataset],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    center_mean: bool = False,
    n_components: int = 2,
) -> TreeNode[Any, tuple[PCAResult, xr.DataArray]]: ...


@overload
def pca_and_hops(
    frames: ShnitselDataset | xr.Dataset,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    center_mean: bool = False,
    n_components: int = 2,
) -> tuple[PCAResult, xr.DataArray]: ...


# TODO: Make signature consistent with `pca()` and standardize extraction of hops mask
@needs(coords_or_vars={'atXYZ', 'astate'})
def pca_and_hops(
    frames: TreeNode[Any, ShnitselDataset | xr.Dataset] | ShnitselDataset | xr.Dataset,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    center_mean: bool = False,
    n_components: int = 2,
) -> TreeNode[Any, tuple[PCAResult, xr.DataArray]] | tuple[PCAResult, xr.DataArray]:
    """
    Get PCA projected data and a mask to provide information on which of the data points represent hopping points.

    Parameters
    ----------
    frames : xr.Dataset | ShnitselDataset | TreeNode[Any, ShnitselDataset | xr.Dataset]
        A Dataset (or tree of those) containing 'atXYZ' and 'astate' variables
    structure_selection: StructureSelection | StructureSelectionDescriptor, optional
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

    # if isinstance(frames, TreeNode):
    #     def tmp_f(x:Frames | Trajectory | xr.Dataset) -> tuple[PCAResult, xr.DataArray]:
    #         return pca_and_hops(
    #                         x,
    #                         feature_selection=feature_selection,
    #                         center_mean=center_mean,
    #                         n_components=n_components,
    #                     )

    #     return frames.map_data(
    #         tmp_f,
    #         keep_empty_branches=True
    #     )

    if isinstance(frames, TreeNode):
        return frames.map_data(
            pca_and_hops,
            structure_selection=structure_selection,
            center_mean=center_mean,
            n_components=n_components,
        )

    wrapped_ds = wrap_dataset(frames, Frames | Trajectory | MultiSeriesDataset)
    assert isinstance(wrapped_ds, (Frames, Trajectory, MultiSeriesDataset)), (
        "provided frames data could not be considered trajectory or frameset data."
    )

    if structure_selection is None:
        # Will default to pairwise distances
        pca_res = pca(
            wrapped_ds,
            dim=None,
            structure_selection=None,
            n_components=n_components,
            center_mean=center_mean,
        )
    else:
        # Perform a pca with feature extraction
        pca_res = pca(
            wrapped_ds,
            dim=None,
            structure_selection=structure_selection,
            n_components=n_components,
            center_mean=center_mean,
        )

    hops_mask = hops_mask_from_active_state(wrapped_ds)

    return pca_res, hops_mask


@overload
def pca(
    data: TreeNode[Any, ShnitselDataset | xr.Dataset],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    dim: None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> TreeNode[Any, PCAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]:
    """Specialization for the pca being mapped over grouped data in a ShnitselDB tree structure"""
    ...


@overload
def pca(
    data: ShnitselDataset | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    dim: None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> PCAResult:
    """Specialization for the pca being applied to only simple structures"""
    ...


@overload
def pca(
    data: ShnitselDataset | xr.Dataset | TreeNode[Any, ShnitselDataset | xr.Dataset],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    dim: None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> (
    PCAResult
    | TreeNode[Any, PCAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
):
    """Perform a PCA decomposition on features derived from `data` using the structural features flagged in `structure_selection`.
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
    structure_selection :  StructureSelection | StructureSelectionDescriptor, optional
        Optional selection of geometric features to include in the PCA. If not provided,
        will fall back to pairwise distances.
    n_components : int, optional
        The number of principal components to be computed, by default 2
    center_mean : bool, optional
        Flag to center data before being passed to the PCA if set to `True`, by default `False`.

    Returns
    -------
    PCAResult
        The result of running the PCA analysis on the features selected in `structure_selection` or a full pairwise distance PCA
        extracted from `data`.
    ShnitselDB[PCAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
        The hierarchical structure of PCA results, where each flat group is used for a PCA analysis.
    """
    ...


@overload
def pca(
    data: xr.DataArray,
    structure_selection: None = None,
    dim: Hashable | None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> PCAResult[xr.DataArray, xr.DataArray]:
    """Perform a PCA decomposition directly on the data in `data` along the dimension `dim`
    without first deriving features from `data` to run the PCA on.

    Parameters
    ----------
    data : xr.DataArray
        The data for which the PCA should be calculated
    dim : Hashable
        The dimension along which the PCA should be performed
    structure_selection : None, optional
        Unused if `dim` is set, by default None
    n_components : int, optional
        The number of principal components to be computed, by default 2
    center_mean : bool, optional
        Flag to center data before being passed to the PCA if set to `True`, by default `False`.

    Returns
    -------
    PCAResult[xr.DataArray, xr.DataArray]
        The result of running the PCA analysis on the `data` array along the dimension `dim`.
    """
    ...


# TODO: FIXME: This should probably be pca_on_features and a separate pca_on_data() function to allo ShnitselDB trees with xr.DataArrays of only features calculated by the user
def pca(
    data: ShnitselDataset
    | xr.Dataset
    | xr.DataArray
    | TreeNode[Any, ShnitselDataset | xr.Dataset],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    dim: Hashable | None = None,
    n_components: int = 2,
    center_mean: bool = False,
) -> (
    PCAResult[xr.DataArray, xr.DataArray]
    | TreeNode[Any, PCAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
):
    """
    Function to perform a PCA decomposition on the `data` of various origins and formats.

    Can accept either full trajectory data in types of `Frames`, `Trajectory` or `ShnitselDB`
    hierarchical formats or as a raw `xr.Dataset`.
    Alternatively, the dataarray


    Parameters
    ----------
    da : xr.DataArray
        A DataArray with at least a dimension with a name matching `dim`
        dtype should be integer or floating with no
        ``nan`` or ``inf`` entries
    structure_selection :  StructureSelection | StructureSelectionDescriptor, optional
        Optional selection of geometric features to include in the PCA. If not provided,
        will fall back to pairwise distances.
    dim
        The name of the array-dimension to reduce (i.e. the axis along which different
        features lie)
    n_components : int, optional
        The number of principal components to return, by default 2
    center_mean : bool, optional
        Flag to center data before being passed to the PCA if set to `True`, by default `False`.

    Returns
    -------
    PCAResult[xr.DataArray, xr.DataArray]
        The full information obtained by the fitting of the result.
        Contains the inputs for the PCA result, the principal components,
        the mapped values for the inputs, the full pipeline to apply the PCA
        transformation again to other data.

        The mapped inputs are a DataArray with the same dimensions as ``da``, except for the dimension
        indicated by `dim`, which is replaced by a dimension ``PC`` of size ``n_components``.

        ``result.principal_components`` holds the fitted principal components.
        ``result.projected_inputs`` provides the PCA projection result when applied to the inputs.
    ShnitselDB[PCAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
        The hierarchical structure of PCA results, where each flat group is used for a PCA analysis.

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
        if isinstance(data, TreeNode):
            # TODO: FIXME: We need to catch the xr.DataArray tree input earlier. Remove dataarray tree input for now

            def traj_to_frame(
                x: ShnitselDataset | xr.Dataset | xr.DataArray,
            ) -> ShnitselDataset | xr.DataArray | None:
                if isinstance(x, xr.DataArray):
                    return x
                x = wrap_dataset(x)
                return x

                if isinstance(x, (Trajectory, Frames)):
                    return x.as_frames
                elif isinstance(x, MultiSeriesDataset):
                    return x.as_stacked
                else:
                    return None

            data_framed: TreeNode[Any, ShnitselDataset | xr.DataArray] = data.map_data(
                traj_to_frame
            )
            data_grouped = data_framed.group_data_by_metadata()
            assert data_grouped is not None

            if structure_selection is not None:

                def extract_features(x: ShnitselDataset | xr.DataArray) -> xr.DataArray:
                    return get_bats(
                        x,
                        structure_selection=structure_selection,  # deg='trig'
                    )
            else:

                def extract_features(x: ShnitselDataset | xr.DataArray) -> xr.DataArray:
                    return (
                        get_standardized_pairwise_dists(x, center_mean=center_mean)
                        # .swap_dims(atomcomb='descriptor')
                        # .rename(atomcomb='descriptor')
                    )

            # We extract the features either with the selection or with the
            # Pairwise distances approach
            feature_data_grouped = data_grouped.map_data(extract_features)

            def filter_flat_group(node: TreeNode) -> bool:
                # We only want to process flat groups
                if isinstance(node, DataGroup):
                    return node.is_flat_group
                return False

            def pca_on_flat_group(
                flat_group: TreeNode[Any, xr.DataArray],
            ) -> DataGroup[PCAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]:
                assert isinstance(flat_group, DataGroup)
                assert flat_group.is_flat_group, (
                    "Something went wrong filtering for only flat groups in PCA"
                )

                inputs: DataGroup[xr.DataArray] = flat_group
                # Extract feature arrays out of leaves
                collected_features = list(flat_group.collect_data())
                if collected_features and 'time' in collected_features[0].sizes:
                    leading_dim = 'time'
                else:
                    leading_dim = 'frame'

                # Concatenate features
                glued_features = (
                    inputs.as_stacked
                )  # xr.concat(collected_features, dim=leading_dim)

                assert isinstance(glued_features, xr.DataArray)

                # Perform concatenated PCA
                tmp_res = pca_direct(
                    glued_features, dim='descriptor', n_components=n_components
                )

                # Calculate hierarchical projected results
                mapped_inputs: DataGroup[xr.DataArray] = inputs.map_data(
                    tmp_res.project_array
                )  # type: ignore # Result should not be none here.

                # Rebuild appropriately shaped results with subtree as input and mapped
                # subtree as output
                full_res = PCAResult(
                    pca_inputs=inputs,
                    pca_dimension=tmp_res.pca_mapped_dimension,
                    pca_pipeline=tmp_res.pca_pipeline,
                    pca_object=tmp_res.fitted_pca_object,
                    pca_projected_inputs=mapped_inputs,
                )

                # Build new subtree with only pca result
                new_leaf = DataLeaf(name='pca', data=full_res)
                new_group: DataGroup[
                    PCAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]
                ] = flat_group.construct_copy(children={'pca': new_leaf})  # type: ignore # with this copy construction we have the right data type
                return new_group

            pca_res = feature_data_grouped.map_filtered_nodes(
                filter_flat_group, pca_on_flat_group
            )
            return pca_res
        else:
            feature_array: xr.DataArray
            # if isinstance(data, ShnitselDataset):
            #     # extract positional data
            #     data = data.positions
            # elif isinstance(data, xr.Dataset):
            #     # Extract positional data from `atXYZ` variable.
            #     data = data.atXYZ

            # At this point, we should have an instance of xr.DataArray in `data` or someone provided us with a weird
            # unsupported input type
            # if isinstance(data, xr.DataArray):
            # We got a `atXYZ` or positions array
            if structure_selection is None:
                # Get array with standardized pairwise distance features.
                # Need to rename to ensure the relevant dimension is called `descriptor` and not `atomcomb`
                feature_array = get_standardized_pairwise_dists(
                    data, center_mean=center_mean
                )  # .rename({'atomcomb': 'descriptor'})
            else:
                feature_array = get_bats(
                    data,
                    structure_selection=structure_selection,
                    # TODO: FIXME: Check if `trig` is the best option for us
                    # deg='trig',
                )

            return pca_direct(
                feature_array, dim='descriptor', n_components=n_components
            )
            # else:
            #     raise ValueError(
            #         "Provided instance of `data` could not be used to extract positional data "
            #         "(or the `atXYZ` variable specifically) required for feature calculation for the PCA."
            #     )


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

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        # TODO: Potentially remove? The new result type holds way more data/information

        coord_initial = [pca_res.coords['PC'], data.coords[dim]]
        # coords = {'PC': pca_res.coords['PC']}
        # coords.update(
        #     {
        #         key: coord
        #         for key, coord in data.coords.items()
        #         if dim in coord.dims and key != dim
        #     }
        # )
        coords = {
            key: coord
            for key, coord in data.coords.items()
            if dim in coord.dims and key != dim
        }

        loadings = xr.DataArray(
            pipeline[-1].components_, coords=coord_initial
        ).assign_coords(coords)

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
