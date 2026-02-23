from functools import partial
import logging
from typing import TYPE_CHECKING, Any, Callable, Generic, TypeVar
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as sk_LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as sk_QDA
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
import xarray as xr

from shnitsel.core.typedefs import DimName
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.data_group import DataGroup
from shnitsel.data.tree.data_leaf import DataLeaf
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.support_functions import tree_zip

from .dim_red_result import PredictorDimRedResult

from shnitsel import _state


OriginType = TypeVar('OriginType')
ResultType = TypeVar('ResultType')
DataType = TypeVar('DataType')


class DAResult(
    Generic[OriginType, ResultType],
    PredictorDimRedResult[OriginType, ResultType, xr.DataArray],
):
    """Class to hold the results of an Discrimant Analysis.

    Also retains input data as well as corresponding results of the LDA/QDA decomposition.
    Input and output types are parametrized to allow for tree structures to be accurately represented.

    Provides accessors for all result meta data as well as the method `project_array(data_array)` to
    project another array of appropriate shape with dimension `mapped_dimension` to the LDA/QDA components.

    Parameters
    ----------
    OriginType
        The type of the original intput data. Should either be `xr.DataArray` for simple types, meaning we were provided a feature array
        or a flat `DataGroup` with `xr.DataArrays` in its leaves for tree types.
    ResultType
        Matching structure to `OriginType` but with the projected input data within it.
        Either an xr.DataArray or a `DataGroup` same as for `OriginType`.
    """

    _da_object: Any
    _da_categories: OriginType

    def __init__(
        self,
        dimred_label: str,
        da_inputs: OriginType,
        da_inputs_categories: OriginType,
        da_dimension: DimName,
        da_pipeline: Pipeline,
        da_object,
        da_projected_inputs: ResultType,
    ):
        if isinstance(da_inputs, xr.DataArray):
            assert isinstance(da_projected_inputs, xr.DataArray), (
                "If inputs are provided as a single data array, the results must also be a single data array"
            )
            coord_initial = {
                'scaling': da_projected_inputs.coords['scaling'],
                da_dimension: da_inputs.coords[da_dimension],
            }

            scalings = xr.DataArray(
                da_object.scalings_,
                dims=(da_dimension, 'scaling'),
                # TODO: FIXME: This is again lacking all other coordinates along the `dim` dimension.
                coords=coord_initial,
            )
            scalings_ids = list(range(scalings.sizes['scaling']))
            _da_components = scalings.assign_coords(
                self.get_extra_coords_for_loadings(da_inputs, da_dimension)
            )

        elif isinstance(da_inputs, TreeNode):
            assert isinstance(da_projected_inputs, TreeNode), (
                "If inputs are provided in tree shape, the projected inputs must also be provided in tree shape."
            )
            inputs_collected = list(da_inputs.collect_data())
            outputs_collected = list(da_projected_inputs.collect_data())
            assert all(isinstance(x, xr.DataArray) for x in inputs_collected), (
                "Tree-shaped inputs of LDA are not of data type xr.DataArray"
            )
            assert all(isinstance(x, xr.DataArray) for x in outputs_collected), (
                "Tree-shaped results of LDA are not of data type xr.DataArray"
            )
            coord_initial = {
                'scaling': outputs_collected[0].coords['scaling'],
                da_dimension: inputs_collected[0].coords[da_dimension],
            }
            scalings = xr.DataArray(
                da_object.scalings_,
                dims=(da_dimension, 'scaling'),
                coords=coord_initial,
            )
            scalings_ids = list(range(scalings.sizes['scaling']))
            _da_components = scalings.assign_coords(
                scalings=('scaling', scalings_ids)
            ).assign_coords(
                self.get_extra_coords_for_loadings(inputs_collected[0], da_dimension)
            )
        else:
            raise ValueError(
                f"Unsupported input format: {type(da_inputs)}. Supported are xarray.DataArray and TreeNode thereof."
            )

        self._da_categories = da_inputs_categories
        self._da_object = da_object
        super().__init__(
            dimred_label=dimred_label,
            inputs=da_inputs,
            mapped_dimension=da_dimension,
            pipeline=da_pipeline,
            loadings=_da_components,
            component_dimension="scaling",
            projected_inputs=da_projected_inputs,
        )

    @property
    def input_labels(self) -> OriginType:
        return self._da_categories

    @property
    def scalings(self) -> xr.DataArray:
        return self.components

    @property
    def fitted_da_object(self):
        return self._da_object


class LDAResult(
    Generic[OriginType, ResultType],
    DAResult[OriginType, ResultType],
):
    """Specialization of the DAResult for the Linear discriminant analysis"""

    _lda_object: sk_LDA

    def __init__(
        self,
        da_inputs: OriginType,
        da_inputs_categories: OriginType,
        da_dimension: DimName,
        da_pipeline: Pipeline,
        da_object: sk_LDA,
        da_projected_inputs: ResultType,
    ):
        self._lda_object = da_object
        super().__init__(
            dimred_label='LDA',
            da_inputs=da_inputs,
            da_inputs_categories=da_inputs_categories,
            da_dimension=da_dimension,
            da_pipeline=da_pipeline,
            da_object=da_object,
            da_projected_inputs=da_projected_inputs,
        )

    @property
    def fitted_lda_object(self) -> sk_LDA:
        return self._lda_object


class QDAResult(
    Generic[OriginType, ResultType],
    DAResult[OriginType, ResultType],
):
    """Specialization of the DAResult for the Quadratic discriminant analysis"""

    _qda_object: sk_QDA

    def __init__(
        self,
        da_inputs: OriginType,
        da_inputs_categories: OriginType,
        da_dimension: DimName,
        da_pipeline: Pipeline,
        da_object: sk_QDA,
        da_projected_inputs: ResultType,
    ):
        self._qda_object = da_object
        super().__init__(
            dimred_label='QDA',
            da_inputs=da_inputs,
            da_inputs_categories=da_inputs_categories,
            da_dimension=da_dimension,
            da_pipeline=da_pipeline,
            da_object=da_object,
            da_projected_inputs=da_projected_inputs,
        )

    @property
    def fitted_qda_object(self) -> sk_QDA:
        return self._qda_object


def lda_direct(
    data_array: xr.DataArray,
    categories: xr.DataArray,
    dim: DimName,
    n_components: int = 2,
) -> LDAResult[xr.DataArray, xr.DataArray]:
    """Linear discriminant analysis performed on the data in `data_array` along `dim` in a total of `n_components

    Parameters
    ----------
    data_array : xr.DataArray
        The data to perform LDA on.
    dim : DimName
        The dimension to perform LDA along. Will be transposed to the back of all dimensions to then perform the LDA.
    categories : xr.DataArray
        Categories, provided as a named xr.DataArray.
        Should have dimension `dim`.
    n_components : int, optional
        The number of best main components to retrieve eventually. Defaults to 2.

    Returns
    -------
    LDAResult
        The results of the LDA as a LDAResult wrapping the xr.DataArray with the categories written to a variable if they weren't there before.
    """

    input_features: xr.DataArray = data_array
    category_flags: xr.DataArray = categories

    input_features = input_features.transpose(..., dim)
    cats_name = categories.name
    category_flags = categories  # .transpose(..., dim) # Has no features dimension. Only category per frame

    num_categories = len(np.unique(category_flags.values))
    num_features = input_features.sizes[dim]

    max_components = min(num_features, num_categories - 1)

    if n_components > max_components:
        logging.error(
            f"LDA can only work on `min(num_features, num_categories-1)={max_components}` compontents, not `{n_components}`. Limiting number to maximum possible."
        )
        n_components = max_components

    # TODO: For the fit, we may need to apply a different approach than for transform. The category data is only needed for fit, not for transform.
    scaler = MinMaxScaler()
    lda_object = sk_LDA(n_components=n_components)
    # Fit the scalers first and get scaled inputs
    scaled_inputs = scaler.fit_transform(input_features.values)

    # Fit with category flags separately
    fitted_lda = lda_object.fit(X=scaled_inputs, y=category_flags)

    pipeline = Pipeline([('scaler', scaler), ('lda', fitted_lda)])

    projected_inputs: xr.DataArray = xr.apply_ufunc(
        pipeline.transform,
        input_features,
        input_core_dims=[[dim]],
        output_core_dims=[['scaling']],
    ).assign_coords({'scaling': list(range(n_components)), cats_name: category_flags})

    lda_res = LDAResult(
        da_inputs=input_features,
        da_inputs_categories=category_flags,
        da_dimension=dim,
        da_pipeline=pipeline,
        da_object=fitted_lda,
        da_projected_inputs=projected_inputs,
    )

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        # TODO: FIXME: Add some explanation of this magic
        accessor_object = getattr(projected_inputs, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.scalings = lda_res.components
        accessor_object.lda_object = lda_object

    return lda_res


def qda_direct(
    data_array: xr.DataArray,
    categories: xr.DataArray,
    dim: DimName,
    n_components: int = 2,
) -> QDAResult[xr.DataArray, xr.DataArray]:
    """Quadratic discriminant analysis performed on the data in `data_array` along `dim` in a total of `n_components

    Parameters
    ----------
    data_array : xr.DataArray
        The data to perform QDA on.
    dim : DimName
        The dimension to perform QDA along. Will be transposed to the back of all dimensions to then perform the LDA.
    categories : xr.DataArray
        Categories, provided as a named xr.DataArray.
        Should have dimension `dim`.
    n_components : int, optional
        The number of best main components to retrieve eventually. Defaults to 2.

    Returns
    -------
    QDAResult
        The results of the QDA as a QDAResult wrapping the `xr.DataArray` with the categories written to a variable if they weren't there before.
    """

    input_features: xr.DataArray = data_array
    category_flags: xr.DataArray = categories

    input_features = input_features.transpose(..., dim)
    cats_name = categories.name
    category_flags = categories  # .transpose(..., dim) # Has no features dimension. Only category per frame

    num_categories = len(np.unique(category_flags.values))
    num_features = input_features.sizes[dim]

    max_components = min(num_features, num_categories - 1)

    if n_components > max_components:
        logging.error(
            f"LDA can only work on `min(num_features, num_categories-1)={max_components}` compontents, not `{n_components}`. Limiting number to maximum possible."
        )
        n_components = max_components

    # TODO: For the fit, we may need to apply a different approach than for transform. The category data is only needed for fit, not for transform.
    scaler = MinMaxScaler()
    qda_object = sk_QDA(n_components=n_components)
    # Fit the scalers first and get scaled inputs
    scaled_inputs = scaler.fit_transform(input_features.values)

    # Fit with category flags separately
    fitted_qda = qda_object.fit(X=scaled_inputs, y=category_flags)

    pipeline = Pipeline([('scaler', scaler), ('lda', fitted_qda)])

    projected_inputs: xr.DataArray = xr.apply_ufunc(
        pipeline.transform,
        input_features,
        input_core_dims=[[dim]],
        output_core_dims=[['scaling']],
    ).assign_coords({'scaling': list(range(n_components)), cats_name: category_flags})

    qda_res = QDAResult(
        da_inputs=input_features,
        da_inputs_categories=category_flags,
        da_dimension=dim,
        da_pipeline=pipeline,
        da_object=fitted_qda,
        da_projected_inputs=projected_inputs,
    )

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        # TODO: FIXME: Add some explanation of this magic
        accessor_object = getattr(projected_inputs, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.scalings = qda_res.components
        accessor_object.qda_object = qda_object

    return qda_res


def _da_shared_core(
    data: TreeNode[Any, ShnitselDataset | xr.Dataset | xr.DataArray]
    | ShnitselDataset
    | xr.Dataset
    | xr.DataArray,
    categories: TreeNode[Any, xr.DataArray] | str | xr.DataArray,
    dim: DimName | None = None,
    features: str | None = None,
    n_components: int = 2,
    da_result_class: type[LDAResult] | type[QDAResult] = LDAResult,
    da_label: str = 'LDA',
    da_direct_func: Callable = lda_direct,
) -> (
    TreeNode[Any, QDAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
    | QDAResult[xr.DataArray, xr.DataArray]
    | TreeNode[Any, LDAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
    | LDAResult[xr.DataArray, xr.DataArray]
    | None
):
    """
    Base function of the shared core code between LDA and QDA fitting.

    Parameters
    ----------
    data_array : xr.Dataset
        The data to perform LDA/QDA on.
    dim : str
        The dimension to perform LDA/QDA along. Will be transposed to the back of all dimensions to then perform the LDA/QDA. this should lead to
    cats : str | xr.DataArray
        Categories, either provided as the name of the variable in dataset where they are stored or as a named
        xr.DataArray. Should have dimension `dim`.
    n_components : int, optional
        The number of best main components to retrieve eventually. Defaults to 2.

    Returns
    -------
    LDA/QDAResult
        The results of the LDA/QDA as a LDA/QDAResult wrapping the xr.DataArray with the categories written to a variable if they weren't there before.
    """
    if isinstance(data, TreeNode):
        assert isinstance(categories, (TreeNode, str)), (
            f"If data is provided as a tree structure to {da_label}, the categories must be provided as a str key or as a Tree of xr.DataArray."
        )
        input_feature_tree: TreeNode[Any, xr.DataArray]
        category_flag_tree: TreeNode[Any, xr.DataArray]

        if isinstance(categories, str):
            cats_name = categories
            category_flag_tree = data.map_data(
                lambda x: x[categories], dtype=xr.DataArray
            )
        else:
            category_flag_tree = categories

        def extract_feature_array(
            entry: ShnitselDataset | xr.Dataset | xr.DataArray,
        ) -> xr.DataArray:
            if isinstance(entry, xr.DataArray):
                return entry

            assert features is not None and (
                features in data.data_vars or features in data.coords
            ), (
                "If provided an `xarray.Dataset` as `data` argument, the `features` key must be set to specify which data in the dataset to use."
            )

            return entry[features]

        # If we have arrays, we use them, otherwise we pick the `features` entry of datasets
        input_feature_tree = data.map_data(extract_feature_array)

        # Now we have a tree of inputs and a tree of flags. Zip them to assure we have the same structure and coorespoding data

        data_flag_tree: TreeNode[Any, tuple[xr.DataArray, xr.DataArray]] = tree_zip(
            input_feature_tree, category_flag_tree
        )  # type: ignore # Type checking with tree_zip is fickle... :/

        def filter_flat_group(node: TreeNode) -> bool:
            # We only want to process flat groups
            if isinstance(node, DataGroup):
                return node.is_flat_group
            return False

        def lda_on_flat_group(
            flat_group: TreeNode[Any, tuple[xr.DataArray, xr.DataArray]],
        ) -> DataGroup:
            assert isinstance(flat_group, DataGroup)
            assert flat_group.is_flat_group, (
                f"Something went wrong filtering for only flat groups in {da_label}"
            )

            group_features: DataGroup[xr.DataArray] = flat_group.map_data(
                lambda x: x[0], dtype=xr.DataArray
            )
            group_categories: DataGroup[xr.DataArray] = flat_group.map_data(
                lambda x: x[1], dtype=xr.DataArray
            )
            # Concatenate features
            glued_features: xr.DataArray = group_features.as_stacked
            glued_categories: xr.DataArray = group_categories.as_stacked

            assert isinstance(glued_features, xr.DataArray)

            # Perform concatenated PCA
            tmp_res = da_direct_func(
                data_array=glued_features,
                categories=glued_categories,
                dim=dim,
                n_components=n_components,
            )

            # Calculate hierarchical projected results
            mapped_inputs: DataGroup[xr.DataArray] = flat_group.map_data(
                lambda x: tmp_res.project_array(x[0])
            )  # type: ignore # Result should not be none here.

            # Rebuild appropriately shaped results with subtree as input and mapped
            # subtree as output
            full_res = da_result_class(
                da_inputs=group_features,
                da_inputs_categories=group_categories,
                da_dimension=tmp_res.mapped_dimension,
                da_pipeline=tmp_res.pipeline,
                da_object=tmp_res.fitted_da_object,
                da_projected_inputs=mapped_inputs,
            )

            # Build new subtree with only pca result
            new_leaf = DataLeaf(name=da_label.lower(), data=full_res)
            new_group = flat_group.construct_copy(children={da_label.lower(): new_leaf})  # type: ignore # with this copy construction we have the right data type
            return new_group

        lda_res = data_flag_tree.map_filtered_nodes(
            filter_flat_group, lda_on_flat_group
        )
        return lda_res
    else:
        assert isinstance(categories, (xr.DataArray, str)), (
            f"If data is provided as a simple structure to {da_label}, the categories must be provided as a str key or as a xr.DataArray."
        )
        input_features: xr.DataArray
        category_flags: xr.DataArray

        if dim is not None:
            data = data.transpose(..., dim)
        else:
            dim = list(data.sizes.keys())[-1]

        if TYPE_CHECKING:
            assert isinstance(data, (xr.DataArray | xr.Dataset))

        if isinstance(data, xr.Dataset):
            assert features is not None and (
                features in data.data_vars or features in data.coords
            ), (
                "If provided an `xarray.Dataset` as `data` argument, the `features` key must be set to specify which data in the dataset to use."
            )

            input_features = data[features]

            if isinstance(categories, str):
                assert features is not None and (
                    features in input_features.data_vars
                    or features in input_features.coords
                ), (
                    "The `categories` key does not identify data within the dataset or array to use for categories flags."
                )
                cats_name = categories
                category_flags = data[categories]
                # TODO: FIXME: I assume we need to delete the cats entry in the dataset here?
            else:
                cats_name = categories.name
                category_flags = categories
        else:
            input_features = data

            if isinstance(categories, str):
                assert categories is not None and (
                    features in input_features.data_vars
                    or features in input_features.coords
                ), (
                    f"The `{categories=}` key must be set to specify which data in the data source to use. "
                    "Was not found in data parameter. Please provide the categories directly as an array of the same size as `data`."
                )
                cats_name = categories
                category_flags = data[categories]
            else:
                cats_name = categories.name
                category_flags = categories

        return da_direct_func(
            data_array=input_features,
            categories=category_flags,
            dim=dim,
            n_components=n_components,
        )


def lda(
    data: TreeNode[Any, ShnitselDataset | xr.Dataset | xr.DataArray]
    | ShnitselDataset
    | xr.Dataset
    | xr.DataArray,
    categories: TreeNode[Any, xr.DataArray] | str | xr.DataArray,
    dim: DimName | None = None,
    features: str | None = None,
    n_components: int = 2,
) -> (
    TreeNode[Any, LDAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
    | LDAResult[xr.DataArray, xr.DataArray]
    | None
):
    """Linear discriminant analysis performed on the data in `data_array` along `dim` in a total of `n_components

    Parameters
    ----------
    data_array : xr.Dataset
        The data to perform LDA on.
    dim : str
        The dimension to perform LDA along. Will be transposed to the back of all dimensions to then perform the LDA. this should lead to
    cats : str | xr.DataArray
        Categories, either provided as the name of the variable in dataset where they are stored or as a named
        xr.DataArray. Should have dimension `dim`.
    n_components : int, optional
        The number of best main components to retrieve eventually. Defaults to 2.

    Returns
    -------
    LDAResult
        The results of the LDA as a LDAResult wrapping the xr.DataArray with the categories written to a variable if they weren't there before.
    """
    return _da_shared_core(
        data=data,
        categories=categories,
        dim=dim,
        features=features,
        n_components=n_components,
        da_result_class=LDAResult,
        da_label='LDA',
        da_direct_func=lda_direct,
    )  # type: ignore # For the LDA case, the appropriate LDA results are returned


# Alternative names for the analysis function
linear_discriminat_analysis = lda
LDA = lda


def qda(
    data: TreeNode[Any, ShnitselDataset | xr.Dataset | xr.DataArray]
    | ShnitselDataset
    | xr.Dataset
    | xr.DataArray,
    categories: TreeNode[Any, xr.DataArray] | str | xr.DataArray,
    dim: DimName | None = None,
    features: str | None = None,
    n_components: int = 2,
) -> (
    TreeNode[Any, QDAResult[DataGroup[xr.DataArray], DataGroup[xr.DataArray]]]
    | QDAResult[xr.DataArray, xr.DataArray]
    | None
):
    """Quadratic discriminant analysis performed on the data in `data_array` along `dim` in a total of `n_components

    Parameters
    ----------
    data_array : xr.Dataset
        The data to perform QDA on.
    dim : str
        The dimension to perform QDA along. Will be transposed to the back of all dimensions to then perform the QDA. this should lead to
    cats : str | xr.DataArray
        Categories, either provided as the name of the variable in dataset where they are stored or as a named
        xr.DataArray. Should have dimension `dim`.
    n_components : int, optional
        The number of best main components to retrieve eventually. Defaults to 2.

    Returns
    -------
    QDAResult
        The results of the QDA as a QDAResult wrapping the xr.DataArray with the categories written to a variable if they weren't there before.
    """
    return _da_shared_core(
        data=data,
        categories=categories,
        dim=dim,
        features=features,
        n_components=n_components,
        da_result_class=QDAResult,
        da_label='QDA',
        da_direct_func=qda_direct,
    )  # type: ignore # For the QDA case, the appropriate QDA results are returned


# Alternative names for the analysis function
quadratic_discriminat_analysis = qda
QDA = qda
