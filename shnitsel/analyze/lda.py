from typing import Generic, TypeVar
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as sk_LDA
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
import xarray as xr

from shnitsel.core.typedefs import DimName
from shnitsel.data.tree.node import TreeNode

from .dim_red_result import DimRedResult

from shnitsel import _state


OriginType = TypeVar('OriginType')
ResultType = TypeVar('ResultType')
DataType = TypeVar('DataType')


class LDAResult(
    Generic[OriginType, ResultType],
    DimRedResult[OriginType, ResultType],
):
    """Class to hold the results of an LDA analysis.

    Also retains input data as well as corresponding results of the LDA decomposition.
    Input and output types are parametrized to allow for tree structures to be accurately represented.

    Provides accessors for all result meta data as well as the method `project_array(data_array)` to
    project another array of appropriate shape with dimension `mapped_dimension` to the LDA components.

    Parameters
    ----------
    OriginType
        The type of the original intput data. Should either be xr.DataArray for simple types, meaning we were provided a feature array
        or a flat DataGroup with xr.DataArrays in its leaves for tree types.
    ResultType
        Matching structure to `OriginType` but with the projected PCA decomposed input data as data within it.
        Either an xr.DataArray or a DataGroup same as for `OriginType`.
    """

    _lda_object: sk_LDA

    def __init__(
        self,
        lda_inputs: OriginType,
        lda_dimension: DimName,
        lda_pipeline: Pipeline,
        lda_object: sk_LDA,
        lda_projected_inputs: ResultType,
    ):
        if isinstance(lda_inputs, xr.DataArray):
            assert isinstance(lda_projected_inputs, xr.DataArray), (
                "If inputs are provided as a single data array, the results must also be a single data array"
            )
            coord_initial = [
                lda_projected_inputs.coords['scalings'],
                lda_inputs.coords[lda_dimension],
            ]
            _lda_components = xr.DataArray(
                lda_object.scalings_,
                coords=coord_initial,
            ).assign_coords(
                LDAResult.get_extra_coords_for_loadings(lda_inputs, lda_dimension)
            )
        elif isinstance(lda_inputs, TreeNode):
            assert isinstance(lda_projected_inputs, TreeNode), (
                "If inputs are provided in tree shape, the projected inputs must also be provided in tree shape."
            )
            inputs_collected = list(lda_inputs.collect_data())
            outputs_collected = list(lda_projected_inputs.collect_data())
            assert all(isinstance(x, xr.DataArray) for x in inputs_collected), (
                "Tree-shaped inputs of LDA are not of data type xr.DataArray"
            )
            assert all(isinstance(x, xr.DataArray) for x in outputs_collected), (
                "Tree-shaped results of LDA are not of data type xr.DataArray"
            )
            coord_initial = [
                outputs_collected[0].coords['scalings'],
                inputs_collected[0].coords[lda_dimension],
            ]

            _lda_components = xr.DataArray(
                lda_object.scalings_,
                coords=coord_initial,
            ).assign_coords(
                LDAResult.get_extra_coords_for_loadings(
                    inputs_collected[0], lda_dimension
                )
            )
        else:
            raise ValueError(
                f"Unsupported input format: {type(lda_inputs)}. Supported are xarray.DataArray and TreeNode thereof."
            )

        self._lda_object = lda_object
        super().__init__(
            inputs=lda_inputs,
            mapped_dimension=lda_dimension,
            pipeline=lda_pipeline,
            loadings=_lda_components,
            component_dimension="scalings",
            projected_inputs=lda_projected_inputs,
        )

    @property
    def fitted_lda_object(self) -> sk_LDA:
        return self._lda_object

    @property
    def scalings(self) -> xr.DataArray:
        return self.components


# TODO: FIXME: We have no example for how this is used, so documentation is purely guesswork.
def lda(
    data: xr.Dataset | xr.DataArray,
    categories: str | xr.DataArray,
    dim: str | None = None,
    features: str | None = None,
    n_components: int = 2,
) -> xr.DataArray:
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
    xr.DataArray
        The results of the LDA as a DataArray with the categories written to a variables if they weren't there before.
    """

    input_features: xr.DataArray
    category_flags: xr.DataArray

    if dim is not None:
        data = data.transpose(..., dim)

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
                "If provided an `xarray.Dataset` as `data` argument, the `features` key must be set to specify which data in the dataset to use."
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
        data,
        input_core_dims=[[dim]],
        output_core_dims=[['scaling']],
    )

    # TODO: FIXME: Build the LDA result and return that.

    projected_inputs[cats_name] = category_flags

    scalings = xr.DataArray(
        lda_object.scalings_,
        dims=(dim, 'scaling'),
        # TODO: FIXME: This is again lacking all other coordinates along the `dim` dimension.
        coords={dim: data.coords[dim], 'scaling': [i for i in range(n_components)]},
    )

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        # TODO: FIXME: Add some explanation of this magic
        accessor_object = getattr(projected_inputs, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.scalings = scalings
        accessor_object.lda_object = lda_object

    return projected_inputs


# Alternative names for the analysis function
linear_discriminat_analysis = lda
LDA = lda
