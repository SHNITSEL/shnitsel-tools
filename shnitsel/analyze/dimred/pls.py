from typing import Generic, TypeVar
from sklearn.cross_decomposition import PLSRegression as sk_PLS
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


class PLSResult(
    Generic[OriginType, ResultType],
    DimRedResult[OriginType, ResultType],
):
    """Class to hold the results of a PLS analysis.

    Also retains input data as well as corresponding results of the PLS decomposition.
    Input and output types are parametrized to allow for tree structures to be accurately represented.

    Provides accessors for all result meta data as well as the method `project_array(data_array)` to
    project another array of appropriate shape with dimension `mapped_dimension` to the PLS components.

    Parameters
    ----------
    OriginType
        The type of the original intput data. Should either be xr.DataArray for simple types, meaning we were provided a feature array
        or a flat DataGroup with xr.DataArrays in its leaves for tree types.
    ResultType
        Matching structure to `OriginType` but with the projected PLS decomposed input data as data within it.
        Either an xr.DataArray or a DataGroup same as for `OriginType`.
    """

    _pls_object: sk_PLS

    def __init__(
        self,
        pls_inputs: OriginType,
        pls_dimension: DimName,
        pls_pipeline: Pipeline,
        pls_object: sk_PLS,
        pls_projected_inputs: ResultType,
    ):
        if isinstance(pls_inputs, xr.DataArray):
            assert isinstance(pls_projected_inputs, xr.DataArray), (
                "If inputs are provided as a single data array, the results must also be a single data array"
            )
            coord_initial = [
                pls_projected_inputs.coords['scalings'],
                pls_inputs.coords[pls_dimension],
            ]
            _pls_components = xr.DataArray(
                # TODO: FIXME: These components do not 1:1 correspond to LDA or PCA loadings due to different input formats (2 sets of data instead of 1)
                [pls_object.x_loadings_, pls_object.y_loadings_],
                coords=coord_initial,
            ).assign_coords(
                PLSResult.get_extra_coords_for_loadings(pls_inputs, pls_dimension)
            )
        elif isinstance(pls_inputs, TreeNode):
            assert isinstance(pls_projected_inputs, TreeNode), (
                "If inputs are provided in tree shape, the projected inputs must also be provided in tree shape."
            )
            inputs_collected = list(pls_inputs.collect_data())
            outputs_collected = list(pls_projected_inputs.collect_data())
            assert all(isinstance(x, xr.DataArray) for x in inputs_collected), (
                "Tree-shaped inputs of LDA are not of data type xr.DataArray"
            )
            assert all(isinstance(x, xr.DataArray) for x in outputs_collected), (
                "Tree-shaped results of LDA are not of data type xr.DataArray"
            )
            coord_initial = [
                outputs_collected[0].coords['scalings'],
                inputs_collected[0].coords[pls_dimension],
            ]

            _pls_components = xr.DataArray(
                # TODO: FIXME: These components do not 1:1 correspond to LDA or PCA loadings due to different input formats (2 sets of data instead of 1)
                [pls_object.x_loadings_, pls_object.y_loadings_],
                coords=coord_initial,
            ).assign_coords(
                PLSResult.get_extra_coords_for_loadings(
                    inputs_collected[0], pls_dimension
                )
            )
        else:
            raise ValueError(
                f"Unsupported input format: {type(pls_inputs)}. Supported are xarray.DataArray and TreeNode thereof."
            )

        self._pls_object = pls_object
        super().__init__(
            inputs=pls_inputs,
            mapped_dimension=pls_dimension,
            pipeline=pls_pipeline,
            loadings=_pls_components,
            component_dimension="scalings",
            projected_inputs=pls_projected_inputs,
        )

    @property
    def fitted_pls_object(self) -> sk_PLS:
        return self._pls_object

    @property
    def scalings(self) -> xr.DataArray:
        return self.components


def pls(
    xdata_array: xr.DataArray,
    ydata_array: xr.DataArray,
    n_components: int = 2,
    common_dim: str | None = None,
) -> xr.Dataset:
    """Performs the partial least squares analysis on the two data arrays provided as arguments and returns the
    requested number of resulting main components.

    Parameters
    ----------
    xdata_array : xr.DataArray
        First set of data. Shape should be (n_samples, n_features) with n_samples being the length of the common dimension with `ydata_array`.
    ydata_array : xr.DataArray
        Second set of data. Shape should be (n_samples, n_targets) with n_samples being the length of the common dimension with `xdata_array`.
    n_components : int, default = 2
        Number of most relevant main components that should be returned. Defaults to 2.
    common_dim : str | None, optional
        The common dimension which should not be reduced in the course of the analysis. Defaults to None and will attempt to find a single common dimension.

    Raises
    ------
    ValueError
        If either `xdata_array` or `ydata_array` do not have exactly 2 dimensions.
    ValueError
        If no common_dim was set and x and y data did not have exactly 1 dimension in common to allow for automatic identification of the common dimension.

    Returns
    -------
    xr.Dataset
        The dataset holding the results of the PLS analysis. Results will either be in variables with the same name as `xdata_array` or `ydata_array` or in variables `x` and `y` if the names on the respective array are not set.
    """
    if len(xdata_array.dims) != 2:
        raise ValueError(
            "xdata_array should have 2 dimensions, in fact it has "
            f"{len(xdata_array.dims)}, namely {xdata_array.dims}"
        )
    if len(ydata_array.dims) != 2:
        raise ValueError(
            "ydata_array should have 2 dimensions, in fact it has "
            f"{len(ydata_array.dims)}, namely {ydata_array.dims}"
        )
    if common_dim is None:
        common_dims = set(xdata_array.dims).intersection(ydata_array.dims)
        if len(common_dims) != 1:
            raise ValueError(
                f"xdata_array and ydata_array have {len(common_dims)} dimension names in "
                f"common, namely {common_dims}. Please specify which of these "
                "should NOT be reduced, using the 'common_dim' parameter."
            )

        common_dim = str(common_dims.pop())

    xdim = (set(xdata_array.dims) - {common_dim}).pop()
    ydim = (set(ydata_array.dims) - {common_dim}).pop()

    xscaled = xr.apply_ufunc(
        MinMaxScaler().fit_transform, xdata_array.transpose(..., xdim)
    )
    yscaled = xr.apply_ufunc(
        MinMaxScaler().fit_transform, ydata_array.transpose(..., ydim)
    )

    # Get the PLS Regression object and perform the regression
    pls_object = sk_PLS(n_components=n_components)
    xres, yres = xr.apply_ufunc(
        pls_object.fit_transform,
        xscaled,  # xdata_array,
        yscaled,  # yda,
        input_core_dims=[[xdim], [ydim]],
        output_core_dims=[['score'], ['score']],
    )
    xname = xdata_array.name or 'x'
    yname = ydata_array.name or 'y'
    pls_res = xr.Dataset({xname: xres, yname: yres})
    # TODO: What are these loadings? and why do we not optionally return them?
    loadings = xr.Dataset(
        {
            xname: ((xdim, 'loading'), pls_object.x_loadings_),
            yname: ((ydim, 'loading'), pls_object.y_loadings_),
        },
        coords={xdim: xdata_array.coords[xdim], ydim: ydata_array.coords[ydim]},
    )

    # TODO: FIXME: Document the purpose of this black magic.
    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        accessor_object = getattr(pls_res, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.loadings = loadings
        accessor_object.pls_object = pls_object

    return pls_res


def pls_ds(
    dataset: xr.Dataset,
    xname: str,
    yname: str,
    n_components: int = 2,
    common_dim: str | None = None,
) -> xr.Dataset:
    """Wrapper function  to perform partial least square analysis on two variables of a dataset

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset holding the variables to apply PLS to
    xname : str
        The name of the variable to use as the x data for the PLS.
    yname : str
        The name of the variable to use as the y data for the PLS.
    n_components : int, optional
        Number of most relevant main components that should be returned. Defaults to 2.
    common_dim : str | None, optional
        The common dimension which should not be reduced in the course of the analysis. Defaults to None and will attempt to find a single common dimension.

    Returns
    -------
    xr.Dataset
        The result of the call to `pls()`. Has the results of the PLS as variables in either the same names as `xname` and `yname` or in `x` and `y`.s
    """
    return pls(
        dataset[xname], dataset[yname], n_components=n_components, common_dim=common_dim
    )


PLS = pls
