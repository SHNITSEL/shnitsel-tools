from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import MinMaxScaler
import xarray as xr

from shnitsel import _state


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
    n_components : int, optional
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
    pls_object = PLSRegression(n_components=n_components)
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
