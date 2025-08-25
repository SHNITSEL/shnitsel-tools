from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import MinMaxScaler
import xarray as xr

from .. import _state

def pca(
    da: xr.DataArray,
    dim: str,
    n_components: int = 2,
    return_pca_object: bool = False
) -> tuple[xr.DataArray, PCA] | xr.DataArray:
    """xarray-oriented wrapper around scikit-learn's PCA

    Parameters
    ----------
    da
        A DataArray with at least a dimension with a name matching `dim`
    dim
        The name of the dimension to reduce
    n_components, optional
        The number of principle components to return, by default 2
    return_pca_object, optional
        Whether to return the scikit-learn `PCA` object as well as the
        transformed data, by default False

    Returns
    -------
    pca_res
        A DataArray with the same dimensions as `da`, except for the dimension
        indicated by `dim`, which is replaced by a dimension `PC` of size `n_components`
    [pca_object]
        The trained PCA object produced by scikit-learn, if return_pca_object=True
    """
    scaled = xr.apply_ufunc(
      MinMaxScaler().fit_transform,
      da.transpose(..., dim)
    )
    
    pca_object = PCA(n_components=n_components)
    pca_object.fit(scaled)
    pca_res: xr.DataArray = xr.apply_ufunc(
        pca_object.transform,
        scaled,
        input_core_dims=[[dim]],
        output_core_dims=[['PC']],
    )
    loadings = xr.DataArray(
        pca_object.components_,
        coords=[pca_res.coords['PC'], da.coords[dim]]
    )
    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        accessor_object = getattr(pca_res, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.loadings = loadings
        accessor_object.pca_object = pca_object

    if return_pca_object:
        return (pca_res, pca_object)
    else:
        return pca_res

def lda(da, dim, cats, n_components=2):
    if isinstance(cats, str):
        cats_name = cats
        cats = da[cats]
    else:
        cats_name = cats.name

    scaled = xr.apply_ufunc(MinMaxScaler().fit_transform, da.transpose(..., dim))
    lda_object = LinearDiscriminantAnalysis(n_components=n_components)

    def fit_transform(X):
        # cats: nonlocal
        return lda_object.fit_transform(X=X, y=cats)

    lda_res: xr.DataArray = xr.apply_ufunc(
        fit_transform,
        scaled,
        input_core_dims=[[dim]],
        output_core_dims=[['PC']],
    )

    lda_res[cats_name] = cats

    scalings = xr.DataArray(
        lda_object.scalings_, coords=[da.coords[dim], lda_res.coords['PC']]
    )

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        accessor_object = getattr(lda_res, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.scalings = scalings
        accessor_object.lda_object = lda_object

    return lda_res



def pls(xda, yda, n_components=2, common_dim=None):
    pls_object = PLSRegression(n_components=n_components)
    if len(xda.dims) != 2:
        raise ValueError(
            "xda should have 2 dimensions, in fact it has "
            f"{len(xda.dims)}, namely {xda.dims}"
        )
    if len(yda.dims) != 2:
        raise ValueError(
            "yda should have 2 dimensions, in fact it has "
            f"{len(yda.dims)}, namely {yda.dims}"
        )
    if common_dim is None:
        common_dims = set(xda.dims).intersection(yda.dims)
        if len(common_dims) != 1:
            raise ValueError(
                f"xda and yda have {len(common_dims)} dimension names in "
                f"common, namely {common_dims}. Please specify which of these "
                "should NOT be reduced, using the 'common_dim' parameter."
            )
        
        common_dim = common_dims.pop()
    
    xdim = (set(xda.dims) - {common_dim}).pop()
    ydim = (set(yda.dims) - {common_dim}).pop()
    xscaled = xr.apply_ufunc(MinMaxScaler().fit_transform, xda.transpose(..., xdim))
    yscaled = xr.apply_ufunc(MinMaxScaler().fit_transform, yda.transpose(..., ydim))
    xres, yres = xr.apply_ufunc(
        pls_object.fit_transform,
        xscaled,  # xda,
        yscaled,  # yda,
        input_core_dims=[[xdim], [ydim]],
        output_core_dims=[['score'], ['score']],
    )
    xname = xda.name or 'x'
    yname = yda.name or 'y'
    pls_res = xr.Dataset({xname: xres, yname: yres})
    loadings = xr.Dataset(
        {
            xname: ((xdim, 'loading'), pls_object.x_loadings_),
            yname: ((ydim, 'loading'), pls_object.y_loadings_),
        },
        coords={xdim: xda.coords[xdim], ydim: yda.coords[ydim]}
    )

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        accessor_object = getattr(pls_res, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.loadings = loadings
        accessor_object.pls_object = pls_object
    
    return pls_res

def pls_ds(ds, xname, yname, n_components=2):
    return pls(ds[xname], ds[yname], n_components=n_components)