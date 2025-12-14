from shnitsel import _state
from shnitsel._contracts import needs
import xarray as xr

from shnitsel.analyze.generic import norm, center, subtract_combinations
from shnitsel.data.multi_indices import mdiff
from sklearn.decomposition import PCA as sk_PCA

from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline

from shnitsel.core.typedefs import AtXYZ


@needs(coords_or_vars={'atXYZ', 'astate'})
def pca_and_hops(frames: xr.Dataset) -> tuple[xr.DataArray, xr.DataArray]:
    """Get PCA points and info on which of them represent hops

    Parameters
    ----------
    frames
        A Dataset containing 'atXYZ' and 'astate' variables

    Returns
    -------
    pca_res
        The PCA-reduced pairwise interatomic distances
    hops_pca_coords
        `pca_res` filtered by hops, to facilitate marking hops when plotting

    """
    pca_res = pairwise_dists_pca(frames['atXYZ'])
    mask = mdiff(frames['astate']) != 0
    hops_pca_coords = pca_res[mask]
    return pca_res, hops_pca_coords


@needs(dims={'atom'})
def pairwise_dists_pca(atXYZ: AtXYZ, **kwargs) -> xr.DataArray:
    """PCA-reduced pairwise interatomic distances

    Parameters
    ----------
    atXYZ
        A DataArray containing the atomic positions;
        must have a dimension called 'atom'

    Returns
    -------
        A DataArray with the same dimensions as `atXYZ`, except for the 'atom'
        dimension, which is replaced by a dimension 'PC' containing the principal
        components (by default 2)
    """
    res = (
        atXYZ.pipe(subtract_combinations, 'atom')
        .pipe(center)
        .pipe(norm)
        .pipe(pca, 'atomcomb', **kwargs)
    )
    assert not isinstance(res, tuple)  # typing
    return res


def pca(
    da: xr.DataArray, dim: str, n_components: int = 2, return_pca_object: bool = False
) -> tuple[xr.DataArray, sk_PCA] | xr.DataArray:
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
    scaler = MinMaxScaler()
    pca_object = sk_PCA(n_components=n_components)

    pipeline = Pipeline([('scaler', scaler), ('pca', pca_object)])

    pca_res: xr.DataArray = xr.apply_ufunc(
        pipeline.fit_transform,
        da,
        input_core_dims=[[dim]],
        output_core_dims=[['PC']],
    )
    loadings = xr.DataArray(
        pipeline[-1].components_, coords=[pca_res.coords['PC'], da.coords[dim]]
    )
    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        accessor_object = getattr(pca_res, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.loadings = loadings
        accessor_object.pca_object = pipeline

    if return_pca_object:
        # Return only PCA part of pipeline for backward-compatibility
        return (pca_res, pipeline)
    else:
        return pca_res


# Alternative names
principal_component_analysis = pca
PCA = pca
