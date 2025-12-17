from shnitsel import _state
from shnitsel._contracts import needs
import xarray as xr

from shnitsel.analyze.generic import get_standardized_pairwise_dists
from shnitsel.data.multi_indices import mdiff
from sklearn.decomposition import PCA as sk_PCA

from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline

from shnitsel.core.typedefs import AtXYZ
from shnitsel.filtering.structure_selection import StructureSelection


@needs(coords_or_vars={'atXYZ', 'astate'})
def pca_and_hops(
    frames: xr.Dataset,
    mean: bool,
    n_components: int = 2,
    return_pca_object: bool = False,
) -> tuple[xr.DataArray | tuple[xr.DataArray, sk_PCA], xr.DataArray]:
    """Get PCA points and info on which of them represent hops

    Parameters
    ----------
    frames
        A Dataset containing 'atXYZ' and 'astate' variables
    mean
        mean center data before pca if true
    n_components, optional
        The number of principle components to return, by default 2
    return_pca_object, optional
        Whether to return the scikit-learn `PCA` object as well as the
        transformed data, by default False

    Returns
    -------
    pca_res
        The PCA-reduced pairwise interatomic distances together with the sklearn.PCA object if `return_pca_object=True`.
    hops_pca_coords
        `pca_res` filtered by hops, to facilitate marking hops when plotting

    """
    pca_res = pairwise_dists_pca(
        frames['atXYZ'], n_components=n_components, return_pca_object=return_pca_object
    )
    mask = mdiff(frames['astate']) != 0
    if isinstance(pca_res, tuple):
        hops_pca_coords = pca_res[0][mask]
    else:
        hops_pca_coords = pca_res[mask]

    return pca_res, hops_pca_coords


@needs(dims={'atom'})
def pairwise_dists_pca(
    atXYZ: AtXYZ,
    mean: bool = False,
    n_components: int = 2,
    return_pca_object: bool = False,
    **kwargs,
) -> xr.DataArray | tuple[xr.DataArray, sk_PCA]:
    """PCA-reduced pairwise interatomic distances

    Parameters
    ----------
    atXYZ
        A DataArray containing the atomic positions;
        must have a dimension called 'atom'
    mean
        mean center data before pca if true
    n_components, optional
        The number of principle components to return, by default 2
    return_pca_object, optional
        Whether to return the scikit-learn `PCA` object as well as the
        transformed data, by default False


    Returns
    -------
        A DataArray with the same dimensions as `atXYZ`, except for the 'atom'
        dimension, which is replaced by a dimension 'PC' containing the principal
        components (by default 2)
        If return_pca_object=True, the result will be a tuple that holds the sklearn.PCA object as second entry
    """

    descr = get_standardized_pairwise_dists(atXYZ, mean=mean)
    res = pca(
        descr, 'atomcomb', n_components=n_components, return_pca_object=True, **kwargs
    )

    if not return_pca_object:
        assert not isinstance(res, tuple)  # typing
    # assert not isinstance(res, tuple)  # typing
    return res


def pca_with_features(
    atXYZ: AtXYZ,
    features: StructureSelection,
    n_components: int = 2,
    return_pca_object: bool = False,
    **kwargs,
) -> xr.DataArray | tuple[xr.DataArray, sk_PCA]:
    """PCA-reduced pairwise interatomic distances

    Parameters
    ----------
    atXYZ (AtXYZ)
        A DataArray containing the atomic positions;
        must have a dimension called 'atom'
    features (StructureSelection)
        The structural feature selection to perform the
    n_components, optional
        The number of principle components to return, by default 2
    return_pca_object, optional
        Whether to return the scikit-learn `PCA` object as well as the
        transformed data, by default False


    Returns
    -------
        A DataArray with the same dimensions as `atXYZ`, except for the 'atom'
        dimension, which is replaced by a dimension 'PC' containing the principal
        components (by default 2)
        If return_pca_object=True, the result will be a tuple that holds the sklearn.PCA object as second entry
    """
    # TODO: FIXME: Implement feature-dependent PCA
    res = (
        atXYZ.pipe(subtract_combinations, 'atom')
        .pipe(norm)
        .pipe(
            pca,
            'atomcomb',
            **{
                'n_components': n_components,
                'return_pca_object': return_pca_object,
                **kwargs,
            },
        )
    )
    if not return_pca_object:
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
        The name of the array-dimension to reduce (i.e. the axis along which different
        features lie)
    n_components, optional
        The number of principle components to return, by default 2
    return_pca_object, optional
        Whether to return the scikit-learn `PCA` object as well as the
        transformed data, by default False

    Returns
    -------
    pca_res
        A DataArray with the same dimensions as ``da``, except for the dimension
        indicated by `dim`, which is replaced by a dimension ``PC`` of size ``n_components``
        If DataArray accessors are active, the following members will be added to
        the accessor of the result:

            - ``pca_res.st.loadings``: The PCA loadings as a DataArray
            - ``pca_res.st.pca_object``: The scikit-learn pipeline used for PCA,
              including the ``MinMaxScaler``
            - ``pca_res_st.use_to_transform(other_da: xr.DataArray)``: A function which
              transforms its argument (other data) using the pipeline that has been
              fitted to the current data.

        (NB. The above assumes that the accessor name used is ``st``, the default)
    [pca_object]
        The trained PCA object produced by scikit-learn, if return_pca_object=True

    Examples:
    ---------
    >>> pca_results1 = data1.st.pca('features')
    >>> pca_results1.st.loadings  # See the loadings
    >>> pca_results2 = pca_results1.st.use_to_transform(data2)
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

    if return_pca_object:
        return (pca_res, pipeline[-1])
    else:
        return pca_res


# Alternative names
principal_component_analysis = pca
PCA = pca
