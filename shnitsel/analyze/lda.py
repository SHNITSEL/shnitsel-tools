from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import MinMaxScaler
import xarray as xr

from shnitsel import _state


# TODO: FIXME: We have no example for how this is used, so documentation is purely guesswork.
def lda(
    dataset: xr.Dataset, dim: str, cats: str | xr.DataArray, n_components: int = 2
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

    if isinstance(cats, str):
        cats_name = cats
        cats = dataset[cats]
        # TODO: FIXME: I assume we need to delete the cats entry in the dataset here?
    else:
        cats_name = cats.name

    # Scale all features/variables in the dataset to be in 0.0-1.0 range
    scaled = xr.apply_ufunc(MinMaxScaler().fit_transform, dataset.transpose(..., dim))
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
        lda_object.scalings_,
        dims=(dim, 'scaling'),
        coords={'dim': dataset.coords[dim]},
    )

    if _state.DATAARRAY_ACCESSOR_REGISTERED:
        # TODO: FIXME: Add some explanation of this magic
        accessor_object = getattr(lda_res, _state.DATAARRAY_ACCESSOR_NAME)
        accessor_object.scalings = scalings
        accessor_object.lda_object = lda_object

    return lda_res


# Alternative names for the analysis function
linear_discriminat_analysis = lda
LDA = lda
