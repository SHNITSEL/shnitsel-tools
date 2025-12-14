import itertools
import logging
from typing import Collection, Union

import numpy as np
import xarray as xr

from shnitsel.data.multi_indices import midx_combs
from shnitsel.core.typedefs import AtXYZ

from ..core.typedefs import DimName


def norm(
    da: xr.DataArray, dim: DimName = 'direction', keep_attrs: bool | str | None = None
) -> xr.DataArray:
    """Calculate the 2-norm of a DataArray, reducing/squeezing the dimension with name `dim`

    Parameters
    ----------
    da
        Array to calculate the norm of
    dim, optional
        Dimension to calculate norm along (and therby reduce), by default 'direction'
    keep_attrs, optional
        How to deal with attributes; passed to xr.apply_ufunc, by default None

    Returns
    -------
        A DataArray with dimension *dim* reduced
    """
    res: xr.DataArray = xr.apply_ufunc(
        np.linalg.norm,
        da,
        input_core_dims=[[dim]],
        on_missing_core_dim='copy',
        kwargs={"axis": -1},
        keep_attrs=keep_attrs,
    )
    return res


def center(
    da: xr.DataArray,
    dim: DimName = 'frame',
    keep_attrs: Union[bool, str, None] = None) -> xr.DataArray:
    """
    Subtract the mean of a DataArray along a specified dimension.

    Parameters
    ----------
    da : DataArray
        Input array to be centered.
    dim : str, optional
        Dimension along which to compute the mean, by default 'frame'.
    keep_attrs : bool or str or None, optional
        How to handle attributes; passed to xr.apply_ufunc, by default None.

    Returns
    -------
    DataArray
        Centered DataArray with the same dimensions as input.
    """
    mean_da = da.mean(dim=dim)
    centered_da = xr.apply_ufunc(
        lambda x, m: x - m,
        da,
        mean_da,
        input_core_dims=[[], []],
        keep_attrs=keep_attrs,
        vectorize=True,
    )
    return centered_da


# @needs(dims={'statecomb'}, coords={'statecomb'})
def subtract_combinations(
    da: xr.DataArray, dim: DimName, labels: bool = False
) -> xr.DataArray:
    """Calculate all possible pairwise differences over a given dimension

    Parameters
    ----------
    da
        Input DataArray; must contain dimension `dim`
    dim
        Dimension (of size $n$) to take pairwise differences over
    labels, optional
        If True, label the pairwise differences based on the index of `dim`, by default False

    Returns
    -------
        A DataArray with the dimension `dim` replaced by a dimension '`dim`comb' of size $n(n-1)/2$
    """

    def midx(da, dim):
        return midx_combs(da.get_index(dim))[f'{dim}comb']

    if dim not in da.dims:
        raise ValueError(f"'{dim}' is not a dimension of the DataArray {da}")

    combination_dimension_name = f"{dim}comb"

    n = da.sizes[dim]
    dim_index = da.get_index(dim)

    coordinates = None
    dims = None
    dims = [combination_dimension_name, dim]

    if combination_dimension_name in da:
        # Don't recalculate the combinations, just take whichever have already been set.
        logging.info(
            f"Dimension {combination_dimension_name} already exists, reusing existing entries."
        )
        # Generate array indices from combination values
        comb_indices = []
        for c_from, c_to in da[combination_dimension_name].values:
            # TODO: Make sure that this is actually unique?
            index_from = dim_index.get_loc(c_from)
            index_to = dim_index.get_loc(c_to)
            comb_indices.append((index_from, index_to))
    else:
        logging.info(f"Dimension {combination_dimension_name} is being generated.")
        da = da.assign_coords()
        comb_indices = list(itertools.combinations(range(n), 2))
        coordinates = {combination_dimension_name: midx(da, dim), dim: dim_index}

    mat = np.zeros((len(comb_indices), n))

    # After matrix multiplication, index r of output vector has value c2 - c1
    for r, (c1, c2) in enumerate(comb_indices):
        mat[r, c1] = -1
        mat[r, c2] = 1

    if labels and coordinates is not None:
        xrmat = xr.DataArray(
            data=mat,
            coords=coordinates,
        )
    else:
        xrmat = xr.DataArray(data=mat, dims=dims)

    newdims = list(da.dims)
    newdims[newdims.index(dim)] = f'{dim}comb'

    res = (xrmat @ da).transpose(*newdims)
    res.attrs = da.attrs
    res.attrs['deltaed'] = set(res.attrs.get('deltaed', [])).union({dim})
    return res


def keep_norming(
    da: xr.DataArray, exclude: Collection[DimName] | None = None
) -> xr.DataArray:
    """Function to calculate the norm of a variable across all dimensions except the ones denoted in `exclude`

    Args:
        da (xr.DataArray): The data array to norm across all non-excluded dimensions
        exclude (Collection[DimName] | None, optional): The dimensions to exclude/retain. Defaults to ['state', 'statecomb', 'frame'].

    Returns:
        xr.DataArray: The resulting, normed array
    """
    if exclude is None:
        exclude = {'state', 'statecomb', 'frame'}

    # Get all non-excluded dimensions
    diff_dims = set(da.dims).difference(exclude)
    for dim in diff_dims:
        da = norm(da, dim, keep_attrs=True)
        da.attrs['norm_order'] = 2
    return da


def replace_total(
    da: xr.DataArray, to_replace: np.ndarray | list, value: np.ndarray | list
):
    """Replaces each occurence of `to_replace` in `da` with the corresponding element of `value`.
    Replacement must be total, i.e. every element of `da` must be in `to_replace`.
    This permits a change of dtype between `to_replace` and `value`.
    This function is based on the snippets at https://github.com/pydata/xarray/issues/6377

    Parameters
    ----------
    da
        An xr.DataArray
    to_replace
        Values to replace
    value
        Values with which to replace them

    Returns
    -------
        An xr.DataArray with dtype matching `value`.
    """
    to_replace = np.array(to_replace)
    value = np.array(value)
    flat = da.values.ravel()

    sorter = np.argsort(to_replace)
    insertion = np.searchsorted(to_replace, flat, sorter=sorter)
    indices = np.take(sorter, insertion, mode='clip')
    replaceable = to_replace[indices] == flat

    out = value[indices[replaceable]]
    return da.copy(data=out.reshape(da.shape))


def relativize(da: xr.DataArray, **sel) -> xr.DataArray:
    """Subtract the minimum of an xr.DataArray from all the array's elements

    Parameters
    ----------
    da
        The xr.DataArray from which to subtract the minimum
    **sel
        If keyword parameters are present, the reference minimum is picked
    from those elements that remain after running :py:meth:`xarray.DataArray.sel`
    using the keyword parameters as arguments.


    Returns
    -------
        The result of subtraction, with ``attrs`` intact.
    """
    res = da - da.sel(**sel).min()
    res.attrs = da.attrs
    return res


def get_standardized_pairwise_dists(
        atXYZ: AtXYZ, 
        mean: bool = True, 
        **kwargs) -> xr.DataArray:
    """
    Compute pairwise distances and standardize it by removing the mean 
    and L2-normalization (if your features are vectors and you want magnitudes only, 
    to lose directional info)

    Parameters
    ----------
    atXYZ
        A DataArray containing the atomic positions;
        must have a dimension called 'atom'
    mean : bool, optional
        If True, subtract the mean (centering), by default True

    Returns
    -------
        A DataArray with the same dimensions as `atXYZ` but transposed
    """

    res = atXYZ.pipe(subtract_combinations, 'atom')

    if mean:
        res = res.pipe(center)
    
    res = res.pipe(norm, **kwargs)

    assert not isinstance(res, tuple)
    return res




