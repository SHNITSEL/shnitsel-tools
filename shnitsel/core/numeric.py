import itertools
import math
from typing import Hashable, TypeAlias

import numpy as np
import xarray as xr

from . import xrhelpers

DimName: TypeAlias = Hashable


def norm(
    da: xr.DataArray, dim: DimName = 'direction', keep_attrs: bool | str | None = None
) -> xr.DataArray:
    """Calculate the 2-norm of a DataArray, reducing the dimension with name *dim*

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
        return xrhelpers.midx_combs(da.get_index(dim))[f'{dim}comb']

    if dim not in da.dims:
        raise ValueError(f"'{dim}' is not a dimension of the DataArray")

    combination_dimension_name = f"{dim}comb"
    if combination_dimension_name in da:
        # TODO: FIXME: Appropriately remove the 'combination_dimension_name' indices and coordinates from the Dataset
        raise ValueError(
            f"'{combination_dimension_name}' is already an index, a variable or a coordinate of the DataArray"
        )

    n = da.sizes[dim]

    mat = np.zeros((math.comb(n, 2), n))
    combs = itertools.combinations(range(n), 2)

    # After matrix multiplication, index r of output vector has value c2 - c1
    for r, (c1, c2) in enumerate(combs):
        mat[r, c1] = -1
        mat[r, c2] = 1

    if labels:
        xrmat = xr.DataArray(
            data=mat, coords={f'{dim}comb': midx(da, dim), dim: da.get_index(dim)}
        )
    else:
        xrmat = xr.DataArray(
            data=mat,
            dims=[f'{dim}comb', dim],
        )

    newdims = list(da.dims)
    newdims[newdims.index(dim)] = f'{dim}comb'

    res = (xrmat @ da).transpose(*newdims)
    res.attrs = da.attrs
    res.attrs['deltaed'] = set(res.attrs.get('deltaed', [])).union({dim})
    return res