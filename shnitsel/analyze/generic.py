import itertools
import logging
from typing import Collection, Hashable, Union

import numpy as np
import xarray as xr

from shnitsel.data.dataset_containers import Frames, Trajectory
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.multi_indices import midx_combs
from shnitsel.core.typedefs import AtXYZ
from shnitsel.filtering.structure_selection import StructureSelection

from ..core.typedefs import DataArrayOrVar, DimName


def norm(
    da: DataArrayOrVar, dim: DimName = 'direction', keep_attrs: bool | str | None = None
) -> DataArrayOrVar:
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
    res: DataArrayOrVar = xr.apply_ufunc(
        np.linalg.norm,
        da,
        input_core_dims=[[dim]],
        on_missing_core_dim='copy',
        kwargs={"axis": -1},
        keep_attrs=keep_attrs,
    )
    return res


def center(
    da: xr.DataArray, dim: DimName = 'frame', keep_attrs: Union[bool, str, None] = None
) -> xr.DataArray:
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
    with xr.set_options(keep_attrs=True):
        mean_da = da.mean(dim=dim)
        centered_da = da - mean_da
    return centered_da


# @needs(dims={'statecomb'}, coords={'statecomb'})
def subtract_combinations(
    da: xr.DataArray, dim: DimName, add_labels: bool = False
) -> xr.DataArray:
    """Calculate all possible pairwise differences over a given dimension

    Parameters
    ----------
    da : xr.DataArray
        Input DataArray; must contain dimension `dim` with an associated coordinate.
    dim : Dimname
        Dimension (of size $n$) to take pairwise differences over.
    add_labels : bool, optional
        If True, label the pairwise differences based on the index of `dim` otherwise do not add labels, by default False.

    Returns
    -------
        A DataArray with the dimension `dim` replaced by a dimension '`dim`comb' of size $n(n-1)/2$
    """

    def midx(da: xr.DataArray, dim: Hashable) -> xr.DataArray:
        """Small helper function to get a Coordinates object from the index of pairwise combinations of index values in `da` along
        dimension `dim`.

        Parameters
        ----------
        da : xr.DataArray
            The data for which to generate the combination index along dimension `dim`
        dim : Hashable
            Dimension along which to get the combinations.

        Returns
        -------
        xr.DataArray
            The DataArray object holding coordinate values for the `{dim}comb` dimension.
        """
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
        # TODO FIXME I don't understand this; if `da` already has a `{dim}comb`
        # dimension, then `xrmat` and `da` will have two dimensions in common
        # and the matrix multiplication will produce strange results or fail.
        # So if anything, shouldn't we raise an exception in that case?

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
        da = da.assign_coords()  # TODO FIXME What does this do?
        comb_indices = list(itertools.combinations(range(n), 2))
        coordinates = {combination_dimension_name: midx(da, dim), dim: dim_index}

    mat = np.zeros((len(comb_indices), n))

    # After matrix multiplication, index r of output vector has value c2 - c1
    for r, (c1, c2) in enumerate(comb_indices):
        mat[r, c1] = -1
        mat[r, c2] = 1

    if add_labels and coordinates is not None:
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
    """ "Function to calculate the norm of a variable across all dimensions except the ones denoted in `exclude`

    Used to obtain scalar representations of vector-values observables for plotting and other calculations
    where only the magnitude of the vector is of relevance.

    Parameters
    ----------
    da : xr.DataArray
        The data array to norm across all non-excluded dimensions
    exclude : Collection[DimName] | None, optional
        The dimensions to exclude/retain. Defaults to ['state', 'statecomb', 'frame', 'time'].

    Returns
    -------
    xr.DataArray
        The resulting, normed array

    Notes
    -----
        The output of keep_norming is not necessarily >= 0; for example, if all dimensions
        in ``da`` are in ``exclude``, the original object possibly containing negative values
        will be returned unaltered.
    """
    if exclude is None:
        exclude = {'state', 'statecomb', 'frame', 'time'}

    # Get all non-excluded dimensions
    # TODO: FIXME: This is not the same as calculating the norm across all non-excluded dimensions. This depends on the order of non-excluded dimensions!
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
    da : xr.DataArray
        An xr.DataArray to replace values within
    to_replace : np.ndarray | list
        Values to search for and replace -- these should be sortable, i.e. each pair of elements
        must be comparable by ``<``
    value : np.ndarray | list
        Values with which to replace the found occurrences of `to_replace` -- the dtype of this argument determines
        the dtype of the result

    Returns
    -------
        An xr.DataArray with dtype matching `value` obtained from `da` by replacing occurrences of `to_replace` with the corresponding values in `value`
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
    da :xr.DataArray
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


def pwdists(
    atXYZ_source: AtXYZ | xr.Dataset | ShnitselDataset, center_mean: bool = False
) -> xr.DataArray:
    """
    Compute pairwise distances and standardize it by removing the mean
    and L2-normalization (if your features are vectors and you want magnitudes only,
    to lose directional info)

    Parameters
    ----------
    atXYZ : xr.DataArray | AtXYZ | xr.Datset | ShnitselDataset
        A DataArray containing the atomic positions;
        Alternatively a dataset, a trajectory or a frameset with the positional data to derive the
        pairwise distances from.
    center_mean : bool, optional
        Subtract mean of calculated pairwise distances if `True`, by default `False`
    Returns
    -------
        A DataArray holding pairwise distance vectors with the same dimensions as
        but with a `descriptor` dimension indexing the pairwise distances.
    """
    from shnitsel.geo.geocalc import get_distances

    atxyz_ds: xr.Dataset
    atxyz_da: xr.DataArray
    if isinstance(atXYZ_source, xr.DataArray):
        atxyz_ds = atXYZ_source.to_dataset()
        atxyz_da = atXYZ_source
    else:
        atxyz_ds = atXYZ_source.dataset
        atxyz_da = atXYZ_source.positions

    struct_selection = StructureSelection.init_from_dataset(
        atxyz_ds, default_selection=['atoms']
    )

    bats_distances = get_distances(
        atxyz_da, structure_selection=struct_selection.select_pw_dists()
    )

    # res = atXYZ.pipe(subtract_combinations, 'atom', add_labels=True)

    # res = norm(res)
    if center_mean:
        bats_distances = center(bats_distances)

    return bats_distances


get_standardized_pairwise_dists = pwdists
