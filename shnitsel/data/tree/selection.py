# WIP file


from typing import Any, Hashable, Mapping, Self
import numpy as np
from xarray import Index
from xarray.core.utils import (
    FilteredMapping,
    Frozen,
    FrozenMappingWarningOnValuesAccess,
    OrderedSet,
    drop_dims_from_indexers,
    either_dict_or_kwargs,
)


def is_fancy_indexer(indexer: Any) -> bool:
    """Return False if indexer is an int, slice, a 1-dimensional list, or a 0 or
    1-dimensional ndarray; in all other cases return True
    """
    if isinstance(indexer, int | slice) and not isinstance(indexer, bool):
        return False
    if isinstance(indexer, np.ndarray):
        return indexer.ndim > 1
    if isinstance(indexer, list):
        return bool(indexer) and not isinstance(indexer[0], int)
    return True


def isel_indexes(
    indexes: Indexes[Index],
    indexers: Mapping[Any, Any],
) -> tuple[dict[Hashable, Index], dict[Hashable, Variable]]:
    # Fast path function _apply_indexes_fast does not work with multi-coordinate
    # Xarray indexes (see https://github.com/pydata/xarray/issues/10063).
    # -> call it only in the most common case where all indexes are default
    # PandasIndex each associated to a single 1-dimensional coordinate.
    if any(type(idx) is not PandasIndex for idx in indexes._indexes.values()):
        return _apply_indexes(indexes, indexers, "isel")
    else:
        return _apply_indexes_fast(indexes, indexers, "isel")


def _apply_indexes(
    indexes: Indexes[Index],
    args: Mapping[Any, Any],
    func: str,
) -> tuple[dict[Hashable, Index], dict[Hashable, Variable]]:
    new_indexes: dict[Hashable, Index] = dict(indexes.items())
    new_index_variables: dict[Hashable, Variable] = {}

    for index, index_vars in indexes.group_by_index():
        index_dims = {d for var in index_vars.values() for d in var.dims}
        index_args = {k: v for k, v in args.items() if k in index_dims}
        if index_args:
            new_index = getattr(index, func)(index_args)
            if new_index is not None:
                new_indexes.update(dict.fromkeys(index_vars, new_index))
                new_index_vars = new_index.create_variables(index_vars)
                new_index_variables.update(new_index_vars)
            else:
                for k in index_vars:
                    new_indexes.pop(k, None)

    return new_indexes, new_index_variables


def _apply_indexes_fast(indexes: Indexes[Index], args: Mapping[Any, Any], func: str):
    # This function avoids the call to indexes.group_by_index
    # which is really slow when repeatedly iterating through
    # an array. However, it fails to return the correct ID for
    # multi-index arrays
    indexes_fast, coords = indexes._indexes, indexes._variables

    new_indexes: dict[Hashable, Index] = dict(indexes_fast.items())
    new_index_variables: dict[Hashable, Variable] = {}
    for name, index in indexes_fast.items():
        coord = coords[name]
        if hasattr(coord, "_indexes"):
            index_vars = {n: coords[n] for n in coord._indexes}
        else:
            index_vars = {name: coord}
        index_dims = {d for var in index_vars.values() for d in var.dims}
        index_args = {k: v for k, v in args.items() if k in index_dims}

        if index_args:
            new_index = getattr(index, func)(index_args)
            if new_index is not None:
                new_indexes.update(dict.fromkeys(index_vars, new_index))
                new_index_vars = new_index.create_variables(index_vars)
                new_index_variables.update(new_index_vars)
            else:
                for k in index_vars:
                    new_indexes.pop(k, None)
    return new_indexes, new_index_variables


def isel(
    self,
    indexers: Mapping[Any, Any] | None = None,
    drop: bool = False,
    missing_dims: ErrorOptionsWithWarn = "raise",
    **indexers_kwargs: Any,
) -> Self:
    """Returns a new dataset with each array indexed along the specified
    dimension(s).

    This method selects values from each array using its `__getitem__`
    method, except this method does not require knowing the order of
    each array's dimensions.

    Parameters
    ----------
    indexers : dict, optional
        A dict with keys matching dimensions and values given
        by integers, slice objects or arrays.
        indexer can be a integer, slice, array-like or DataArray.
        If DataArrays are passed as indexers, xarray-style indexing will be
        carried out. See :ref:`indexing` for the details.
        One of indexers or indexers_kwargs must be provided.
    drop : bool, default: False
        If ``drop=True``, drop coordinates variables indexed by integers
        instead of making them scalar.
    missing_dims : {"raise", "warn", "ignore"}, default: "raise"
        What to do if dimensions that should be selected from are not present in the
        Dataset:
        - "raise": raise an exception
        - "warn": raise a warning, and ignore the missing dimensions
        - "ignore": ignore the missing dimensions

    **indexers_kwargs : {dim: indexer, ...}, optional
        The keyword arguments form of ``indexers``.
        One of indexers or indexers_kwargs must be provided.

    Returns
    -------
    obj : Dataset
        A new Dataset with the same contents as this dataset, except each
        array and dimension is indexed by the appropriate indexers.
        If indexer DataArrays have coordinates that do not conflict with
        this object, then these coordinates will be attached.
        In general, each array's data will be a view of the array's data
        in this dataset, unless vectorized indexing was triggered by using
        an array indexer, in which case the data will be a copy.

    Examples
    --------

    >>> dataset = xr.Dataset(
    ...     {
    ...         "math_scores": (
    ...             ["student", "test"],
    ...             [[90, 85, 92], [78, 80, 85], [95, 92, 98]],
    ...         ),
    ...         "english_scores": (
    ...             ["student", "test"],
    ...             [[88, 90, 92], [75, 82, 79], [93, 96, 91]],
    ...         ),
    ...     },
    ...     coords={
    ...         "student": ["Alice", "Bob", "Charlie"],
    ...         "test": ["Test 1", "Test 2", "Test 3"],
    ...     },
    ... )

    # A specific element from the dataset is selected

    >>> dataset.isel(student=1, test=0)
    <xarray.Dataset> Size: 68B
    Dimensions:         ()
    Coordinates:
        student         <U7 28B 'Bob'
        test            <U6 24B 'Test 1'
    Data variables:
        math_scores     int64 8B 78
        english_scores  int64 8B 75

    # Indexing with a slice using isel

    >>> slice_of_data = dataset.isel(student=slice(0, 2), test=slice(0, 2))
    >>> slice_of_data
    <xarray.Dataset> Size: 168B
    Dimensions:         (student: 2, test: 2)
    Coordinates:
        * student         (student) <U7 56B 'Alice' 'Bob'
        * test            (test) <U6 48B 'Test 1' 'Test 2'
    Data variables:
        math_scores     (student, test) int64 32B 90 85 78 80
        english_scores  (student, test) int64 32B 88 90 75 82

    >>> index_array = xr.DataArray([0, 2], dims="student")
    >>> indexed_data = dataset.isel(student=index_array)
    >>> indexed_data
    <xarray.Dataset> Size: 224B
    Dimensions:         (student: 2, test: 3)
    Coordinates:
        * student         (student) <U7 56B 'Alice' 'Charlie'
        * test            (test) <U6 72B 'Test 1' 'Test 2' 'Test 3'
    Data variables:
        math_scores     (student, test) int64 48B 90 85 92 95 92 98
        english_scores  (student, test) int64 48B 88 90 92 93 96 91

    See Also
    --------
    :func:`Dataset.sel <Dataset.sel>`
    :func:`DataArray.isel <DataArray.isel>`

    :doc:`xarray-tutorial:intermediate/indexing/indexing`
        Tutorial material on indexing with Xarray objects

    :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
        Tutorial material on basics of indexing

    """
    indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "isel")
    if any(is_fancy_indexer(idx) for idx in indexers.values()):
        return self._isel_fancy(indexers, drop=drop, missing_dims=missing_dims)

    # Much faster algorithm for when all indexers are ints, slices, one-dimensional
    # lists, or zero or one-dimensional np.ndarray's
    indexers = drop_dims_from_indexers(indexers, self.dims, missing_dims)

    variables = {}
    dims: dict[Hashable, int] = {}
    coord_names = self._coord_names.copy()

    indexes, index_variables = isel_indexes(self.xindexes, indexers)

    for name, var in self._variables.items():
        # preserve variable order
        if name in index_variables:
            var = index_variables[name]
        else:
            var_indexers = {k: v for k, v in indexers.items() if k in var.dims}
            if var_indexers:
                var = var.isel(var_indexers)
                if drop and var.ndim == 0 and name in coord_names:
                    coord_names.remove(name)
                    continue
        variables[name] = var
        dims.update(zip(var.dims, var.shape, strict=True))

    return self._construct_direct(
        variables=variables,
        coord_names=coord_names,
        dims=dims,
        attrs=self._attrs,
        indexes=indexes,
        encoding=self._encoding,
        close=self._close,
    )


def _isel_fancy(
    self,
    indexers: Mapping[Any, Any],
    *,
    drop: bool,
    missing_dims: ErrorOptionsWithWarn = "raise",
) -> Self:
    valid_indexers = dict(self._validate_indexers(indexers, missing_dims))

    variables: dict[Hashable, Variable] = {}
    indexes, index_variables = isel_indexes(self.xindexes, valid_indexers)

    for name, var in self.variables.items():
        if name in index_variables:
            new_var = index_variables[name]
        else:
            var_indexers = {k: v for k, v in valid_indexers.items() if k in var.dims}
            if var_indexers:
                new_var = var.isel(indexers=var_indexers)
                # drop scalar coordinates
                # https://github.com/pydata/xarray/issues/6554
                if name in self.coords and drop and new_var.ndim == 0:
                    continue
            else:
                new_var = var.copy(deep=False)
            if name not in indexes:
                new_var = new_var.to_base_variable()
        variables[name] = new_var

    coord_names = self._coord_names & variables.keys()
    selected = self._replace_with_new_dims(variables, coord_names, indexes)

    # Extract coordinates from indexers
    coord_vars, new_indexes = selected._get_indexers_coords_and_indexes(indexers)
    variables.update(coord_vars)
    indexes.update(new_indexes)
    coord_names = self._coord_names & variables.keys() | coord_vars.keys()
    return self._replace_with_new_dims(variables, coord_names, indexes=indexes)


def map_index_queries(
    obj: T_Xarray,
    indexers: Mapping[Any, Any],
    method=None,
    tolerance: int | float | Iterable[int | float] | None = None,
    **indexers_kwargs: Any,
) -> IndexSelResult:
    """Execute index queries from a DataArray / Dataset and label-based indexers
    and return the (merged) query results.

    """
    from xarray.core.dataarray import DataArray

    # TODO benbovy - flexible indexes: remove when custom index options are available
    if method is None and tolerance is None:
        options = {}
    else:
        options = {"method": method, "tolerance": tolerance}

    indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "map_index_queries")
    grouped_indexers = group_indexers_by_index(obj, indexers, options)

    results = []
    for index, labels in grouped_indexers:
        if index is None:
            # forward dimension indexers with no index/coordinate
            results.append(IndexSelResult(labels))
        else:
            results.append(index.sel(labels, **options))

    merged = merge_sel_results(results)

    # drop dimension coordinates found in dimension indexers
    # (also drop multi-index if any)
    # (.sel() already ensures alignment)
    for k, v in merged.dim_indexers.items():
        if isinstance(v, DataArray):
            if k in v._indexes:
                v = v.reset_index(k)
            drop_coords = [name for name in v._coords if name in merged.dim_indexers]
            merged.dim_indexers[k] = v.drop_vars(drop_coords)

    return merged


def group_indexers_by_index(
    obj: T_Xarray,
    indexers: Mapping[Any, Any],
    options: Mapping[str, Any],
) -> list[tuple[Index, dict[Any, Any]]]:
    """Returns a list of unique indexes and their corresponding indexers."""
    unique_indexes = {}
    grouped_indexers: Mapping[int | None, dict] = defaultdict(dict)

    for key, label in indexers.items():
        index: Index = obj.xindexes.get(key, None)

        if index is not None:
            index_id = id(index)
            unique_indexes[index_id] = index
            grouped_indexers[index_id][key] = label
        elif key in obj.coords:
            raise KeyError(f"no index found for coordinate {key!r}")
        elif key not in obj.dims:
            raise KeyError(
                f"{key!r} is not a valid dimension or coordinate for "
                f"{obj.__class__.__name__} with dimensions {obj.dims!r}"
            )
        elif len(options):
            raise ValueError(
                f"cannot supply selection options {options!r} for dimension {key!r}"
                "that has no associated coordinate or index"
            )
        else:
            # key is a dimension without a "dimension-coordinate"
            # failback to location-based selection
            # TODO: depreciate this implicit behavior and suggest using isel instead?
            unique_indexes[None] = None
            grouped_indexers[None][key] = label

    return [(unique_indexes[k], grouped_indexers[k]) for k in unique_indexes]


def sel(
    self,
    indexers: Mapping[Any, Any] | None = None,
    method: str | None = None,
    tolerance: int | float | Iterable[int | float] | None = None,
    drop: bool = False,
    **indexers_kwargs: Any,
) -> Self:
    """Returns a new dataset with each array indexed by tick labels
    along the specified dimension(s).

    In contrast to `Dataset.isel`, indexers for this method should use
    labels instead of integers.

    Under the hood, this method is powered by using pandas's powerful Index
    objects. This makes label based indexing essentially just as fast as
    using integer indexing.

    It also means this method uses pandas's (well documented) logic for
    indexing. This means you can use string shortcuts for datetime indexes
    (e.g., '2000-01' to select all values in January 2000). It also means
    that slices are treated as inclusive of both the start and stop values,
    unlike normal Python indexing.

    Parameters
    ----------
    indexers : dict, optional
        A dict with keys matching dimensions and values given
        by scalars, slices or arrays of tick labels. For dimensions with
        multi-index, the indexer may also be a dict-like object with keys
        matching index level names.
        If DataArrays are passed as indexers, xarray-style indexing will be
        carried out. See :ref:`indexing` for the details.
        One of indexers or indexers_kwargs must be provided.
    method : {None, "nearest", "pad", "ffill", "backfill", "bfill"}, optional
        Method to use for inexact matches:

        * None (default): only exact matches
        * pad / ffill: propagate last valid index value forward
        * backfill / bfill: propagate next valid index value backward
        * nearest: use nearest valid index value
    tolerance : optional
        Maximum distance between original and new labels for inexact
        matches. The values of the index at the matching locations must
        satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
    drop : bool, optional
        If ``drop=True``, drop coordinates variables in `indexers` instead
        of making them scalar.
    **indexers_kwargs : {dim: indexer, ...}, optional
        The keyword arguments form of ``indexers``.
        One of indexers or indexers_kwargs must be provided.

    Returns
    -------
    obj : Dataset
        A new Dataset with the same contents as this dataset, except each
        variable and dimension is indexed by the appropriate indexers.
        If indexer DataArrays have coordinates that do not conflict with
        this object, then these coordinates will be attached.
        In general, each array's data will be a view of the array's data
        in this dataset, unless vectorized indexing was triggered by using
        an array indexer, in which case the data will be a copy.

    See Also
    --------
    :func:`Dataset.isel <Dataset.isel>`
    :func:`DataArray.sel <DataArray.sel>`

    :doc:`xarray-tutorial:intermediate/indexing/indexing`
        Tutorial material on indexing with Xarray objects

    :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
        Tutorial material on basics of indexing

    """
    indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "sel")
    query_results = map_index_queries(
        self, indexers=indexers, method=method, tolerance=tolerance
    )

    if drop:
        no_scalar_variables = {}
        for k, v in query_results.variables.items():
            if v.dims:
                no_scalar_variables[k] = v
            elif k in self._coord_names:
                query_results.drop_coords.append(k)
        query_results.variables = no_scalar_variables

    result = self.isel(indexers=query_results.dim_indexers, drop=drop)
    return result._overwrite_indexes(*query_results.as_tuple()[1:])
