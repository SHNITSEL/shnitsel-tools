from contextlib import suppress
from dataclasses import dataclass
from typing import Any, Hashable, Iterable, Literal, Mapping, Self
import rdkit
import rdkit.Chem  # Avoid import error seen in RDKit 2025.09.3, Python 3.12.7
import xarray as xr

from ..xr_io_compatibility import (
    MetaData,
    ResType,
    SupportsFromXrConversion,
    SupportsToXrConversion,
)
from shnitsel.units.conversion import (
    convert_all_units_to_shnitsel_defaults,
    convert_to_target_units,
)


@dataclass
class ShnitselDataset(SupportsFromXrConversion, SupportsToXrConversion):
    _raw_dataset: xr.Dataset

    def __init__(self, ds: xr.Dataset):
        self._raw_dataset = ds

    @property
    def dataset(self) -> xr.Dataset:
        return self._raw_dataset

    @property
    def leading_dimension(self) -> str:
        if "frame" in self.dataset.dims:
            return "frame"
        elif "time" in self.dataset.dims:
            return "time"
        else:
            raise ValueError(
                "Unknown leading dimension of the contained dataset. The Dataset may have been misconstructed or loaded from malformed data."
            )

    @property
    def state_ids(self):
        if "state" not in self.dataset.coords:
            raise KeyError("No coordinate `state` provided for the trajectory")
        return self.dataset.coords["state"]

    @property
    def state_names(self):
        if "state_names" not in self.dataset.coords:
            raise KeyError("No coordinate `state_names` provided for the trajectory")
        return self.dataset.coords["state_names"]

    @property
    def state_types(self):
        if "state_types" not in self.dataset.coords:
            raise KeyError("No coordinate `state_types` provided for the trajectory")
        return self.dataset.coords["state_types"]

    @property
    def state_magnetic_number(self):
        if "state_magnetic_number" not in self.dataset.coords:
            raise KeyError(
                "No coordinate `state_magnetic_number` provided for the trajectory"
            )
        return self.dataset.coords["state_magnetic_number"]

    @property
    def state_degeneracy_group(self):
        if "state_degeneracy_group" not in self.dataset.coords:
            raise KeyError(
                "No coordinate `state_degeneracy_group` provided for the trajectory"
            )
        return self.dataset.coords["state_degeneracy_group"]

    @property
    def state_charges(self):
        if "state_charges" not in self.dataset.coords:
            raise KeyError("No coordinate `state_charges` provided for the trajectory")
        return self.dataset.coords["state_charges"]

    @property
    def active_state(self):
        if "astate" not in self.coords:
            if "astate" not in self.data_vars:
                raise KeyError(
                    "No coordinate `astate` holding the active state id provided for the trajectory"
                )
            return self.data_vars['astate']
        return self.coords["astate"]

    @property
    def state_diagonal(self):
        if "sdiag" not in self.dataset.coords:
            raise KeyError(
                "No coordinate `sdiag` holding the active state id provided for the trajectory"
            )
        return self.dataset.coords["sdiag"]

    @property
    def atom_names(self):
        if "atom_names" not in self.dataset.coords:
            raise KeyError("No coordinate `atom_names` provided for the trajectory")
        return self.dataset.coords["atom_names"]

    @property
    def atom_numbers(self):
        if "atom_numbers" not in self.dataset.coords:
            raise KeyError("No coordinate `atom_numbers` provided for the trajectory")
        return self.dataset.coords["atom_numbers"]

    @property
    def charge(self) -> float:
        """The charge of the molecule if set on the trajectory data.
        Loaded from `charge` attribute (or variable) or `state_charges` coordinate
        if provided.

        If no information is found, 0 is returned."""
        charge = self._param_from_vars_or_attrs('charge')
        if charge is None:
            if 'state_charges' in self._raw_dataset:
                return float(self._raw_dataset['state_charges'][0].item())
            return 0
        return float(charge)

    @charge.setter
    def charge(self, value):
        """Setter for the charge of a trajectory

        Parameters
        ----------
        float : float
            The charge in units of elementary charges.
        """
        # TODO: FIXME: We probably do not want to support setters here!
        self._raw_dataset.attrs['charge'] = value
        self._raw_dataset['state_charges'].values[:] = value

    def set_charge(self, value: float | xr.DataArray) -> Self:
        """Method to set the charge on a dataset, clear conflicting positions
        of charge info on the dataset and return a new instance of the wrapped dataset.


        Parameters
        ----------
        value : float | xr.DataArray
            Either a single value (optionally wrapped in a DataArray already) to indicate
            the charge of the full molecule in all states (will be set to coordinate `charge`) or a DataArray that represents
            state-dependent charges (which will be set to `state_charges`)

        Returns
        -------
        Self
            The updated object as a copy.

        Raises
        ------
        ValueError
            If an unsupported `value` was provided.
        """
        new_attrs = dict(self._raw_dataset.attrs)
        if 'charge' in new_attrs:
            del new_attrs['charge']

        tmp_ds = (
            self._raw_dataset.drop_attrs(deep=False)
            .assign_attrs(new_attrs)
            .drop_vars('charge', errors='ignore')
        )
        if isinstance(value, (float, int)):
            # Create a dummy, zero-dim
            charge_da = xr.DataArray(
                float(value),
                dims=[],
                attrs={'units': 'e', 'unitdim': 'charge'},
                name='charge',
            )

            new_ds = tmp_ds.assign_coords(charge=charge_da)
        elif isinstance(value, xr.DataArray):
            if len(value.sizes) == 0:
                # a scalar variable
                new_ds = tmp_ds.assign_coords(charge=value)
            else:
                # This is a state_charge contender:
                new_ds = tmp_ds.assign_coords(state_charges=value)
        else:
            raise ValueError("Unknown type for charge value: %s" % type(value))

        return type(self)(new_ds)

    # TODO: Forward all unmet requests to dataset.

    @property
    def dims(self):
        return self.dataset.dims

    @property
    def coords(self):
        return self.dataset.coords

    @property
    def sizes(self):
        return self.dataset.sizes

    @property
    def data_vars(self):
        return self.dataset.data_vars

    def has_variable(self, name: str) -> bool:
        return name in self.data_vars

    def has_dimension(self, name: str) -> bool:
        return name in self.dims

    def has_coordinate(self, name: str) -> bool:
        return name in self.coords

    def has_data(self, name: str) -> bool:
        return self.has_variable(name) or self.has_coordinate(name)

    def has(self, name: str) -> bool:
        return self.has_data(name) or self.has_dimension(name)

    @property
    def mol(self) -> rdkit.Chem.Mol:
        """Helper method to get a representative molecule object for the geometry within this dataset.

        Returns
        -------
        rdkit.Chem.Mol
            Either a copy of a cached mol object (for partial substructures) or a newly constructed default object
        """
        from shnitsel.bridges import construct_default_mol

        for key in ["__mol", "_mol", "mol"]:
            if key in self.attrs:
                return rdkit.Chem.Mol(self.attrs[key])

        mol_constr = construct_default_mol(self)
        # Cache the molecule
        self.attrs["__mol"] = mol_constr
        return mol_constr

    def sel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        method: str | None = None,
        tolerance: int | float | Iterable[int | float] | None = None,
        drop: bool = False,
        **indexers_kwargs: Any,
    ) -> "Self | ShnitselDataset":
        """Returns a new dataset with each data array indexed by tick labels
        along the specified dimension(s).

        In contrast to `.isel`, indexers for this method should use
        labels (i.e. explicit values in that dimension) instead of integers.

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
        dataset : Self
            A new Dataset with the same contents as this dataset, except each
            variable and dimension is indexed by the appropriate indexers.
            If indexer DataArrays have coordinates that do not conflict with
            this object, then these coordinates will be attached.
            In general, each array's data will be a view of the array's data
            in this dataset, unless vectorized indexing was triggered by using
            an array indexer, in which case the data will be a copy.

        See Also
        --------
        :func:`ShnitselDataset.isel <ShnitselDataset.isel>`
        :func:`Dataset.sel <Dataset.sel>`
        :func:`Dataset.isel <Dataset.isel>`
        :func:`DataArray.sel <DataArray.sel>`

        :doc:`xarray-tutorial:intermediate/indexing/indexing`
            Tutorial material on indexing with Xarray objects

        :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
            Tutorial material on basics of indexing

        """
        from shnitsel.data.dataset_containers import wrap_dataset

        # TODO: FIXME: Also select derived, base and cached variables?
        selres = self._raw_dataset.sel(
            indexers=indexers,
            method=method,
            tolerance=tolerance,
            drop=drop,
            **indexers_kwargs,
        )
        try:
            return type(self)(selres)
        except:
            return wrap_dataset(selres)  # type: ignore # Without a target type, datasets will be at least a ShnitselDataset

    def isel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        drop: bool = False,
        missing_dims: Literal["raise", "warn", "ignore"] = "raise",
        **indexers_kwargs: Any,
    ) -> "Self | ShnitselDataset":
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

        # A specific element from the dataset is selected

        >>> dataset.isel(atom=1, time=0)
        <xarray.Dataset> Size:
        Dimensions:         (direction: 3)
        Coordinates:
            atom        int16 2B 1
            time        float64 8B 0.0
            direction   (direction) <U1 3B 'x' 'y' 'z'
        Data variables:
            energy  float64 8B -238.2
            forces  (direction) float64 24B 1.2 -0.2 0.1

        # Indexing with a slice using isel

        >>> slice_of_data = dataset.isel(atom=slice(0, 2), time=slice(0, 2))
        >>> slice_of_data
        <xarray.Dataset> Size:
        Dimensions:         (atom: 2, time: 2, direction: 3)
        Coordinates:
            * atom         (atom) int16 2B 1
            * time         (time) float64 16B 0.0 0.5
            * direction    <U1 3B 'x' 'y' 'z'
        Data variables:
            energy      (time) float64 24B -238.2
            forces      (time, atom, direction) float64 96B -0.5 -0.4 0.4 ...

        >>> index_array = xr.DataArray([0, 2], dims="atom")
        >>> indexed_data = dataset.isel(atom=index_array)
        >>> indexed_data
        <xarray.Dataset> Size:
        Dimensions:         (atom: 2, time: 3, direction: 3)
        Coordinates:
          * atom            (atom) int16 4B 1 3
          * time            (time) float64 16B 0.0 0.5 1.0
          * direction       <U1 3B 'x' 'y' 'z'
        Data variables:
            energy      (time) float64 24B -238.2 -238.4 -237.9
            forces      (time, atom, direction) float64 96B -0.5 -0.4 0.4 ...

        See Also
        --------
        :func:`ShnitselDataset.sel <Dataset.sel>`
        :func:`Dataset.sel <Dataset.sel>`
        :func:`Dataset.isel <Dataset.sel>`
        :func:`DataArray.isel <DataArray.isel>`

        :doc:`xarray-tutorial:intermediate/indexing/indexing`
            Tutorial material on indexing with Xarray objects

        :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
            Tutorial material on basics of indexing

        """
        from shnitsel.data.dataset_containers import wrap_dataset
        # TODO: FIXME: Also select derived, base and cached variables?

        selres = self._raw_dataset.isel(
            indexers=indexers,
            drop=drop,
            missing_dims=missing_dims,
            **indexers_kwargs,
        )
        try:
            return type(self)(selres)
        except:
            return wrap_dataset(selres)  # type: ignore # Without a target type, the wrap result will be at least ShnitselDataset

    # @property
    # def _attr_sources(self) -> Iterable[Mapping[Hashable, Any]]:
    #     """Places to look-up items for attribute-style access"""
    #     yield from ()

    # @property
    # def _item_sources(self) -> Iterable[Mapping[Hashable, Any]]:
    #     """Places to look-up items for key-autocompletion"""
    #     yield from ()

    @property
    def _attr_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for attribute-style access"""
        yield from self._item_sources
        yield self._raw_dataset.attrs

    @property
    def _item_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for key-completion"""
        # print(hasattr(self, '_raw_dataset'))
        yield from self._raw_dataset._item_sources

    def __getattr__(self, name: str) -> Any:
        if name.startswith("__"):
            raise AttributeError(name)
        # print(name)
        if name not in {"__dict__", "__setstate__", "_raw_dataset"}:
            # this avoids an infinite loop when pickle looks for the
            # __setstate__ attribute before the object is initialized
            for source in self._attr_sources:
                with suppress(KeyError):
                    return source[name]

        if not name.startswith("_") and name in dir(self._raw_dataset):
            return getattr(self._raw_dataset, name)

        raise AttributeError(
            f"{type(self).__name__!r} object has no attribute {name!r}"
        )

    def __contains__(self, a):
        return self._raw_dataset.__contains__(a)

    def _repr_html_(self) -> Any:
        from .dataset_vis import shnitsel_dataset_repr

        return shnitsel_dataset_repr(self)

    def __getitem__(self, key):
        return self._raw_dataset[key]

    # This complicated two-method design boosts overall performance of simple operations
    # - particularly DataArray methods that perform a _to_temp_dataset() round-trip - by
    # a whopping 8% compared to a single method that checks hasattr(self, "__dict__") at
    # runtime before every single assignment. All of this is just temporary until the
    # FutureWarning can be changed into a hard crash.
    # def _setattr_dict(self, name: str, value: Any) -> None:
    #     """Deprecated third party subclass (see ``__init_subclass__`` above)"""
    #     object.__setattr__(self, name, value)
    #     if name in self.__dict__:
    #         # Custom, non-slotted attr, or improperly assigned variable?
    #         warnings.warn(
    #             f"Setting attribute {name!r} on a {type(self).__name__!r} object. Explicitly define __slots__ "
    #             "to suppress this warning for legitimate custom attributes and "
    #             "raise an error when attempting variables assignments.",
    #             FutureWarning,
    #             stacklevel=2,
    #         )

    # def __setattr__(self, name: str, value: Any) -> None:
    #     # """Objects with ``__slots__`` raise AttributeError if you try setting an
    #     # undeclared attribute. This is desirable, but the error message could use some
    #     # improvement.
    #     # """
    #     # try:
    #     #     object.__setattr__(self, name, value)
    #     # except AttributeError as e:
    #     #     # Don't accidentally shadow custom AttributeErrors, e.g.
    #     #     # DataArray.dims.setter
    #     #     if str(e) != f"{type(self).__name__!r} object has no attribute {name!r}":
    #     #         raise
    #     #     raise AttributeError(
    #     #         f"cannot set attribute {name!r} on a {type(self).__name__!r} object. Use __setitem__ style"
    #     #         "assignment (e.g., `ds['name'] = ...`) instead of assigning variables."
    #     #     ) from e
    #     # No support for arbitrary assignments currently
    #     pass

    def __dir__(self) -> list[str]:
        """Provide method name lookup and completion. Only provide 'public'
        methods.
        """
        extra_attrs = {
            item
            for source in self._attr_sources
            for item in source
            if isinstance(item, str)
        }
        ds_attrs = {item for item in dir(self._raw_dataset) if not item.startswith("_")}
        return sorted(set(dir(type(self))) | extra_attrs | ds_attrs)

    def _ipython_key_completions_(self) -> list[str]:
        """Provide method for the key-autocompletions in IPython.
        See https://ipython.readthedocs.io/en/stable/config/integrating.html#tab-completion
        For the details.
        """
        items = {
            item
            for source in self._item_sources
            for item in source
            if isinstance(item, str)
        }
        ds_completion = self._raw_dataset._ipython_key_completions_()
        items |= set(x for x in ds_completion if not x.startswith("_"))
        return list(items)

    def convert(self, varname: str | None = None, unit: str | None = None) -> Self:
        """Convert an entry in this dataset to a specific unit.

        Returns a copy of the dataset with the entry updated.

        Parameters
        ----------
        varname : str, optional
            Optionally the name of a single variable. If not provided, will apply to all variables.
        unit : str | None
            The target unit to convert to.
            If not set, Will convert to default shnitsel units.

        Returns
        -------
        Self
            The updated dataset with converted units.
        """
        if varname is None and unit is None:
            return type(self)(convert_all_units_to_shnitsel_defaults(self._raw_dataset))
        if varname is None and unit is not None:
            return type(self)(convert_to_target_units(self._raw_dataset, unit))
        if varname is not None:
            return type(self)(
                convert_to_target_units(self._raw_dataset, {varname: unit})
            )
        return self

    def as_xr_dataset(self) -> tuple[str | None, xr.Dataset, MetaData]:
        return self.get_type_marker(), self.dataset, dict()

    @classmethod
    def get_type_marker(cls) -> str:
        return "shnitsel::ShnitselDataset"

    @classmethod
    def from_xr_dataset(
        cls: type[ResType], dataset: xr.Dataset, metadata: MetaData
    ) -> ResType:
        return cls(dataset)

    # @overload
    # def __getitem__(self, key: Hashable) -> DataArray: ...

    # # Mapping is Iterable
    # @overload
    # def __getitem__(self, key: Iterable[Hashable]) -> Self: ...

    # def __getitem__(
    #     self, key: Mapping[Any, Any] | Hashable | Iterable[Hashable]
    # ) -> Self | DataArray:
    #     """Access variables or coordinates of this dataset as a
    #     :py:class:`~xarray.DataArray` or a subset of variables or a indexed dataset.

    #     Indexing with a list of names will return a new ``Dataset`` object.
    #     """
    #     from xarray.core.formatting import shorten_list_repr

    #     if utils.is_dict_like(key):
    #         return self.isel(**key)
    #     if utils.hashable(key):
    #         try:
    #             return self._construct_dataarray(key)
    #         except KeyError as e:
    #             message = f"No variable named {key!r}."

    #             best_guess = utils.did_you_mean(key, self.variables.keys())
    #             if best_guess:
    #                 message += f" {best_guess}"
    #             else:
    #                 message += f" Variables on the dataset include {shorten_list_repr(list(self.variables.keys()), max_items=10)}"

    #             # If someone attempts `ds['foo' , 'bar']` instead of `ds[['foo', 'bar']]`
    #             if isinstance(key, tuple):
    #                 message += f"\nHint: use a list to select multiple variables, for example `ds[{list(key)}]`"
    #             raise KeyError(message) from e

    #     if utils.iterable_of_hashable(key):
    #         return self._copy_listed(key)
    #     raise ValueError(f"Unsupported key-type {type(key)}")

    # def __setitem__(
    #     self, key: Hashable | Iterable[Hashable] | Mapping, value: Any
    # ) -> None:


@dataclass
class ShnitselDerivedDataset(
    ShnitselDataset, SupportsFromXrConversion, SupportsToXrConversion
):
    _base_dataset: xr.Dataset | None

    def __init__(self, base_ds: xr.Dataset | None, derived_ds: xr.Dataset):
        super().__init__(derived_ds)
        self._base_dataset = base_ds

    @property
    def base(self) -> xr.Dataset | None:
        return self.base

    @property
    def _item_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for key-completion"""
        # print(hasattr(self, '_raw_dataset'))
        yield from self._raw_dataset._item_sources
        if self._base_dataset is not None:
            yield from self._base_dataset._item_sources

    # TODO: Forward all unmet requests to dataset.

    def as_xr_dataset(self) -> tuple[str | None, xr.Dataset, MetaData]:
        raise NotImplementedError

    @classmethod
    def get_type_marker(cls) -> str:
        raise NotImplementedError

    @classmethod
    def from_xr_dataset(
        cls: type[ResType], dataset: xr.Dataset, metadata: MetaData
    ) -> ResType:
        raise NotImplementedError
