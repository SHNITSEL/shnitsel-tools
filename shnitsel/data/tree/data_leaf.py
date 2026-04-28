from dataclasses import dataclass
from types import UnionType
from typing import (
    Any,
    Callable,
    Generic,
    Hashable,
    Literal,
    Mapping,
    Protocol,
    Self,
    TypeVar,
    overload,
)
from typing_extensions import TypeForm

from .datatree_level import DataTreeLevelMap
from .node import TreeNode

DataType = TypeVar("DataType", covariant=True)
ResType = TypeVar("ResType")
NewChildType = TypeVar("NewChildType", bound=TreeNode)


@dataclass
class DataLeaf(Generic[DataType], TreeNode[None, DataType]):
    """Class to represent a leaf node holding data in the ShnitselDB tree hierarchy.

    May be inherited from to provide leaves with more advanced features like provision
    of delayed results for support of parallel processing or delayed loading from disc, etc.
    """

    def __init__(
        self, *, name: str | None = None, data: DataType | None = None, **kwargs
    ):
        super().__init__(
            name=name, data=data, level_name=DataTreeLevelMap['data'], **kwargs
        )

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, None] | None = None,
        dtype: None = None,
        data: DataType | None = ...,
        **kwargs,
    ) -> Self: ...

    @overload
    def construct_copy(
        self,
        children: None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | None = ...,
        **kwargs,
    ) -> "DataLeaf[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: None = None,
        **kwargs,
    ) -> "DataLeaf[ResType]": ...

    # def construct_copy_(
    #     self,
    #     children: Mapping[Hashable, CompoundGroup[ResType]] | None = None,
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> Self | "ShnitselDBRoot[ResType]":

    def construct_copy(
        self,
        children: Mapping[Hashable, None]
        | Mapping[Hashable, NewChildType]
        | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | None = ...,
        **kwargs,
    ) -> Self | "DataLeaf[ResType]":
        """Helper function to create a copy of this tree structure, but with potential changes to metadata or data

        Parameters:
        -----------
        data: ResType | None, optional
            Data to replace the current data in the copy of this node
        children: None, optional
            Parameter not supported by this type of node.
        dtype: type[ResType] | TypeForm[ResType], optional
            The data type of the data in the copy constructed tree.

        Raises
        -----------
        AssertionError
            If dtype is set without a new `data` entry being provided

        Returns:
        -----------
            Self
                A copy of this node with recursively copied children if `data` is not set .
            DataLeaf[ResType]
                A new leaf with a new data type if `data` is provided.
        """
        assert children is None, "No children can be provided for a Leaf node"

        if 'name' not in kwargs:
            kwargs['name'] = self._name
        if 'attrs' not in kwargs:
            kwargs['attrs'] = self._attrs

        if data is ...:
            assert dtype is None, (
                "Cannot reassign data type if new data entry is not provided"
            )
            return type(self)(
                data=self._data,
                dtype=self._dtype,
                **kwargs,
            )
        else:
            return DataLeaf(
                data=data,
                dtype=dtype,
                **kwargs,
            )

    def group_children_by(
        self,
        key_func: Any = None,
        group_leaves_only: bool = False,
        recurse: bool = True,
    ) -> Self:
        """Specialization of the grouping operation for leaf nodes.

        Simply returns a copy of the current node.

        Parameters
        ----------
        key_func : Any, optional
            Unused, by default None
        group_leaves_only : bool, optional
            Unused, by default False
        recurse : bool, optional
            Unused, by default True

        Returns
        -------
        Self
            A copy of the current node. No further grouping possible at the leaf layer.
        """
        return self.construct_copy()


class ProvidesDelayedData(Generic[DataType], Protocol):
    """Helper class to encapsulate data that is not immediately accessible but may be provided by asynchronous loading or by
    computation in parallel/another process."""

    def get_data(self) -> DataType:
        """May be blocking if async processing is required"""
        ...

    def data_ready(self) -> bool:
        """Should not be blocking"""
        ...


@dataclass
class DelayedDataLeaf(DataLeaf[DataType]):
    """Class to hold data in a leaf of the tree structure,
    where the data is not immediately accessible but may be the delayed
    result of asynchronous processing.
    """

    _data_provider: ProvidesDelayedData[DataType] | None
    _data_is_set: bool

    def __init__(
        self,
        *,
        base_leaf: DataLeaf | None = None,
        name: str | None = None,
        data: ProvidesDelayedData[DataType] | None = None,
        attrs: Mapping[str, Any] | None = None,
        **kwargs,
    ):
        """Construct a delayed leaf either directly from the metadata and the data provider or from
        a base leaf with a data provider for the new data.

        Parameters
        ----------
        base_leaf : DataLeaf | None, optional
            Leaf to copy the metadata from, by default None
        name : str | None, optional
            Name of the new delayed leaf if not provided a `base_leaf`, by default None
        data : ProvidesDelayedData[DataType] | None, optional
            The provider of delayed data. Needs to implement the protocol `ProvidesDelayedData`, by default None
        attrs : Mapping[str, Any] | None, optional
            Attribute mapping of the new delayed leaf if not provided a `base_leaf`, by default None
        """
        # TODO: FIXME: The way that the type derivation of the tree works the delayed execution breaks the resulting tree dtype.
        data_dummy = None
        if base_leaf is not None:
            if 'name' not in kwargs:
                kwargs['name'] = base_leaf._name
            if 'attrs' not in kwargs:
                kwargs['attrs'] = base_leaf._attrs
            super().__init__(
                name=base_leaf.name, data=data_dummy, attrs=base_leaf.attrs, **kwargs
            )
        else:
            if 'name' not in kwargs:
                kwargs['name'] = self._name
            if 'attrs' not in kwargs:
                kwargs['attrs'] = self._attrs
            super().__init__(name=name, data=data_dummy, attrs=attrs, **kwargs)

        if data is not None:
            self._data_is_set = True
            self._data_provider = data
        else:
            self._data_is_set = False
            self._data_provider = None

    @property
    def has_data(self) -> bool:
        # TODO: FIXME: We may want to instead use the `._data_is_set` flag to check for data?
        try:
            return self.data is not None
        except ValueError:
            return False

    @property
    def data_ready(self) -> bool:
        """Property to check whether the data in this delayed node is ready to be accessed

        Returns
        -------
        bool
            _description_
        """
        try:
            # We have data already set
            if self.has_data:
                return True

            # Data is set to a provider but not yet retrieved.
            if self._data_is_set:
                if self._data_provider is not None and self._data_provider.data_ready():
                    return True
        except:
            return False

        return False

    @property
    def data(self) -> DataType:
        if not self._data_is_set:
            if self._data_provider is not None:
                new_data = self._data_provider.get_data()
                self._data = new_data
                self._data_is_set = True
            else:
                raise ValueError(
                    "Object has been created without any Data Provider. Data cannot be retrieved."
                )
        # Perform parent class data accession
        return super().data

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, None] | None = None,
        dtype: None = None,
        data: DataType | ProvidesDelayedData[ResType] | None = ...,
        **kwargs,
    ) -> Self: ...

    @overload
    def construct_copy(
        self,
        children: None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | ProvidesDelayedData[ResType] | None = ...,
        **kwargs,
    ) -> "DataLeaf[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: None = None,
        **kwargs,
    ) -> "DataLeaf[ResType]": ...

    # def construct_copy_(
    #     self,
    #     children: Mapping[Hashable, CompoundGroup[ResType]] | None = None,
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> Self | "ShnitselDBRoot[ResType]":

    def construct_copy(
        self,
        children: Mapping[Hashable, None]
        | Mapping[Hashable, NewChildType]
        | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | ProvidesDelayedData[ResType] | None = ...,
        **kwargs,
    ) -> Self | "DataLeaf[ResType]":
        """Helper function to create a copy of this tree structure, but with potential changes to metadata or data.

        If the data in this delayed DataLeaf is available, this will instead return a normal DataLeaf.

        Parameters:
        -----------
        data: ResType | ProvidesDelayedData[ResType] | None, optional
            Data to replace the current data in the copy of this node
        children: None, optional
            Parameter not supported by this type of node.
        dtype: type[ResType] | TypeForm[ResType], optional
            The data type of the data in the copy constructed tree.

        Raises
        -----------
        AssertionError
            If dtype is set without a new `data` entry being provided

        Returns:
        -----------
            Self
                A copy of this node with recursively copied children if `data` is not set .
            DataLeaf[ResType]
                A new leaf with a new data type if `data` is provided.
        """
        assert children is None, "No children can be provided for a Leaf node"

        if 'name' not in kwargs:
            kwargs['name'] = self._name
        if 'attrs' not in kwargs:
            kwargs['attrs'] = self._attrs

        if data is ...:
            assert dtype is None, (
                "Cannot reassign data type if new data entry is not provided"
            )
            return type(self)(
                data=self._data,
                dtype=self._dtype,
                **kwargs,
            )
        else:
            return DataLeaf(
                data=data,
                dtype=dtype,
                **kwargs,
            )
