from dataclasses import dataclass
from typing import (
    Any,
    Callable,
    Generic,
    Hashable,
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

    def __init__(self, name: str | None = None, data: DataType | None = None, **kwargs):
        super().__init__(
            name=name, data=data, level_name=DataTreeLevelMap['data'], **kwargs
        )

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, None] | None = None,
        dtype: None = None,
        data: DataType | None = None,
        **kwargs,
    ) -> Self: ...

    @overload
    def construct_copy(
        self,
        children: None = None,
        dtype: type[ResType] | TypeForm[ResType] | None = None,
        data: ResType | None = None,
        **kwargs,
    ) -> "DataLeaf[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | TypeForm[ResType] | None = None,
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
        dtype: type[ResType] | TypeForm[ResType] | None = None,
        data: ResType | None = None,
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

        if data is None:
            assert dtype is None, (
                "Cannot reassign data type if new data entry is not provided"
            )
            return type(self)(
                name=self._name,
                data=self._data,
                dtype=self._dtype,
                **kwargs,
            )
        else:
            return DataLeaf[ResType](
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

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
        dtype: type[ResType] | TypeForm[ResType] | None = None,
    ) -> "DataLeaf[ResType] | None":
        """Specialization of the `map_data()` method for leaf nodes.

        Applies the `func()` function to available data in this node.

        Parameters
        ----------
        func : Callable[[DataType], ResType  |  None]
            The function to apply to the data in this node
        recurse : bool, optional
            Ignored by this level, by default True
        keep_empty_branches : bool, optional
            If set to True and the result of the mapping is `None` for this node, it will be dropped and `None` returned, by default False
        dtype : type[ResType] | TypeForm[ResType] | None, optional
            Optional specific type parameter for the result of this mapping, by default None

        Returns
        -------
        DataLeaf[ResType] | None
            The new DataLeaf with the updated data or None if no data could be obtained from the mapping and `keep_empty_branches=True`.
        """
        if self.has_data:
            new_data = func(self.data)
        else:
            new_data = None

        if not keep_empty_branches and new_data is None:
            return None
        else:
            # This yields a different kind of tree.
            return DataLeaf[ResType](
                name=self.name, data=new_data, attrs=dict(self.attrs), dtype=dtype
            )


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
        name: str | None = None,
        data_provider: ProvidesDelayedData[DataType] | None = None,
        attrs: Mapping[str, Any] | None = None,
    ):
        data_dummy = None
        super().__init__(name, data=data_dummy, attrs=attrs)

        if data_provider is not None:
            self._data_is_set = True
            self._data_provider = None
        else:
            self._data_is_set = False
            self._data_provider = data_provider

    @property
    def has_data(self) -> bool:
        try:
            return self.data is not None
        except ValueError:
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
