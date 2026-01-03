from dataclasses import dataclass
from typing import Any, Callable, Generic, Mapping, Protocol, Self, TypeVar
from .node import TreeNode

DataType = TypeVar("DataType", covariant=True)
ResType = TypeVar("ResType")


@dataclass
class DataLeaf(Generic[DataType], TreeNode[None, DataType]):
    def __init__(
        self,
        name: str | None = None,
        data: DataType | None = None,
        attrs: Mapping[str, Any] | None = None,
    ):
        super().__init__(
            name=name,
            data=data,
            attrs=attrs,
            level_name=self.__class__.__name__,
        )

    def construct_copy(self, **kwargs) -> Self:
        return super().construct_copy(**kwargs)

    def group_children_by(
        self,
        key_func: Any = None,
        group_leaves_only: bool = False,
        recurse: bool = True,
    ) -> Self:
        return self.construct_copy()

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "DataLeaf[ResType] | None":
        if self.has_data:
            new_data = func(self.data)
        else:
            new_data = None

        if not keep_empty_branches and new_data is None:
            return None
        else:
            # This yields a different kind of tree.
            return DataLeaf[ResType](
                name=self.name, data=new_data, attrs=dict(self.attrs)
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
class DelayedDataLeaf(Generic[DataType], DataLeaf[DataType]):
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
