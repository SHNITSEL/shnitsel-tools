from dataclasses import dataclass
from typing import Any, Callable, Generic, Mapping, TypeVar
from .node import TreeNode

DataType = TypeVar("DataType")
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
