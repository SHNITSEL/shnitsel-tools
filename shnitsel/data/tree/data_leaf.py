from dataclasses import dataclass
from typing import Any, Generic, Mapping, TypeVar
from .node import TreeNode

DataType = TypeVar("DataType")


@dataclass
class DataLeaf(Generic[DataType], TreeNode[None, DataType]):
    def __init__(
        self,
        name: str | None = None,
        data: DataType | None = None,
        attrs: Mapping[str, Any] | None = None,
    ):
        super().__init__(
            name=name, data=data, attrs=attrs, level_name=self.__class__.__name__
        )
