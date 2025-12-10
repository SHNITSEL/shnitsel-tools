from dataclasses import dataclass
from typing import Any, Hashable, Mapping, TypeVar, Generic


ChildType = TypeVar("ChildType")
DataType = TypeVar("DataType")


@dataclass
class TreeNode(Generic[ChildType, DataType]):
    name: str
    data: DataType | None = None
    children: Mapping[Hashable, ChildType] = dict()
    attrs: Mapping[str, Any] = dict()

    
@dataclass
class LeafNode(Generic[ChildType, DataType]):
    name: str
    data: DataType | None = None
    children: Mapping[Hashable, ChildType] = dict()
    attrs: Mapping[str, Any] = dict()
