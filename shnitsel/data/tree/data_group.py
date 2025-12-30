from dataclasses import dataclass
from typing import Any, Generic, Hashable, Mapping, TypeVar
from .node import TreeNode
from .data_leaf import DataLeaf

DataType = TypeVar("DataType")


@dataclass
class GroupInfo:
    """Class to hold auxiliary info of a group/collection of Data in ShnitselDB"""

    group_name: str
    group_attributes: dict[str, Any] | None = None
    grouped_properties: dict[str, float | str | int] | None = None


@dataclass
class DataGroup(
    Generic[DataType], TreeNode["DataGroup[DataType]|DataLeaf[DataType]", DataType]
):
    def __init__(
        self,
        name: str | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[
            Hashable,
            "DataGroup[DataType]|DataLeaf[DataType]",
        ]
        | None = None,
        attrs: Mapping[str, Any] | None = None,
        level_name: str | None = None,
    ):
        if name is None and group_info is not None:
            name = group_info.group_name

        super().__init__(
            name=name, data=None, children=children, attrs=attrs, level_name=level_name
        )
