from dataclasses import dataclass
from typing import Any, Callable, Generic, Hashable, Mapping, TypeVar
from .node import TreeNode
from .data_leaf import DataLeaf

DataType = TypeVar("DataType")
ResType = TypeVar("ResType")


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
    _group_info: GroupInfo | None = None

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
            name=name,
            data=None,
            children=children,
            attrs=attrs,
            level_name=level_name,
        )
        self._group_info = group_info

    def collect_data_nodes(self) -> list[DataLeaf[DataType]]:
        """Function to retrieve all nodes with data in this subtree

        Returns:
            list[DataLeaf[DataType]]: List of all nodes with DataLeaf Type in this tree.
        """
        res = []

        for x in self.children.values():
            if isinstance(x, DataGroup):
                res += x.collect_data_nodes()
            elif isinstance(x, DataLeaf):
                res.append(x)

        return res

    @property
    def group_info(self) -> GroupInfo:
        if self._group_info is None:
            raise ValueError("No group info set")
        else:
            return self._group_info

    @property
    def subgroups(self) -> Mapping[Hashable, "DataGroup[DataType]"]:
        from .data_group import DataGroup

        return {k: v for k, v in self._children.items() if isinstance(v, DataGroup)}

    @property
    def subleaves(self) -> Mapping[Hashable, "DataLeaf[DataType]"]:
        from .data_leaf import DataLeaf

        return {k: v for k, v in self._children.items() if isinstance(v, DataLeaf)}

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "DataGroup[ResType]|None":
        new_children: dict[Hashable, DataGroup[ResType] | DataLeaf[ResType]] | None = (
            None
        )
        if recurse:
            new_children = {
                k: res
                for k, v in self._children.items()
                if v is not None
                and (res := v.map_data(func, recurse, keep_empty_branches)) is not None
            }

            if len(new_children) == 0 and not keep_empty_branches:
                new_children = None

        if not keep_empty_branches and new_children is None:
            return None
        else:
            return DataGroup[ResType](
                name=self._name,
                group_info=self._group_info,
                children=new_children,
                level_name=self._level_name,
                attrs=dict(self.attrs),
            )
