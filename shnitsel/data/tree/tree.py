from typing import Any, Callable, Generic, Hashable, Mapping, TypeVar

from shnitsel.data.tree.data_leaf import DataLeaf

from .data_group import DataGroup, GroupInfo
from .node import TreeNode

from .compound import CompoundGroup, CompoundInfo


DataType = TypeVar("DataType")
ResType = TypeVar("ResType")


class ShnitselDB(Generic[DataType], TreeNode[CompoundGroup[DataType], DataType]):
    def __init__(
        self, compounds: Mapping[Hashable, CompoundGroup[DataType]] | None = None
    ):
        super().__init__(
            name="ROOT",
            data=None,
            children=compounds or {},
        )

    def add_compound(
        self,
        name: str | None = None,
        compound_info: CompoundInfo | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[Hashable, DataGroup[DataType] | DataLeaf[DataType]]
        | None = None,
        attrs: Mapping[str, Any] | None = None,
    ) -> CompoundGroup[DataType]:
        # TODO: FIXME: Should add compound to tree
        new_compound = CompoundGroup[DataType](
            name=name,
            compound_info=compound_info,
            group_info=group_info,
            children=children,
            attrs=attrs,
        )
        return new_compound

    @property
    def compounds(self) -> Mapping[Hashable, CompoundGroup[DataType]]:
        return self.children

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "ShnitselDB[ResType]":
        new_children = None
        if recurse:
            new_children = {
                k: res
                for k, v in self._children.items()
                if v is not None
                and (res := v.map_data(func, recurse, keep_empty_branches)) is not None
            }

            if len(new_children) == 0 and not keep_empty_branches:
                new_children = None

        return ShnitselDB[ResType](new_children)
