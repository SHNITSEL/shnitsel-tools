from typing import Any, Callable, Generic, Hashable, Mapping, Self, TypeVar

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.trajectory_grouping_params import TrajectoryGroupingMetadata
from shnitsel.data.tree.data_leaf import DataLeaf

from .data_group import DataGroup, GroupInfo
from .node import TreeNode

from .compound import CompoundGroup, CompoundInfo


DataType = TypeVar("DataType")
ResType = TypeVar("ResType")
KeyType = TypeVar("KeyType")


class ShnitselDBRoot(Generic[DataType], TreeNode[CompoundGroup[DataType], DataType]):
    def __init__(
        self, compounds: Mapping[Hashable, CompoundGroup[DataType]] | None = None
    ):
        super().__init__(
            name="ROOT",
            data=None,
            children=compounds or {},
        )

    def construct_copy(self, **kwargs) -> Self:
        return super().construct_copy(**kwargs)

    def add_compound(
        self,
        name: str | None = None,
        compound_info: CompoundInfo | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[Hashable, DataGroup[DataType] | DataLeaf[DataType]]
        | None = None,
        attrs: Mapping[str, Any] | None = None,
    ) -> Self:
        new_compound = CompoundGroup[DataType](
            name=name,
            compound_info=compound_info,
            group_info=group_info,
            children=children,
            attrs=attrs,
        )
        return self.add_child(name, new_compound)

    @property
    def compounds(self) -> Mapping[Hashable, CompoundGroup[DataType]]:
        return self.children

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "ShnitselDBRoot[ResType]":
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

        return ShnitselDBRoot[ResType](new_children)

    def group_children_by(
        self, key_func: Callable[[TreeNode], KeyType], group_leaves_only: bool = False
    ) -> Self:
        new_children = {
            k: res
            for k, v in self._children.items()
            if (res := v.group_children_by(key_func, group_leaves_only)) is not None
        }
        new_node = self.construct_copy()
        new_node._children = new_children
        return new_node

    def group_data_by_metadata(self) -> Self:
        return self.group_children_by(
            key_func=_trajectory_key_func, group_leaves_only=True
        )


def _trajectory_key_func(node: TreeNode) -> None | str | TrajectoryGroupingMetadata:
    if isinstance(node, DataLeaf):
        if not node.has_data:
            # Do not group empty data
            return None
        else:
            # Get grouping metadata
            if isinstance(node.data, Trajectory) or isinstance(node.data, Frames):
                return node.data.get_grouping_metadata()
        # Don't attempt to group weird data types
        return None
    else:
        return node.name
