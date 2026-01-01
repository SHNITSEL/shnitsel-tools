from typing import Any, Generic, Hashable, Mapping, TypeVar

from shnitsel.data.tree.data_leaf import DataLeaf

from ..dataset_containers.frames import Frames
from ..dataset_containers.trajectory import Trajectory
from .data_group import DataGroup, GroupInfo
from .node import TreeNode

from .compound import CompoundGroup, CompoundInfo


DataType = TypeVar("DataType")


class ShnitselDB(Generic[DataType], TreeNode[CompoundGroup[DataType], DataType]):
    def __init__(self, dtype=Trajectory | Frames):
        super().__init__(
            name="ROOT",
            data=None,
            children={},
            dtype=dtype,
        )

    def add_compound(
        self,
        name: str | None = None,
        compound_info: CompoundInfo | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[Hashable, DataGroup[DataType] | DataLeaf[DataType]]
        | None = None,
        attrs: Mapping[str, Any] | None = None,
        dtype=None,
    ) -> CompoundGroup[DataType]:
        if dtype is None:
            if self._dtype is not None:
                dtype = self._dtype
        else:
            if self._dtype is None:
                self._dtype = dtype
            else:
                if not issubclass(dtype, self._dtype):
                    raise ValueError(
                        "Cannot assign data type %s to a tree with dtype %s"
                        % (dtype, self._dtype)
                    )

        new_compound = CompoundGroup[DataType](
            name=name,
            compound_info=compound_info,
            group_info=group_info,
            children=children,
            attrs=attrs,
            dtype=dtype,
        )
        return new_compound

    @property
    def compounds(self) -> Mapping[Hashable, CompoundGroup[DataType]]:
        return self.children
