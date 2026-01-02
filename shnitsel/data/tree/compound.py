from dataclasses import dataclass
from typing import Any, Callable, Generic, Hashable, Mapping, Self, TypeVar
from .data_group import DataGroup, GroupInfo
from .data_leaf import DataLeaf

DataType = TypeVar("DataType")
ResType = TypeVar("ResType")


@dataclass
class CompoundInfo:
    """Class to hold identifying and auxiliary info of a compound type in ShnitselDB"""

    compound_name: str = "unknown"
    compound_smiles: str | None = None


class CompoundGroup(Generic[DataType], DataGroup[DataType]):
    """DataTree node to keep track of all data associated with a common compound within the datatree"""

    _compound_info: CompoundInfo

    def __init__(
        self,
        name: str | None = None,
        compound_info: CompoundInfo | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[
            Hashable,
            DataGroup[DataType] | DataLeaf[DataType],
        ]
        | None = None,
        level_name: str | None = None,
        attrs: Mapping[str, Any] | None = None,
    ):
        if name is None:
            if compound_info is not None:
                name = compound_info.compound_name
            elif group_info is not None:
                name = group_info.group_name

        super().__init__(
            name,
            group_info=group_info,
            attrs=attrs,
            level_name=level_name,
            children=children,
        )

        self._compound_info = (
            compound_info if compound_info is not None else CompoundInfo()
        )

    def construct_copy(self, **kwargs) -> Self:
        if 'compound_info' not in kwargs:
            kwargs['compound_info'] = self._compound_info
        return super().construct_copy(**kwargs)

    @property
    def compound_info(self) -> CompoundInfo:
        """Get the stored compound info of this Compound group.

        Returns:
            CompoundInfo: The metadata for the compound in this compound group
        """
        return self._compound_info

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "CompoundGroup[ResType]|None":
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
            return CompoundGroup[ResType](
                name=self._name,
                compound_info=self._compound_info,
                group_info=self._group_info,
                children=new_children,
                level_name=self._level_name,
                attrs=dict(self.attrs),
            )
