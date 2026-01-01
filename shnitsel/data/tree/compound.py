from dataclasses import dataclass
from typing import Any, Generic, Hashable, Mapping, TypeVar
from .data_group import DataGroup, GroupInfo
from .data_leaf import DataLeaf

DataType = TypeVar("DataType")


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
        dtype: type[DataType] | None = None,
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
            dtype=dtype,
        )

        self._compound_info = (
            compound_info if compound_info is not None else CompoundInfo()
        )

    @property
    def compound_info(self) -> CompoundInfo:
        """Get the stored compound info of this Compound group.

        Returns:
            CompoundInfo: The metadata for the compound in this compound group
        """
        return self._compound_info
