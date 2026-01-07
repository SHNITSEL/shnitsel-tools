from dataclasses import dataclass
from typing import Any, Callable, Generic, Hashable, Mapping, Self, TypeVar, overload
from typing_extensions import TypeForm

from .datatree_level import DataTreeLevelMap
from .data_group import DataGroup, GroupInfo
from .data_leaf import DataLeaf
from .node import TreeNode

DataType = TypeVar("DataType", covariant=True)
ResType = TypeVar("ResType")
NewChildType = TypeVar("NewChildType", bound=TreeNode)


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
        **kwargs,
    ):
        if name is None:
            if compound_info is not None:
                name = compound_info.compound_name
            elif group_info is not None:
                name = group_info.group_name

        if level_name is None:
            level_name = DataTreeLevelMap['compound']

        super().__init__(
            name,
            group_info=group_info,
            attrs=attrs,
            level_name=level_name,
            children=children,
            **kwargs,
        )

        self._compound_info = (
            compound_info if compound_info is not None else CompoundInfo()
        )

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, "DataGroup[DataType]|DataLeaf[DataType]"]
        | None = None,
        dtype: None = None,
        data: DataType | None = None,
        **kwargs,
    ) -> Self: ...

    @overload
    def construct_copy(
        self,
        children: None = None,
        dtype: type[ResType] | TypeForm[ResType] | None = None,
        data: ResType | None = None,
        **kwargs,
    ) -> "CompoundGroup[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | TypeForm[ResType] | None = None,
        data: None = None,
        **kwargs,
    ) -> "CompoundGroup[ResType]": ...

    # def construct_copy_(
    #     self,
    #     children: Mapping[Hashable, CompoundGroup[ResType]] | None = None,
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> Self | "ShnitselDBRoot[ResType]":

    def construct_copy(
        self,
        children: Mapping[Hashable, "DataGroup[DataType]|DataLeaf[DataType]"]
        | Mapping[Hashable, NewChildType]
        | None = None,
        dtype: type[ResType] | TypeForm[ResType] | None = None,
        data: ResType | None = None,
        **kwargs,
    ) -> Self | "CompoundGroup[ResType]":
        """Helper function to create a copy of this tree structure, but with potential changes to metadata, data or children

        Parameters:
        -----------
        data: None, optional
            Data setting not supported on this type of node.
        children: Mapping[Hashable, CompoundGroup[ResType]], optional
            The mapping of children with a potentially new `DataType`. If not provided, will be copied from the current node's child nodes.
        dtype: type[ResType] | TypeForm[ResType], optional
            The data type of the data in the copy constructed tree.

        Raises
        -----------
        AssertionError
            If dtype is provided but children parameter not set and node has children, indicating an issue with a type update without setting the new children

        Returns:
        -----------
            Self: A copy of this node with recursively copied children if `children` is not set with an appropriate mapping.
        """
        assert data is None, "No data must be set on a root node"

        if 'name' not in kwargs:
            kwargs['name'] = self._name
        if 'compound_info' not in kwargs:
            kwargs['compound_info'] = self._compound_info
        if 'group_info' not in kwargs:
            kwargs['group_info'] = self._group_info

        if 'attrs' not in kwargs:
            kwargs['attrs'] = dict(self._attrs)
        if 'level_name' not in kwargs:
            kwargs['level_name'] = str(self._level_name)

        if children is None:
            return type(self)(
                children={
                    # TODO: FIXME: Figure out this typing issue
                    k: v.construct_copy()
                    for k, v in self._children.items()
                    if v is not None
                },
                dtype=self._dtype,
                **kwargs,
            )
        else:
            assert all(
                isinstance(child, (DataGroup, DataLeaf)) for child in children
            ), (
                "Children provided to `construct_copy` for compound group are not of type `DataGroup` or `DataLeaf"
            )
            new_dtype: type[ResType] | TypeForm[ResType] | None = dtype

            return CompoundGroup[ResType](
                children=children,
                dtype=new_dtype,
                **kwargs,
            )

    # def construct_copy(
    #     self,
    #     data: DataType | ResType | None = None,
    #     dtype: type[DataType]
    #     | TypeForm[DataType]
    #     | type[ResType]
    #     | TypeForm[ResType]
    #     | None = None,
    #     **kwargs,
    # ) -> Self | "CompoundGroup[ResType]":
    #     """Function to generate a copy of a `CompoundGroup` instance with either a new data type
    #     or the same datatype but potentially new children.

    #     Parameters
    #     ----------
    #     data: None
    #         No data must be set on a `CompoundGroup` node. Do not set this parameter
    #     dtype: type[DataType] | TypeForm[DataType] | type[ResType] | TypeForm[ResType], optional
    #         The data type of the data in the copy constructed group.
    #     **kwargs, optional
    #         Keyword arguments for the constructor of this type. If parameters for the constructor are not set, will be populated with the relevant values
    #         set in this instance.

    #     Returns
    #     -------
    #     CompoundGroup[DataType]
    #         If no new type was provided and the children were also not set with a new dtype.
    #     CompoundGroup[ResType]
    #         If a new dtype was set, the duplicate will have a new `DataType` set.
    #     """
    #     assert data is None, "No data must be set on a compound group node"
    #     if 'compound_info' not in kwargs:
    #         kwargs['compound_info'] = self._compound_info

    #     return super().construct_copy(dtype=dtype, **kwargs)

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
        dtype: type[ResType] | TypeForm[ResType] | None = None,
    ) -> "CompoundGroup[ResType]|None":
        """Function to map the data in this subtree using a map function `func` and build a new tree
        with the same structure (if `keep_empty_branches=True`) but mapped data.

        This is the Specialization for nodes of type `CompoundGroup`

        Parameters
        ----------
        func : Callable[[DataType], ResType  |  None]
            The mapping function to be applied to data in this tree. The data is generally limited to the Leaf-level in our data structure.
        recurse : bool, optional
            Flag to map the function automatically across child levels of this subtree, by default True
        keep_empty_branches : bool, optional
            Flag to control whether empty branches should be omitted. If set to False, branches with no data and only empty children will be
            truncated. Only nodes with at least one child with mapped data will be kept. If True, the resulting tree will have the same
            structure as the input tree, by default False
        dtype : type[ResType] | TypeForm[ResType] | None, optional
            Optional specific type parameter for the result of this mapping, by default None

        Returns
        -------
        DataGroup[ResType]|None
            A new tree with all data mapped by `func` to a new type and optionally empty branches truncated.
        """
        new_children: dict[Hashable, DataGroup[ResType] | DataLeaf[ResType]] | None = (
            None
        )
        if recurse:
            new_children = {
                k: res
                for k, v in self._children.items()
                if v is not None
                and (res := v.map_data(func, recurse, keep_empty_branches, dtype=dtype))
                is not None
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
                dtype=dtype,
            )
