from dataclasses import dataclass
from types import UnionType
from typing import Any, Callable, Generic, Hashable, Mapping, Self, TypeVar, overload
from typing_extensions import TypeForm

from shnitsel.data.tree.child_support_functions import find_child_key

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
        *,
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
        if compound_info is None:
            compound_info = CompoundInfo()

        if name is None:
            if compound_info is not None:
                name = compound_info.compound_name
            elif group_info is not None:
                name = group_info.group_name

        if level_name is None:
            level_name = DataTreeLevelMap['compound']

        super().__init__(
            name=name,
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
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | None = None,
        **kwargs,
    ) -> "CompoundGroup[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | UnionType | None = None,
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
        dtype: type[ResType] | UnionType | None = None,
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
                isinstance(child, (DataGroup, DataLeaf)) for child in children.values()
            ), (
                "Children provided to `construct_copy` for compound group are not of type `DataGroup` or `DataLeaf"
            )
            new_dtype: type[ResType] | UnionType | None = dtype

            return CompoundGroup(
                children=children,
                dtype=new_dtype,
                **kwargs,
            )

    @property
    def compound_info(self) -> CompoundInfo:
        """Get the stored compound info of this Compound group.

        Returns:
            CompoundInfo: The metadata for the compound in this compound group
        """
        return self._compound_info

    def add_data_group(
        self,
        group_info: GroupInfo,
        filter_func_data: Callable[[DataGroup | DataLeaf], bool] | None = None,
        flatten_data=False,
        **kwargs,
    ) -> Self:
        """Function to add trajectories within this compound subtree to a `TrajectoryGroup` of trajectories.

        The `group_name` will be set as the name of the group in the tree.
        If `flatten_trajectories=True` all existing groups will be dissolved before filtering and the children will be turned into an ungrouped list of trajectories.
        The `filter_func_trajectories` will either be applied to only the current groups and trajectories immediately beneath this compound or to the flattened list of all child directories.

        Args:
            group_name (str): The name to be set for the TrajectoryGroup object
            filter_func_Trajectories (Callable[[Trajectory|GroupInfo], bool] | None, optional): A function to return true for Groups and individual trajectories that should be added to the new group. Defaults to None.
            flatten_trajectories (bool, optional): A flag whether all descendant groups should be dissolved and flattened into a list of trajectories first before applying a group. Defaults to False.

        Returns:
            CompoundGroup: The restructured Compound with a new added group if at least one trajectory has satisfied the filter condition.
        """
        base_children: dict[Hashable, DataGroup[DataType] | DataLeaf[DataType]]
        # Either get the children or a map of all data nodes
        if flatten_data:
            all_data = self.collect_data_nodes()
            base_children = {}

            for x in all_data:
                x_key = find_child_key(base_children.keys(), x, "data")
                base_children[x_key] = x
        else:
            base_children = dict(self.children)

        # Filter out the nodes that should be in this group
        grouped_data: dict[Hashable, DataGroup[DataType] | DataLeaf[DataType]]
        ungrouped_data: dict[Hashable, DataGroup[DataType] | DataLeaf[DataType]]

        if filter_func_data is None:
            # All nodes should be in the group
            grouped_data = {
                str(k): v.construct_copy() for k, v in base_children.items()
            }
            ungrouped_data = {}
        else:
            # Split into data in new group and other data
            grouped_data = {}
            ungrouped_data = {}
            for key, child in base_children.items():
                if filter_func_data(child):
                    grouped_data[key] = child.construct_copy()
                else:
                    ungrouped_data[key] = child.construct_copy()

        # Build new group
        new_group = DataGroup(group_info=group_info, children=grouped_data)
        new_group_key = find_child_key(ungrouped_data.keys(), new_group, "group")
        ungrouped_data[new_group_key] = new_group
        return self.construct_copy(children=ungrouped_data)
