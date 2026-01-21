from dataclasses import asdict
import logging
from types import UnionType
from typing import (
    Any,
    Callable,
    Generic,
    Hashable,
    Iterable,
    Mapping,
    Self,
    TypeVar,
    overload,
)
from typing_extensions import TypeForm

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.trajectory_grouping_params import TrajectoryGroupingMetadata
from shnitsel.data.tree.data_leaf import DataLeaf
from .datatree_level import DataTreeLevelMap

from .data_group import DataGroup, GroupInfo
from .node import TreeNode

from .compound import CompoundGroup, CompoundInfo


DataType = TypeVar("DataType", covariant=True)
ResType = TypeVar("ResType")
KeyType = TypeVar("KeyType")
NewChildType = TypeVar("NewChildType", bound=TreeNode)


class ShnitselDBRoot(Generic[DataType], TreeNode[CompoundGroup[DataType], DataType]):
    """Class to use as a root for a ShnitselDB tree structure with specific Node types at different layer depths.

    Will always have `CompoundGroup` entries on the layer underneath the root.
    Will only have data in `DataLeaf` instances.
    Between leaf and compound nodes, there may be arbitrary `DataGroup` layers to allow for hiearchical structuring.

    Parameters
    ----------
    DataType: TypeVar
        A covariant template type parameter describing the kind of data that may be located in the leaves of this tree.
    TreeNode[CompoundGroup[DataType], DataType]
        The basic tree node type that this root node represents. Allows for sharing of functions between different levels of the tree.
    """

    def __new__(
        cls,
        *,
        compounds: Mapping[Hashable, CompoundGroup[DataType]] | None = None,
        **kwargs,
    ):
        if compounds is not None:
            kwargs['children'] = compounds

        kwargs['data'] = None

        return TreeNode.__new__(
            cls,
            **kwargs,
        )

    def __init__(
        self,
        *,
        compounds: Mapping[Hashable, CompoundGroup[DataType]] | None = None,
        **kwargs,
    ):
        if 'name' not in kwargs or kwargs['name'] is None:
            kwargs['name'] = "ROOT"
        if 'level_name' not in kwargs or kwargs['level_name'] is None:
            kwargs['level_name'] = DataTreeLevelMap['root']

        if compounds is not None:
            kwargs['children'] = compounds

        super().__init__(
            data=None,
            **kwargs,
        )

    # @overload
    # def construct_copy(
    #     self,
    #     children: Mapping[Hashable, CompoundGroup[DataType]] | None = None,
    #     dtype: None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> Self: ...

    # @overload
    # def construct_copy(
    #     self,
    #     children: Mapping[Hashable, CompoundGroup[ResType]],
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> "ShnitselDBRoot[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, CompoundGroup[DataType]] | None = None,
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
    ) -> "ShnitselDBRoot[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: None = None,
        **kwargs,
    ) -> "ShnitselDBRoot[ResType]": ...

    # def construct_copy_(
    #     self,
    #     children: Mapping[Hashable, CompoundGroup[ResType]] | None = None,
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> Self | "ShnitselDBRoot[ResType]":
    def construct_copy(
        self,
        children: Mapping[Hashable, CompoundGroup[DataType]]
        | Mapping[Hashable, NewChildType]
        | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: ResType | None = None,
        **kwargs,
    ) -> Self | "ShnitselDBRoot[ResType]":
        """Helper function to create a copy of this tree structure, but with potential changes to metadata, data or children

        Parameters:
        -----------
        children: Mapping[Hashable, CompoundGroup[DataType]] Mapping[Hashable, CompoundGroup[ResType]], optional
            The mapping of children with a potentially new `DataType`. If not provided, will be copied from the current node's child nodes.
        dtype: type[ResType] | UnionType, optional
            The data type of the data in the copy constructed tree.
        data: None, optional
            Data setting not supported on this type of node.

        Returns:
        -----------
            Self: A copy of this node with recursively copied children if `children` is not set with an appropriate mapping.
        """
        assert data is None, "No data must be set on a root node"

        if 'attrs' not in kwargs:
            kwargs['attrs'] = dict(self._attrs)
        if 'level_name' not in kwargs:
            kwargs['level_name'] = str(self._level_name)

        new_dtype: type[ResType] | UnionType | type[DataType] | None = dtype

        if children is not None:
            if 'compounds' in kwargs and kwargs["compounds"] is not None:
                raise KeyError(
                    "Provided both `compounds` and `children` argument to `construct_copy`"
                )
            kwargs['compounds'] = children

        if 'compounds' not in kwargs or kwargs["compounds"] is None:
            assert dtype is None, (
                "Cannot cast the data type of the tree without reassigning children/compounds of appropriate new type."
            )
            kwargs["compounds"] = {
                # TODO: FIXME: Figure out this typing issue
                k: v.construct_copy()
                for k, v in self._children.items()
                if v is not None
            }
            return type(self)(
                **kwargs,
            )
        else:
            compound_candidates = kwargs["compounds"]
            assert all(
                isinstance(child, CompoundGroup)
                for child in compound_candidates.values()
            ), (
                f"Children provided to `construct_copy` for tree root are not of type `CompoundGroup`: {compound_candidates}"
            )
            return ShnitselDBRoot[ResType](
                dtype=new_dtype,
                **kwargs,
            )

    def add_compound(
        self,
        name: str | None = None,
        compound_info: CompoundInfo | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[Hashable, DataGroup[DataType] | DataLeaf[DataType]]
        | None = None,
        attrs: Mapping[str, Any] | None = None,
    ) -> Self:
        """Helper function to add a new compound to this data structure without manually
        creating a `CompoundGroup` instance

        A compound is provided with a name used as an identifier for the compound and
        optionally a more in-depth `CompoundInfo` object.
        Due to compounds also being a `DataGroup`, group information can optionally be set.
        Similarly, children and attributes for the compound can be provided.

        Parameters
        ----------
        name : str | None, optional
            The compound identifier under which to register the compound, by default None, meaning it will be taken from `compound_info`.
            If no name can be extracted, a random name may be assigned.
        compound_info : CompoundInfo | None, optional
            Optional data structure to provide Compound meta data, by default None.
        group_info : GroupInfo | None, optional
            Optional data structure to set grouping information on the compound, by default None.
        children : Mapping[Hashable, DataGroup[DataType]  |  DataLeaf[DataType]] | None, optional
            Optionally a mapping of children (e.g. Trajectories) to use in the CompoundGroup creation, by default None
        attrs : Mapping[str, Any] | None, optional
            A mapping of keys to attribute values to set on the CompoundGroup, by default None

        Returns
        -------
        Self
            A new tree structure with the CompoundGroup inserted.
        """
        new_compound = CompoundGroup[DataType](
            name=name,
            compound_info=compound_info,
            group_info=group_info,
            children=children,
            attrs=attrs,
        )
        return self.add_child(name, new_compound)

    def add_data_group(
        self,
        group_info: GroupInfo,
        filter_func_compound: Callable[[CompoundInfo], bool] | None = None,
        filter_func_data: Callable[[DataLeaf | DataGroup], bool] | None = None,
        flatten_compound_data: bool = False,
        **kwargs,
    ) -> Self:
        """
        Function to add a group under the compound level for arbitrary compounds.
        The group is inserted at the top level underneath `CompoundGroup` nodes.

        `filter_func_compound` can be used to only generate the group for certain compounds.
        This parameter should be a function that only returns True if the group should be created underneath this comound.
        `filter_func_data` can be used to select only specific groups and leaves out of the children of a compound to be part of this group.
        `flatten_compound_data` can be set to `True` if existing groups within a compound are supposed to be dissolved (i.e. all data leaves gathered and put directly as children of the Compound)

        Parameters
        ----------
        group_info : GroupInfo
            The name and optionally additional metadata of the group to be created
        filter_func_compound : Callable[[CompoundInfo], bool] | None, optional
            Filter function that should return True if the group should be created for this compound, by default None, meaning all compounds will be filtered.
        filter_func_data : Callable[[DataLeaf | DataGroup], bool] | None, optional
            Filter function to determine whether a group or data leaf should be included in the new group, by default None
        flatten_compound_data : bool, optional
            Flag to determine whether all trajectories under selected compounds should be ungrouped before selecting for the new group, by default False

        Returns
        -------
        Self
             A resulting ShnitselDB structure with the grouping applied.
        """
        new_children = {}

        for child_key, child in self.children.items():
            if filter_func_compound is None or filter_func_compound(
                child.compound_info
            ):
                new_children[child_key] = child.add_data_group(
                    group_info,
                    filter_func_data,
                    flatten_compound_data,
                    **kwargs,
                )
            else:
                new_children[child_key] = child.construct_copy()
        return self.construct_copy(children=new_children)

    def set_compound_info(
        self, compound: str | CompoundInfo, overwrite_all: bool = False
    ) -> Self:
        """Function to set the compound information on either all unknown compounds (`overwrite_all=False`) or for all trajectories in the tree
        creating a new CompoundGroup holding all trajectories. (if `overwrite_all=True`).

        By default, the compound info will only be applied to trajectories with unknown compounds.
        If all compounds are merged or a compound info is assigned that is already in use, the concerned compound subtrees will be merged
        before the new `compound_info` is applied.

        Parameters
        ----------
        compound : str | CompoundInfo
            Either the compound name as a string or the compound information to apply to either the unknown compounds or all data in the tree.
        overwrite_all : bool, optional
            Flag to control whether the compound group of all data should be overwritten, by default False

        Returns
        -------
        Self
            The updated database
        """
        from .support_functions import tree_merge

        if isinstance(compound, str):
            compound_info = CompoundInfo(compound_name=compound)
        else:
            compound_info = compound

        if overwrite_all:
            new_compound: CompoundGroup[DataType] | None = tree_merge(
                *self.children.values(), res_data_type=self._dtype
            )
            if new_compound is not None:
                return self.construct_copy(
                    compounds={
                        compound_info.compound_name: new_compound.construct_copy(
                            compound_info=compound_info, dtype=self._dtype
                        )
                    }
                )
            else:
                # No compounds, so we can just return a copy
                return self.construct_copy()
        else:
            if "unknown" not in self.children:
                logging.warning(
                    "Non `unknown` compounds in tree to assign the compound info to."
                )
                return self.construct_copy()
            else:
                new_children = {
                    k: v.construct_copy()
                    for k, v in self.children.items()
                    if k != "unknown" and v is not None
                }
                unknown_child = self.children["unknown"]

                renamed_child = unknown_child.construct_copy(
                    name=compound_info.compound_name, compound_info=compound_info
                )
                res_name = renamed_child.name
                if res_name in new_children:
                    merged_child = tree_merge(renamed_child, new_children[res_name])
                    assert merged_child is not None, (
                        "Something went wrong with the merge of at least 2 trees."
                    )
                    new_children[res_name] = merged_child
                else:
                    new_children[res_name] = renamed_child

                # print(new_children)

                return self.construct_copy(compounds=new_children)

    @property
    def compounds(self) -> Mapping[Hashable, CompoundGroup[DataType]]:
        """The `compounds` held within this `ShnitselDB` structure.

        Auxiliary function to get the `children` property with a more domain-specific attribute name.

        Returns
        -------
        Mapping[Hashable, CompoundGroup[DataType]]
            The mapping of compound identifiers to the Compounds within this structure.
        """
        return self.children

    def group_children_by(
        self, key_func: Callable[[TreeNode], KeyType], group_leaves_only: bool = True
    ) -> Self:
        """This function creates a tree with likely a new structure having several desireable properties like groups
        either only having leaves or other groups underneath them and leaves within the same group having identical group keys.

        Specifically the grouping will generate a tree with the following properties:
        - CompoundGroup layer is left mostly untouched
        - DataGroup layers are refactored such that all leaves (or groups) within the same group have the same key resulting from `key_func`
        - If children with different `key_func` results are under the same group, a new group will be created to hold children with the same `key_func` result.
        - Nodes for which `key_func` yields `None` will not be retained.
        - if `group_leaves_only=True`, existing subgroups will be kept without invoking `key_func` and only leaves under the same group will be partitioned
          according to their `key_func` result.
        - If all children of an existing group yield the same `key` (NOTE: not `None`) result, then the group properties will be updated but the group will retain the same children.

        Parameters
        ----------
        key_func : Callable[[TreeNode], KeyType]
            A function to map all TreeNodes to a certain key that allows grouping by comparison and must be hashable. Ideally a dataclass result that allows the invocation of `as_dict()` to
            set group properties after grouping.
        group_leaves_only : bool, optional
            A flag whether grouping should only performed for `DataLeaf` type nodes, by default True.

        Returns
        -------
        Self
            A new tree with grouping performed across all `DataGroup` levels.
        """
        new_children = {
            k: res
            for k, v in self._children.items()
            if (res := v.group_children_by(key_func, group_leaves_only)) is not None
        }
        new_node = self.construct_copy(children=new_children)
        return new_node


ShnitselDB = ShnitselDBRoot
