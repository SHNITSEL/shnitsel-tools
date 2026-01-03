from typing import Any, Callable, Generic, Hashable, Iterable, Mapping, Self, TypeVar

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.trajectory_grouping_params import TrajectoryGroupingMetadata
from shnitsel.data.tree.data_leaf import DataLeaf

from .data_group import DataGroup, GroupInfo
from .node import TreeNode

from .compound import CompoundGroup, CompoundInfo


DataType = TypeVar("DataType", covariant=True)
ResType = TypeVar("ResType")
KeyType = TypeVar("KeyType")


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

    def __init__(
        self, compounds: Mapping[Hashable, CompoundGroup[DataType]] | None = None
    ):
        super().__init__(
            name="ROOT",
            data=None,
            children=compounds or {},
        )

    def construct_copy(self, **kwargs) -> Self:
        """Helper function to create a copy of this node of the same type, but with potential changes to metadata, data or children

        Parameters:
        -----------
        **kwargs, optional
            Keyword arguments for the constructor of this type. If parameters for the constructor are not set, will be populated with the relevant values
            set in this instance.


        Returns:
        -----------
            Self: A copy of this node with recursively copied children if `children` is not set with an appropriate mapping.
        """
        return super().construct_copy(**kwargs)  # type: ignore # function is built such that the main type is preserved, only template arguments may change

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

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = True,
    ) -> "ShnitselDBRoot[ResType]":
        """Function to map the data in this subtree using a map function `func` and build a new tree
        with the same structure (if `keep_empty_branches=True`) but mapped data.

        Parameters
        ----------
        func : Callable[[DataType], ResType  |  None]
            The mapping function to be applied to data in this tree. The data is generally limited to the Leaf-level in our data structure.
        recurse : bool, optional
            Flag to map the function automatically across child levels of this subtree, by default True
        keep_empty_branches : bool, optional
            Flag to control whether empty branches should be omitted. If set to False, branches with no data and only empty children will be
            truncated. Only nodes with at least one child with mapped data will be kept. If True, the resulting tree will have the same
            structure as the input tree, by default True.

        Returns
        -------
        ShnitselDBRoot[ResType]
            A new tree with all data mapped by `func` to a new type and optionally empty branches truncated.
        """
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

    def map_flat_group_data(
        self, map_func: Callable[[Iterable[DataType]], ResType]
    ) -> "ShnitselDBRoot[ResType]":
        """Helper function to apply a mapping function to all flat group nodes.

        Will only apply the mapping function to nodes of type `DataGroup` and only those who have exclusively `DataLeaf` children.


        Args:
            map_func (Callable[[Iterable[DataType]], ResType]): Function mapping the data in the flat groups to a new result type

        Returns:
            ShnitselDBRoot[ResType]: A new tree structure, which will hold leaves with ResType data underneath each mapped group.
        """

        def extended_mapper(
            flat_group: TreeNode[Any, DataType],
        ) -> TreeNode[Any, ResType]:
            assert isinstance(flat_group, DataGroup)
            assert len(flat_group.subgroups) == 0
            child_data = {k: v.data for k, v in flat_group.subleaves.items()}
            # Actually perform the mapping over child data
            res = map_func(child_data.values())

            new_leaf = DataLeaf[ResType](name="reduced", data=res)
            new_group = flat_group.construct_copy(children={new_leaf.name: new_leaf})
            return new_group

        def filter_flat_groups(node: TreeNode[Any, DataType]) -> bool:
            return isinstance(node, DataGroup) and node.is_flat_group

        return self.map_filtered_nodes(
            filter_func=filter_flat_groups, map_func=extended_mapper
        )

    def map_filtered_nodes(
        self,
        filter_func: Callable[[TreeNode[Any, DataType]], bool],
        map_func: Callable[[TreeNode[Any, DataType]], TreeNode[Any, ResType]],
    ) -> "ShnitselDBRoot[ResType]":
        """_summary_

        _extended_summary_

        Parameters
        ----------
        filter_func : Callable[[TreeNode[Any, DataType]], bool]
            _description_
        map_func : Callable[[TreeNode[Any, DataType]], TreeNode[Any, ResType]]
            _description_

        Returns
        -------
        ShnitselDBRoot[ResType]
            _description_
        """
        return super().map_filtered_nodes(filter_func, map_func)  # type: ignore # The mapped node of a ShnitelDBRoot will be of the same type

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
        new_node = self.construct_copy()
        new_node._children = new_children
        return new_node

    def group_data_by_metadata(self) -> Self:
        """Helper function to allow for grouping of data within the tree by the metadata
        extracted from Trajectories.

        Should only be called on trees where `DataType=Trajectory` or `DataType=Frames` or subtypes thereof.
        Will fail due to an attribute error or yield an empty tree otherwise.

        Returns
        -------
        Self
            A tree where leaves are grouped to have similar metadata and only leaves with the same metadata are within the same gorup.
        """
        return self.group_children_by(
            key_func=_trajectory_key_func, group_leaves_only=True
        )


def _trajectory_key_func(node: TreeNode) -> None | str | TrajectoryGroupingMetadata:
    """Helper function to extract trajectory metadata of leaf nodes for trees with
    appropriate data types.

    If applied to other nodes may yield a `None` key or just their `name` attribute as a `str`.

    Parameters
    ----------
    node : TreeNode
        The node to extract the `TrajectoryGroupingMetadata` metadata from.
        See `Trajectory.get_grouping_metadata()` for creation of the meta data
        instance.

    Returns
    -------
    None | str | TrajectoryGroupingMetadata
        The key to use for the grouping of this node.
    """
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
