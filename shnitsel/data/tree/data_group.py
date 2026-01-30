from dataclasses import dataclass, asdict
from types import UnionType
from typing import Any, Callable, Generic, Hashable, Mapping, Self, TypeVar, overload
from typing_extensions import TypeForm

from .datatree_level import DataTreeLevelMap
from .node import TreeNode
from .data_leaf import DataLeaf

DataType = TypeVar("DataType", covariant=True)
ResType = TypeVar("ResType")
KeyType = TypeVar("KeyType")
NewChildType = TypeVar("NewChildType", bound=TreeNode)


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
        *,
        name: str | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[
            Hashable,
            "DataGroup[DataType]|DataLeaf[DataType]",
        ]
        | None = None,
        attrs: Mapping[str, Any] | None = None,
        level_name: str | None = None,
        **kwargs,
    ):
        if name is None and group_info is not None:
            name = group_info.group_name

        if level_name is None:
            level_name = DataTreeLevelMap['group']

        super().__init__(
            name=name,
            data=None,
            children=children,
            attrs=attrs,
            level_name=level_name,
            **kwargs,
        )
        self._group_info = group_info

    # @overload
    # def construct_copy(
    #     self,
    #     children: None = None,
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> Self: ...

    # @overload
    # def construct_copy(
    #     self,
    #     children: Mapping[Hashable, "DataGroup[ResType] | DataLeaf[ResType]"],
    #     dtype: type[ResType] | TypeForm[ResType] | None = None,
    #     data: None = None,
    #     **kwargs,
    # ) -> "DataGroup[ResType]": ...

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
    ) -> "DataGroup[ResType]": ...

    @overload
    def construct_copy(
        self,
        children: Mapping[Hashable, NewChildType] | None = None,
        dtype: type[ResType] | UnionType | None = None,
        data: None = None,
        **kwargs,
    ) -> "DataGroup[ResType]": ...

    # def construct_copy_(
    #     self,
    #     children: Mapping[Hashable, CompoundGroup[ResType]] | None = None,
    #     dtype: type[ResType] | UnionType | None = None,
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
    ) -> Self | "DataGroup[ResType]":
        """Helper function to create a copy of this tree structure, but with potential changes to metadata, data or children

        Parameters:
        -----------
        data: None, optional
            Data setting not supported on this type of node.
        children: Mapping[Hashable, DataGroup[ResType]], optional
            The mapping of children with a potentially new `DataType`. If not provided, will be copied from the current node's child nodes.
        dtype: type[ResType] | TypeForm[ResType], optional
            The data type of the data in the copy constructed tree.

        Returns:
        -----------
            Self: A copy of this node with recursively copied children if `children` is not set with an appropriate mapping.
        """
        assert data is None, "No data must be set on a root node"
        if 'name' not in kwargs:
            kwargs['name'] = self._name
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
                "Children provided to `construct_copy` for datagroup are not of type `DataGroup` or `DataLeaf"
            )
            # We have new children and can extract the ResType from them
            new_dtype: type[ResType] | TypeForm[ResType] | None = dtype

            return DataGroup(
                children=children,
                dtype=new_dtype,
                **kwargs,
            )

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
    def is_flat_group(self) -> bool:
        """Boolean flag that is true if there are no more sub-groups beneath this group, thus making the children of this group exclusively data-nodes."""
        return len(self.subgroups) == 0

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

    def group_children_by(
        self,
        key_func: Callable[["TreeNode"], KeyType | None],
        group_leaves_only: bool = False,
    ) -> Self | None:
        """Specialization of the `group_children_by` function for group nodes, where grouping may need to be
        performed on subsets of their children.

        Returns
        -------
        Self
            Generally returns the same node type, potentially with updated children and an additional layer of `DataGroup` nodes underneath
        """
        # At the end of this, we should have either only sub-groups or only sub-leaves
        num_categories = 0
        key_set: set[KeyType | str] = set()
        member_children: Mapping[
            KeyType | str,
            list[tuple[Hashable, DataGroup[DataType] | DataLeaf[DataType]]],
        ] = {}

        res_children: Mapping[Hashable, DataGroup[DataType] | DataLeaf[DataType]] = {}

        for k, child in self.children.items():
            # If we recurse, group the child first.
            child = child.group_children_by(
                key_func=key_func, group_leaves_only=group_leaves_only
            )

            if isinstance(child, DataGroup):
                if group_leaves_only:
                    res_children[k] = child
                    num_categories += 1
                else:
                    key = key_func(child)
                    if key is None:
                        continue

                    if key not in key_set:
                        key_set.add(key)
                        member_children[key] = []
                        num_categories += 1
                    member_children[key].append((k, child))
            elif isinstance(child, DataLeaf):
                key = key_func(child)
                if key is None:
                    continue
                if key not in key_set:
                    key_set.add(key)
                    member_children[key] = []
                    num_categories += 1
                member_children[key].append((k, child))

        new_children = res_children
        base_group_info = (
            self._group_info
            if self._group_info is not None
            else GroupInfo(group_name=self._name or "group")
        )

        # TODO: FIXME: Make key to group info more straightforward

        for key, group in member_children.items():
            try:
                # First try to treat it as a dataclass
                key_dict = asdict(key)  # type: ignore # We know that the type does not need to fit, but it will raise a TypeError that we handle right afterwards if we are wrong
            except TypeError:
                # It wasn't a dataclass apparently
                key_dict = {'key': key}

            group_child_dict: dict[
                Hashable, DataGroup[DataType] | DataLeaf[DataType]
            ] = {e[0]: e[1] for e in group}

            if num_categories == 1:
                # Only one category, update the group info and return full node
                base_group_info.group_name = str(key)
                base_group_info.group_attributes = key_dict
                new_children.update(group_child_dict)
            else:
                # Generate new group for this category
                new_group_info = GroupInfo(str(key), group_attributes=key_dict)
                new_group = DataGroup(
                    group_info=new_group_info,
                    children=group_child_dict,
                    dtype=self._dtype,
                )
                for i in range(10000):
                    group_name_try = f"group_{i}"
                    if group_name_try not in new_children:
                        new_children[group_name_try] = new_group

        return self.construct_copy(children=new_children, group_info=base_group_info)
