import abc
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Callable, Hashable, Mapping, Self, TypeVar, Generic, overload


ChildType = TypeVar("ChildType", bound="TreeNode|None", covariant=True)
DataType = TypeVar("DataType", covariant=True)
NewDataType = TypeVar("NewDataType")
NewChildType = TypeVar("NewChildType", bound="TreeNode|None")
ResType = TypeVar("ResType")
KeyType = TypeVar("KeyType")


@dataclass
class TreeNode(Generic[ChildType, DataType], abc.ABC):
    _name: str | None

    _data: DataType | None
    _children: Mapping[Hashable, ChildType]
    _attrs: Mapping[str, Any]

    _parent: Self | None
    _level_name: str | None

    def __init__(
        self,
        name: str | None,
        data: DataType | None = None,
        children: Mapping[Hashable, ChildType] | None = None,
        attrs: Mapping[str, Any] | None = None,
        level_name: str | None = None,
    ):
        self._name = name
        self._data = data
        self._children = children if children is not None else dict()
        self._attrs = attrs if attrs is not None else dict()
        self._parent = None
        self._level_name = (
            level_name if level_name is not None else self.__class__.__qualname__
        )

    def construct_copy(
        self,
        **kwargs,
    ) -> Self | "TreeNode[Any, ResType]":
        if 'name' not in kwargs:
            kwargs['name'] = self._name
        if 'data' not in kwargs:
            kwargs['data'] = self._data
        if 'children' not in kwargs:
            kwargs['children'] = {
                k: v.construct_copy()
                for k, v in self._children.items()
                if v is not None
            }
        if 'attrs' not in kwargs:
            kwargs['attrs'] = dict(self._attrs)
        if 'level_name' not in kwargs:
            kwargs['level_name'] = str(self._level_name)

        return type(self)(
            **kwargs,
        )

    @property
    def is_leaf(self) -> bool:
        return len(self._children) == 0

    @property
    def has_data(self) -> bool:
        return self._data is not None

    @property
    def data(self) -> DataType:
        if self._data is None:
            raise ValueError("Object has no data to be retrieved")
        else:
            return self._data

    @property
    def children(self) -> Mapping[Hashable, ChildType]:
        return self._children

    @property
    def attrs(self) -> Mapping[str, Any]:
        return self._attrs

    @property
    def name(self) -> str:
        if self._name is None:
            raise KeyError("Node has no `name` attribute set.")
        else:
            return self._name

    # def map_over_child_nodes(
    #     self, func: Callable[[Self], ResType | None]
    # ) -> Mapping[Hashable, ResType]:
    #     new_children = {
    #         k: res
    #         for k, v in self._children.items()
    #         if v is not None and (res := v.map_node(func)) is not None
    #     }
    #     return new_children

    def map_subtree(self, func: Callable[[Self], ResType]) -> ResType:
        return func(self)

    @abc.abstractmethod
    def group_children_by(
        self,
        key_func: Callable[["TreeNode"], KeyType],
        group_leaves_only: bool = False,
    ) -> Self | None: ...

    @abc.abstractmethod
    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "TreeNode|None": ...

    def map_filtered_nodes(
        self,
        filter_func: Callable[["TreeNode[Any, DataType]"], bool],
        map_func: Callable[["TreeNode[Any, DataType]"], "TreeNode[Any, ResType]"],
    ) -> "TreeNode[Any, ResType]":
        """Map nodes if the filter function picks them as relevant to this run.
        If the node is not picked by `filter_func` a copy will be created with its children being recursively mapped
        according to the same rule.
        If a node is mapped, the `map_func` must take care of potential mapping over children."""

        if filter_func(self):
            new_node = map_func(self)
        else:
            new_children = {
                k: res
                for k, v in self.children.items()
                if v is not None
                and (res := v.map_filtered_nodes(filter_func, map_func)) is not None
            }

            new_node: TreeNode[Any, ResType] = self.construct_copy(
                children=new_children
            )  # type: ignore # By mapping the children, we are sure that they now also hold the resulting datatype.

        return new_node

    def filter_nodes(
        self,
        filter_func: Callable[..., bool],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> Self | None:
        keep_self = filter_func(self)

        # Stop if the node is not kept.
        if not keep_self:
            return None

        new_children = None
        if recurse:
            new_children = {
                k: res
                for k, v in self._children.items()
                if v is not None and (res := v.filter_nodes(filter_func)) is not None
            }

            if len(new_children) == 0:
                new_children = None

        if not keep_empty_branches and not keep_self and new_children is None:
            return None
        else:
            tmp_res = self.construct_copy(children=new_children)
            return tmp_res

    def add_child(self, child_name: str | None, child: ChildType) -> Self:
        """Add a new child node with a preferred name in the mapping of children.
        If the child name is already in use, will attempt to find a collision-free alternative name.

        Parameters
        -----------
            child_name (str | None): The preferred name under which the child should be registered.
                    To avoid overriding, a different name will be chosen if the name is in use.
            child (ChildType): Object to register as the child-subtree

        Raises
        -----------
            OverflowError: If the attempts to find a new collision-free name have exceeded 1000.

        Returns
        -----------
            Self: The new instance of a subtree
        """
        new_children = dict(self._children)
        if child_name is not None and child_name not in new_children:
            new_children[child_name] = child
        else:
            if child_name is None:
                child_name = type(child).__name__

            found = False
            for i in range(1000):
                tmp_name = child_name + "_" + str(i)
                if tmp_name not in new_children:
                    found = True
                    new_children[tmp_name] = child
                    break

            if not found:
                raise OverflowError(
                    "Could not patch child name without name collision after 1000 modifications"
                )
        return self.construct_copy(children=new_children)

    def assign_children(self, new_children: Mapping[Hashable, ChildType]) -> Self:
        """Construct a new instance of a subtree with new children assigned to this tree"""
        # TODO: FIXME: Implement
        all_children = dict(self._children)
        all_children.update(new_children)
        return self.construct_copy(children=all_children)

    def is_level(self, target_level: str) -> bool:
        """Check whether we are at a certain level

        Args:
            target_level (str): Desired level to check for

        Returns:
            bool: True if this level satisfies the requirements
        """
        return self._level_name == target_level

    def collect_data(self) -> Iterable[DataType]:
        """Get an iterator over all data in this subtree. Does not preserve the tree structure but can be helpful for aggregation functions."""
        if self.has_data:
            yield self.data

        for child in self.children.values():
            if child is not None:
                yield from child.collect_data()
