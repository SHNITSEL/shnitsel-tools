from dataclasses import dataclass
from typing import Any, Callable, Hashable, Mapping, Self, TypeVar, Generic


ChildType = TypeVar("ChildType", bound="TreeNode|None")
DataType = TypeVar("DataType")
ResType = TypeVar("ResType")


@dataclass
class TreeNode(Generic[ChildType, DataType]):
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
            level_name if level_name is not None else self.__class__.__name__
        )

    def copy(self, copy_children: bool = False, new_type=None) -> Self:
        new_children = None
        if copy_children:
            new_children = {}
            for key, child in self._children.items():
                assert isinstance(child, TreeNode), (
                    "Can only copy children of TreeNode type."
                )
                new_children[key] = child.copy(copy_children=True)

        new_type = type(self) if new_type is None else new_type
        return new_type(
            name=self._name,
            data=self._data,
            children=new_children,
            attrs=dict(self._attrs),
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

    def map_over_child_nodes(
        self, func: Callable[[Self], ResType | None]
    ) -> Mapping[Hashable, ResType]:
        new_children = {
            k: res
            for k, v in self._children.items()
            if v is not None and (res := v.map_node(func)) is not None
        }
        return new_children

    def map_node(self, func: Callable[[Self], ResType]) -> ResType:
        return func(self)

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "TreeNode[ChildType, ResType] | None":
        new_data = self._data
        if self.has_data:
            new_data = func(self.data)

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

        if not keep_empty_branches and new_children is None and new_data is None:
            return None
        else:
            # This yields a different kind of tree.
            tmp_res = self.copy()
            tmp_res._data = new_data
            if recurse:
                tmp_res._children = new_children if new_children is not None else {}
            return tmp_res

    def filter_nodes(
        self,
        filter_func: Callable[..., bool],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> Self | None:
        keep_self = filter_func(self)

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
            tmp_res = self.copy()
