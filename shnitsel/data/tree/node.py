import abc
from dataclasses import dataclass
from typing import Any, Callable, Hashable, Mapping, Self, TypeVar, Generic, get_args


ChildType = TypeVar("ChildType", bound="TreeNode|None")
DataType = TypeVar("DataType")
NewDataType = TypeVar("NewDataType")
NewChildType = TypeVar("NewChildType", bound="TreeNode|None")
ResType = TypeVar("ResType")


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

    def construct_copy(self, **kwargs) -> Self:
        return type(self)(
            name=self._name,
            data=self.data,
            children=self.children,
            attrs=dict(self._attrs),
            level_name=self._level_name,
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

    # def map_node(self, func: Callable[[Self], ResType]) -> ResType:
    #     return func(self)

    @abc.abstractmethod
    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "TreeNode|None": ...

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
            tmp_res = self.construct_copy()
            return tmp_res

    def is_level(self, target_level: str) -> bool:
        """Check whether we are at a certain level

        Args:
            target_level (str): Desired level to check for

        Returns:
            bool: True if this level satisfies the requirements
        """
        return self._level_name == target_level
