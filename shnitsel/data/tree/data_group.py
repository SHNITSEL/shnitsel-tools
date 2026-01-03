from dataclasses import dataclass, asdict
from functools import cached_property
from typing import Any, Callable, Generic, Hashable, Mapping, Self, TypeVar
from .node import TreeNode
from .data_leaf import DataLeaf

DataType = TypeVar("DataType", covariant=True)
ResType = TypeVar("ResType")
KeyType = TypeVar("KeyType")


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
        name: str | None = None,
        group_info: GroupInfo | None = None,
        children: Mapping[
            Hashable,
            "DataGroup[DataType]|DataLeaf[DataType]",
        ]
        | None = None,
        attrs: Mapping[str, Any] | None = None,
        level_name: str | None = None,
    ):
        if name is None and group_info is not None:
            name = group_info.group_name

        super().__init__(
            name=name,
            data=None,
            children=children,
            attrs=attrs,
            level_name=level_name,
        )
        self._group_info = group_info

    def construct_copy(self, **kwargs) -> Self:
        if 'group_info' not in kwargs:
            kwargs['group_info'] = self._group_info
        return super().construct_copy(**kwargs)

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

    @cached_property
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
                key_dict = asdict(key)
            except:
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
                new_group = DataGroup[DataType](
                    group_info=new_group_info, children=group_child_dict
                )
                for i in range(10000):
                    group_name_try = f"group_{i}"
                    if group_name_try not in new_children:
                        new_children[group_name_try] = new_group

        return self.construct_copy(children=new_children, group_info=base_group_info)

    def map_data(
        self,
        func: Callable[[DataType], ResType | None],
        recurse: bool = True,
        keep_empty_branches: bool = False,
    ) -> "DataGroup[ResType]|None":
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
            return DataGroup[ResType](
                name=self._name,
                group_info=self._group_info,
                children=new_children,
                level_name=self._level_name,
                attrs=dict(self.attrs),
            )
