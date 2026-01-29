from types import UnionType
from typing import Any, Hashable, Iterable, Sequence, TypeVar, Union, overload
from typing_extensions import TypeForm

from .data_leaf import DataLeaf
from .node import TreeNode
from .tree import ShnitselDBRoot
from .data_group import DataGroup
from .compound import CompoundGroup

DataType = TypeVar("DataType")
DataType1 = TypeVar("DataType1")
DataType2 = TypeVar("DataType2")
DataType3 = TypeVar("DataType3")
DataType4 = TypeVar("DataType4")
ResDataType = TypeVar("ResDataType")

NodeType = TypeVar("NodeType", bound=TreeNode)
ChildType = TypeVar("ChildType", bound=TreeNode)


@overload
def tree_zip(
    *trees: DataLeaf, res_data_type: type[ResDataType] | TypeForm[ResDataType]
) -> DataLeaf[ResDataType] | None: ...


@overload
def tree_zip(*trees: DataLeaf, res_data_type: None = None) -> DataLeaf | None: ...


@overload
def tree_zip(
    *trees: CompoundGroup, res_data_type: type[ResDataType] | TypeForm[ResDataType]
) -> CompoundGroup[ResDataType] | None: ...


@overload
def tree_zip(
    *trees: CompoundGroup, res_data_type: None = None
) -> CompoundGroup | None: ...


@overload
def tree_zip(
    *trees: DataGroup, res_data_type: type[ResDataType] | TypeForm[ResDataType]
) -> DataGroup[ResDataType] | None: ...


@overload
def tree_zip(*trees: DataGroup, res_data_type: None = None) -> DataGroup | None: ...


@overload
def tree_zip(
    *trees: ShnitselDBRoot, res_data_type: type[ResDataType] | TypeForm[ResDataType]
) -> ShnitselDBRoot[ResDataType] | None: ...


@overload
def tree_zip(
    *trees: ShnitselDBRoot, res_data_type: None = None
) -> ShnitselDBRoot | None: ...


@overload
def tree_zip(
    *trees: TreeNode,
    res_data_type: type[ResDataType] | TypeForm[ResDataType] | None = None,
) -> TreeNode | TreeNode[Any, ResDataType] | None: ...


def tree_zip(
    *trees: TreeNode,
    res_data_type: type[ResDataType] | TypeForm[ResDataType] | None = None,
) -> TreeNode | TreeNode[Any, ResDataType] | None:
    """Helper function to allow zipping of multiple trees into a single tree with tuples of data for
    its data.

    The zipping is only performed on the data, metadata will be taken from the tree provided first.
    If provided with a `res_data_type`, the data type for the resulting tree will be set accordingly

    The resulting data tuples will hold data from the various trees in order.

    Parameters
    ----------
    *trees: TreeNode
        An arbitrary positional list of trees to use for the zipping.
    res_data_type : type[ResDataType] | TypeForm[ResDataType] | None, optional
        Optional datatype for the resulting tree, by default None, which means, it will be inferred.

    Returns
    -------
    TreeNode | TreeNode[Any, ResDataType] | None
        The tree node of the same type as the root in the first provided tree but with an updated
        DataType.
        If no zipping was possible, because no trees were provided, None is returned.

    Raises
    ------
    ValueError
        If trees with inconsistent structure were provided
    """
    tree_list: list[TreeNode] = list(trees)

    if len(tree_list) == 0:
        return None

    if not has_same_structure(*tree_list):
        raise ValueError(
            "Trees provided to `zip` were not of same structure. Zipping impossible."
        )

    res_data_entries = []

    # TODO: Build tuple of child types and explictly set dtype of new tree?

    child_keys: set[str] | None = None
    has_data: bool | None = None

    for tree in tree_list:
        child_keys = set(str(k) for k, v in tree.children.items() if v is not None)
        has_data = tree.has_data
        break

    if isinstance(tree_list[0], DataLeaf):
        if has_data:
            # We are a data leaf:
            data_types = []
            for tree in tree_list:
                tree_data = tree.data
                res_data_entries.append(tree_data)
                data_types.append(type(tree_data))

            # TODO: FIXME: Figure out how to do correct tuple typing.
            # if res_data_type is None:
            #     new_data_type = tuple[*data_types]
            # else:
            #     new_data_type = res_data_type

            new_data = tuple(res_data_entries)
            return DataLeaf(
                name=tree_list[0].name,
                data=new_data,
                # dtype=new_data_type,
                attrs=tree_list[0]._attrs,
            )
        else:
            return tree_list[0].construct_copy()

    assert (
        isinstance(tree_list[0], ShnitselDBRoot)
        or isinstance(tree_list[0], CompoundGroup)
        or isinstance(tree_list[0], DataGroup)
    ), "Unsupported node type provided to `tree_zip() : %s" % type(tree_list[0])

    new_children = {}
    if child_keys:
        for key in child_keys:
            if res_data_type is not None:
                new_children[key] = tree_zip(
                    *[tree.children[key] for tree in tree_list],
                    res_data_type=res_data_type,
                )
            else:
                # TODO: Figure out, why this gives us a type check warning.
                new_children[key] = tree_zip(
                    *[tree.children[key] for tree in tree_list]
                )
    if res_data_type is not None:
        return tree_list[0].construct_copy(children=new_children, dtype=res_data_type)
    else:
        return tree_list[0].construct_copy(children=new_children)


def has_same_structure(*trees: TreeNode) -> bool:
    """Function to check whether a set of trees has the same overall structure

    This means, they must have same keys to not-None children at every level and data in nodes along the same path.

    Returns
    -------
    bool
        True if all tree structures match, False otherwise.
    """
    child_keys: set[str] | None = None
    has_data: bool | None = None

    tree_list = list(trees)

    for tree in tree_list:
        tree_keys = set(str(k) for k, v in tree.children.items() if v is not None)
        if child_keys is None:
            child_keys = tree_keys
        else:
            if child_keys.symmetric_difference(tree_keys):
                return False
        if has_data is None:
            has_data = tree.has_data
        else:
            if has_data != tree.has_data:
                return False

    if child_keys is None:
        return True

    for child_key in child_keys:
        child_nodes = [tree.children[child_key] for tree in tree_list]
        if not has_same_structure(*child_nodes):
            return False
    return True


@overload
def tree_merge(
    *trees: ShnitselDBRoot[DataType],
    res_data_type: type[DataType] | UnionType | None = None,
) -> ShnitselDBRoot[DataType] | None: ...


@overload
def tree_merge(
    *trees: CompoundGroup[DataType],
    res_data_type: type[DataType] | UnionType | None = None,
) -> CompoundGroup[DataType] | None: ...


@overload
def tree_merge(
    *trees: DataGroup[DataType],
    res_data_type: type[DataType] | UnionType | None = None,
) -> DataGroup[DataType] | None: ...


@overload
def tree_merge(
    *trees: TreeNode[Any, DataType],
    res_data_type: type[DataType] | UnionType | None = None,
) -> TreeNode[Any, DataType] | None: ...


def tree_merge(
    *trees: ShnitselDBRoot[DataType]
    | CompoundGroup[DataType]
    | DataGroup[DataType]
    | TreeNode[Any, DataType],
    res_data_type: type[DataType] | UnionType | None = None,
) -> (
    ShnitselDBRoot[DataType]
    | CompoundGroup[DataType]
    | DataGroup[DataType]
    | TreeNode[Any, DataType]
    | None
):
    """Helper function to merge two trees at the same level.
    Data leaves on the same level will all be retained.
    Data Group children of the roots will be merged recursively.


    Parameters
    ----------
    *trees: ShnitselDBRoot[DataType] | CompoundGroup[DataType] | DataGroup[DataType] | TreeNode[Any, DataType]
        Compatible roots at the same level that represent a group of children.
        If inconsistent types are provided, the merge may fail.
    res_data_type : type[DataType] | TypeForm[DataType] | None, optional
        An explicit indicator of which type we expect the merged tree to have, by default None

    Returns
    -------
    ShnitselDBRoot[DataType] | CompoundGroup[DataType] | DataGroup[DataType] | TreeNode[Any, DataType] | None
        The merged tree of the same level as the input tree roots.
        Specifically, the same level as `trees[0]`.
        If there are no `trees`, then `None` is returned.
        If a single `trees` parameter is provided, then a copy of that tree is returned.

    Raises
    ------
    ValueError
        _description_
    """
    from .child_support_functions import find_child_key

    if len(trees) == 0:
        return None

    if len(trees) == 1:
        return trees[0].construct_copy(dtype=res_data_type)

    data_children: Sequence[DataLeaf[DataType]] = []

    children_to_merge: dict[Hashable, list[DataGroup[DataType]]] = dict()

    for tree in trees:
        data_children = data_children + [
            child.construct_copy(dtype=res_data_type)
            for child in tree.children.values()
            if isinstance(child, DataLeaf)
        ]
        for key, child in tree.children.items():
            if isinstance(child, DataGroup):
                if key not in children_to_merge:
                    children_to_merge[key] = [child]
                else:
                    children_to_merge[key].append(child)

    children_res: dict[Hashable, DataGroup[DataType] | DataLeaf[DataType]] = {
        key: res
        for key, candidates in children_to_merge.items()
        if (res := tree_merge(*candidates, res_data_type=res_data_type)) is not None
        and isinstance(res, DataGroup)
    }

    for data_child in data_children:
        child_key = find_child_key(children_res.keys(), data_child)
        children_res[child_key] = data_child

    if isinstance(trees[0], ShnitselDBRoot):
        return ShnitselDBRoot[DataType](compounds=children_res, dtype=res_data_type)  # type: ignore # If the prior Compound groups were right, this is too
    elif isinstance(trees[0], CompoundGroup):
        return trees[0].construct_copy(children=children_res, dtype=res_data_type)
    elif isinstance(trees[0], DataGroup):
        return trees[0].construct_copy(children=children_res, dtype=res_data_type)
    else:
        raise ValueError("Cannot merge type of provided root: %s" % type(trees[0]))
