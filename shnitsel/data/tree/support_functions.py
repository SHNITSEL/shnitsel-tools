from types import UnionType
from typing import Any, TypeVar, Union, overload
from typing_extensions import TypeForm

from .data_leaf import DataLeaf
from .node import TreeNode
from .tree import ShnitselDBRoot
from .data_group import DataGroup
from .compound import CompoundGroup

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

            if res_data_type is None:
                new_data_type = tuple(data_types)
            else:
                new_data_type = res_data_type

            new_data = tuple(res_data_entries)
            return DataLeaf(
                name=tree_list[0].name,
                data=new_data,
                dtype=new_data_type,
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
