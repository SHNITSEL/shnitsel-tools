from types import UnionType
from typing import Any, Callable, List, Literal, TypeVar
import xarray as xr

from shnitsel.data.tree import DataLeaf, TreeNode
from shnitsel.data.trajectory_format import Trajectory
from .tree.datatree_level import (
    _datatree_level_attribute_key,
    DataTreeLevelMap,
)

DataType = TypeVar("DataType")
AggType = TypeVar("AggType")
R = TypeVar("R", bound=xr.Dataset)


def unwrap_single_entry_in_tree(
    tree: TreeNode[Any, DataType],
) -> DataType | TreeNode[Any, DataType]:
    """Attempts to unwrap a single dataset from a tree.

    If multiple or none are found, it will return the original tree
    If a single entry was found, will return the dataset

    Parameters:
    root : TreeNode[Any, DataType]
        Root of the subtree to parse

    Returns
    -------
    DataType|TreeNode[Any, DataType]
        Returns either the single point of data in the subtree or the full tree, if unwrapping would be unfeasible.
    """

    data = list(tree.collect_data())
    if len(data) == 1:
        return data[0]
    else:
        return tree


def aggregate_xr_over_levels(
    tree: TreeNode[Any, DataType],
    func: Callable[[TreeNode[Any, DataType]], R],
    level_name: Literal['root', 'compound', 'group', 'data'],
    dtype: type[R] | UnionType | None = None
) -> TreeNode[Any, R] | None:
    """Apply an aggregation function to every node at a level of a db structure

    Parameters
    ----------
    tree : TreeNode[Any, DataType]
        The tree to aggregate at the specific level
    func : Callable[[TreeNode[Any, DataType]], R]
        The function to apply to the subtree at the specified level
    level_name : Literal['root', 'compound', 'group', 'data']
        The target level to apply the function `func` to.
        See `tree.datatree_level` for values.
    dtype : type | UnionType, optional
        The dtype of the resulting tree after aggregation.

    Returns
    -------
    TreeNode[Any, R]
        The resulting tree after applying the transform `func` to the subtrees.
    """
    new_children: dict[str, TreeNode[Any, R]] = {}
    drop_keys = []

    if tree.is_level(DataTreeLevelMap[level_name]):
        tmp_aggr: R = func(tree)
        tmp_label = f"aggregate of subtree({tree.name})"
        new_node = DataLeaf(name=tmp_label, data=tmp_aggr)
        new_children[tmp_label] = new_node

    for k, child in tree.children.items():
        child_res = aggregate_xr_over_levels(child, func, level_name)
        if child_res is not None:
            new_children[k] = child_res
        else:
            drop_keys.append(k)

    if len(new_children) > 0:
        return tree.construct_copy(children=new_children, dtype=dtype)
    else:
        return None


# TODO: FIXME: currently no path logic in tree.
# def get_trajectories_with_path(subtree: T) -> List[tuple[str, Trajectory]]:
#     """Function to get a list of all datasets in the tree with their respective path

#     Args:
#         subtree (xr.DataTree): The subtree to generate the collection for.

#     Returns:
#         List[tuple[str, Trajectory]]: A list of tuples (path, dataset at that path) for all datasets in the respective subtree.
#     """
#     # TODO: FIXME: This needs to be a bit more generalized for trees with arbitrary data

#     res = []
#     if subtree.has_data:
#         # the tree will give us empty datasets instead of none if an attribute on the node has been set.
#         res.append((subtree.path, subtree.dataset))

#     for key, child in subtree.children.items():
#         child_res = get_trajectories_with_path(child)
#         if child_res is not None and len(child_res) > 0:
#             res = res + child_res

#     return res
