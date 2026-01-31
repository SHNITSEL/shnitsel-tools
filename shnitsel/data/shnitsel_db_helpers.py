from types import UnionType
from typing import Any, Callable, List, Literal, TypeVar
import xarray as xr

from shnitsel.data.dataset_containers import Frames, Trajectory
from shnitsel.data.traj_combiner_methods import concat_trajs, layer_trajs
from shnitsel.data.tree import DataLeaf, TreeNode
from .tree.datatree_level import (
    _datatree_level_attribute_key,
    DataTreeLevelMap,
)

DataType = TypeVar("DataType")
AggType = TypeVar("AggType")
R = TypeVar("R", bound=xr.Dataset)
TrajType = TypeVar("TrajType", bound=Trajectory | Frames | xr.Dataset)


def concat_subtree(
    subtree: TreeNode[Any, TrajType], only_direct_children: bool = False
) -> xr.Dataset:
    """Helper function to concatenate the trajectories in a subtree into a multi-trajetctory dataset.

    The resulting dataset has a new `frame` dimension along which we can iterate through all individual frames of all trajectories.

    Parameters
    ----------
    subtree : TreeNode[Any, TrajType]
        The subtree of the ShnitselDB datastructure
    only_direct_children : bool, optional
        Whether to only gather trajectories from direct children of this subtree.

    Returns
    -------
    xr.Dataset
        The resulting multi-trajectory dataset
    """
    if not only_direct_children:
        trajectories = list(subtree.collect_data())
    else:
        trajectories = [
            x.data
            for x in subtree.children.values()
            if x.is_leaf and x.has_data and x.data is not None
        ]

    return concat_trajs(trajectories)


def layer_subtree(
    subtree: TreeNode[Any, Trajectory | Frames | xr.Dataset],
    only_direct_children: bool = False,
) -> xr.Dataset:
    """Helper function to layer the trajectories in a subtree into a multi-trajetctory dataset.

    The resulting dataset has a new `trajectory` dimension along which we can iterate through all individual frames of all trajectories.

    Parameters
    ----------
    subtree : TreeNode[Any, TrajTrajectoryType]
        The subtree of the ShnitselDB datastructure
    only_direct_children : bool, optional
        Whether to only gather trajectories from direct children of this subtree, by default false.

    Returns
    -------
    xr.Dataset
        The resulting multi-trajectory dataset
    """
    if not only_direct_children:
        trajectories = list(subtree.collect_data())
    else:
        trajectories = [
            x.data for x in subtree.children.values() if x.is_leaf and x.has_data
        ]
    return layer_trajs(trajectories)


def list_subtree(
    subtree: TreeNode[Any, DataType], only_direct_children: bool = False
) -> list[DataType]:
    """Helper function to collect the data in a subtree into a list.

    Parameters
    ----------
    subtree : TreeNode[Any, DataType]
        The subtree of the ShnitselDB datastructure holding the DataType entries to collect
    only_direct_children : bool, optional
        Whether to only gather trajectories from direct children of this subtree, by default false

    Returns
    -------
    list[DataType]
        The resulting data/trajectory list
    """
    if not only_direct_children:
        data = list(subtree.collect_data())
    else:
        data = [x.data for x in subtree.children.values() if x.is_leaf and x.has_data]
    return data


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
    dtype: type[R] | UnionType | None = None,
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


def get_data_with_path(subtree: TreeNode[Any, DataType]) -> List[tuple[str, DataType]]:
    """Function to get a list of all data in a tree with their respective path

    Parameters
    ----------
    subtree : TreeNode[Any, DataType]
        The subtree to generate the collection for.

    Returns
    -------
    List[tuple[str, DataType]]
        A list of tuples (path, data at that path) for all data in the subtree.
    """

    res = []
    if subtree.has_data:
        # the tree will give us empty datasets instead of none if an attribute on the node has been set.
        res.append((subtree.path, subtree.data))

    for key, child in subtree.children.items():
        child_res = get_data_with_path(child)
        if child_res is not None and len(child_res) > 0:
            res = res + child_res

    return res
