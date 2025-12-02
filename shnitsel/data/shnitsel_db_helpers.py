from typing import Any, Callable, List, TypeVar
import xarray as xr

from shnitsel.data.shnitsel_db.db_trajectory_data import TrajectoryData
from shnitsel.data.trajectory_format import Trajectory
from .shnitsel_db.datatree_level import (
    _datatree_level_attribute_key,
    DataTreeLevel_keys,
)

T = TypeVar("T", bound=xr.DataTree)
R = TypeVar("R", bound=xr.Dataset)


def unwrap_single_entry_in_tree(tree: T) -> Trajectory | T:
    """Attempts to unwrap a single dataset from a tree.

    If multiple or none are found, it will return the original tree
    If a single entry was found, will return the dataset

    Args:
        root (xr.DataTree): Root of the subtree to parse

    Returns:
        xr.Dataset|List[Any]|None: Returns None if no entry was found, a list instance if multiple entries were found or a single dataset if a single entry was found in the subtree
    """

    def collect_single_data(root: xr.DataTree) -> Trajectory | List[Any] | None:
        """Returns the single dataset in the subtree

        Args:
            root (xr.DataTree): Root of the subtree to parse

        Returns:
            xr.Dataset|List[Any]|None: Returns None if no entry was found, a list instance if multiple entries were found or a single dataset if a single entry was found in the subtree
        """
        res = None

        if root.dataset is not None:
            res = root.dataset

        for c_k, child in root.children.items():
            child_res = collect_single_data(child)
            if child_res is None:
                continue
            elif isinstance(child_res, list):
                return child_res
            elif isinstance(child_res, Trajectory):
                if res is None:
                    res = child_res
                else:
                    return []

        return res

    res_unwrap = collect_single_data(tree)
    if isinstance(res_unwrap, Trajectory):
        return res_unwrap
    else:
        return tree


def aggregate_xr_over_levels(tree: T, func: Callable[[T], R], level: str) -> T | None:
    """Apply an aggregation function to every node at a level of a db structure

    Args:
        tree (T): The tree to aggregate at the specific level
        func (callable): The function to apply to that subtree
        level (str): The target level to apply the function `func` to. See `shnitsel_db.datatree_level.py` for values.

    Returns:
        T: The resulting tree after applying the transform `func` to the subtrees.
    """
    new_children = {}
    drop_keys = []

    if tree.attrs[_datatree_level_attribute_key] == DataTreeLevel_keys[level]:
        tmp_aggr: R = func(tree)
        tmp_label = f"aggregate of subtree({tree.name})"
        new_node = TrajectoryData(tmp_aggr, tmp_label)
        new_children[tmp_label] = new_node

    for k, child in tree.children.items():
        child_res = aggregate_xr_over_levels(child, func, level)
        if child_res is not None:
            new_children[k] = child_res
        else:
            drop_keys.append(k)

    if len(new_children) > 0:
        return tree.copy().drop_nodes(drop_keys).assign(new_children)
    else:
        return None


def get_trajectories_with_path(subtree: xr.DataTree) -> List[tuple[str, Trajectory]]:
    """Function to get a list of all datasets in the tree with their respective path

    Args:
        subtree (xr.DataTree): The subtree to generate the collection for.

    Returns:
        List[tuple[str, Trajectory]]: A list of tuples (path, dataset at that path) for all datasets in the respective subtree.
    """
    # TODO: FIXME: This needs to be a bit more generalized for trees with arbitrary data
    
    res = []
    if subtree.dataset is not None:
        res.append((subtree.path, subtree.dataset))

    for key, child in subtree.children.items():
        child_res = get_trajectories_with_path(child)
        if child_res is not None and len(child_res) > 0:
            res = res + child_res

    return res
