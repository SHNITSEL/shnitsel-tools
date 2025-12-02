from dataclasses import dataclass
import logging
from typing import Callable, List, Mapping, Self, TypeVar

from .helpers import traj_list_to_child_mapping
from shnitsel.data.trajectory_format import Trajectory

from .db_trajectory_group import GroupInfo, TrajectoryGroup

from .db_trajectory_data import TrajectoryData
import xarray as xr

from .datatree_level import _datatree_level_attribute_key

T = TypeVar("T")


@dataclass
class CompoundInfo:
    """Class to hold identifying and auxiliary info of a compound type in ShnitselDB"""

    compound_name: str = "unknown"


class CompoundGroup(xr.DataTree):
    """DataTree node to keep track of all data associated with a common compound within the datatree"""

    def __init__(
        self,
        compound_info: CompoundInfo | None = None,
        children: Mapping[
            str,
            TrajectoryGroup | TrajectoryData,
        ]
        | None = None,
    ):
        super().__init__(
            None,
            children,
            compound_info.compound_name if compound_info is not None else None,
        )
        self.attrs[_datatree_level_attribute_key] = "CompoundGroup"
        if compound_info is not None:
            self.attrs["compound_info"] = compound_info.__dict__

    def is_level(self, target_level: str) -> bool:
        """Check whether we are at a certain level

        Args:
            target_level (str): Desired level to check for

        Returns:
            bool: True if this level satisfies the requirements
        """
        return target_level == "TrajectoryGroup" or target_level == "CompoundGroup"

    def get_compound_info(self) -> CompoundInfo:
        """Get the store compound info of this Compound group.

        Returns:
            CompoundInfo: _description_
        """
        if "compound_info" in self.attrs:
            compound_data = self.attrs["compound_info"]
            return CompoundInfo(**compound_data)
        else:
            return CompoundInfo()

    def collect_trajectories(self) -> List[TrajectoryData]:
        """Function to retrieve all trajectories in this subtree

        Returns:
            List[TrajectoryData]: List of all nodes with TrajectoryData type
        """
        res = []

        for x in self.children.values():
            res += x.collect_trajectories()

        return res

    def merge_with(self, other: Self) -> Self:
        """Function to merge two compound groups into one.

        Called when merging two database states. Will fail if compound_info differs between compounds to avoid loss of information.

        Args:
            other (CompoundGroup): The other CompoundGroup to be merged

        Raises:
            ValueError: Raised if the compound_info differs.

        Returns:
            CompoundGroup: A CompoundGroup object holding the entire merged subtree
        """
        own_info = self.get_compound_info()
        other_info = other.get_compound_info()
        if other_info != own_info:
            message = f"Cannot merge compounds with conflicting compound information: {other_info} vs. {own_info}"
            logging.error(message)
            raise ValueError(message)

        res_children = []

        for k, v in self.children.items():
            key_str = str(k)
            if key_str in other.children:
                other_v = other.children[key_str]
                if isinstance(v, TrajectoryGroup) and isinstance(
                    other_v, TrajectoryGroup
                ):
                    res_children.append(v.merge_with(other_v))
                else:
                    res_children.append(other_v.copy())
            else:
                res_children.append(v.copy())

        for k, v in other.children.items():
            key_str = str(k)
            if key_str in self.children:
                continue
            else:
                res_children.append(v.copy())

        res_children = traj_list_to_child_mapping(res_children)

        return type(self)(own_info, res_children)  # type: ignore

    def filter_trajectories(
        self,
        filter_func: Callable[[xr.Dataset], bool] | None = None,
        est_level: str | List[str] | None = None,
        basis_set: str | List[str] | None = None,
        **kwargs,
    ) -> Self | None:
        """Function to filter trajectories based on their attributes.

        Args:
            filter_func (Callable[[xr.Dataset], bool] | None, optional): A function to evaluate whether a trajectory should be retained. Should return True if the trajectory should stay in the filtered set. Defaults to None.
            est_level (str | List[str] | None, optional): Option to filter for a certain level of electronic structure theory/calculation method. Can be a single key value or a set of values to retain. Defaults to None.
            basis_set (str | List[str] | None, optional): Option to filter for a certain basis set. Can be a single key value or a set of values to retain. Defaults to None.
            **kwargs: Key-value pairs, where the key denotes an attribute

        Returns:
            Compoundgroup|None: Either returns the CompoundGroup with the remaining set of trajectories or None if the group would be empty.
        """

        if filter_func is None:
            filter_func = lambda x: True

        if isinstance(est_level, str):
            est_level = list(est_level)

        if isinstance(basis_set, str):
            basis_set = list(basis_set)

        filter_vals = {
            k: list(v) if isinstance(v, str) else v for k, v in kwargs.items()
        }
        filter_vals["est_level"] = est_level
        filter_vals["basis_set"] = basis_set

        def composed_filter(data: xr.Dataset) -> bool:
            if not filter_func(data):
                return False

            for k, v in filter_vals.items():
                if k not in data.attrs or data.attrs[k] not in v:
                    return False

            return True

        filtered_traj = [
            t
            for t in self.collect_trajectories()
            if t.has_data and composed_filter(t.dataset)
        ]

        if len(filtered_traj) > 0:
            return type(self)(
                self.attrs["compound_info"],
                {f"{i}": v for i, v in enumerate(filtered_traj)},
            )
        else:
            return None

    def add_trajectory_group(
        self,
        group_name: str,
        filter_func_trajectories: Callable[[Trajectory | GroupInfo], bool]
        | None = None,
        flatten_trajectories=False,
        **kwargs,
    ) -> Self:
        """Function to add trajectories within this compound subtree to a `TrajectoryGroup` of trajectories.

        The `group_name` will be set as the name of the group in the tree.
        If `flatten_trajectories=True` all existing groups will be dissolved before filtering and the children will be turned into an ungrouped list of trajectories.
        The `filter_func_trajectories` will either be applied to only the current groups and trajectories immediately beneath this compound or to the flattened list of all child directories.

        Args:
            group_name (str): The name to be set for the TrajectoryGroup object
            filter_func_Trajectories (Callable[[Trajectory|GroupInfo], bool] | None, optional): A function to return true for Groups and individual trajectories that should be added to the new group. Defaults to None.
            flatten_trajectories (bool, optional): A flag whether all descendant groups should be dissolved and flattened into a list of trajectories first before applying a group. Defaults to False.

        Returns:
            CompoundGroup: The restructured Compound with a new added group if at least one trajectory has satisfied the filter condition.
        """
        if flatten_trajectories:
            all_traj = self.collect_trajectories()

            grouped_traj = []
            ungrouped_traj = []

            if filter_func_trajectories is None:
                # Group all
                grouped_traj = all_traj
            else:
                for x in all_traj:
                    if (
                        isinstance(x, TrajectoryData)
                        and filter_func_trajectories(x.dataset)
                        or isinstance(x, TrajectoryGroup)
                        and filter_func_trajectories(x.get_group_info())
                    ):
                        grouped_traj.append(x)
                    else:
                        ungrouped_traj.append(x)
        else:
            grouped_traj: List[TrajectoryGroup | TrajectoryData] = []
            ungrouped_traj: List[TrajectoryGroup | TrajectoryData] = []

            if filter_func_trajectories is None:
                # Group all
                grouped_traj = [x.copy() for x in self.children.values()]  # type: ignore
            else:
                for x in self.children.values():
                    if (
                        isinstance(x, TrajectoryData)
                        and filter_func_trajectories(x.dataset)
                        or isinstance(x, TrajectoryGroup)
                        and filter_func_trajectories(x.get_group_info())
                    ):
                        grouped_traj.append(x.copy())
                    else:
                        ungrouped_traj.append(x.copy())

        group_children = dict(traj_list_to_child_mapping(grouped_traj))
        res_children = dict(traj_list_to_child_mapping(ungrouped_traj))
        new_group = TrajectoryGroup(GroupInfo(group_name, kwargs), group_children)
        res_children[group_name] = new_group

        res = type(self)(self.get_compound_info(), res_children)
        return res

    def map_over_trajectories(
        self,
        map_func: Callable[[Trajectory], T],
        result_as_dict=False,
        result_var_name: str = 'result',
    ) -> Self | dict:
        """Method to apply a function to all trajectories in this subtree.

        Args:
            map_func (Callable[[Trajectory], T]): Function to be applied to each individual trajectory in this database structure.
            result_as_dict (bool, optional): Whether to return the result as a dict or as a CompoundGroup structure. Defaults to False which yields a CompoundGroup.
            result_var_name (str,optional): The name of the result variable to be assigned in either the result dataset or in the result dict.

        Returns:
            CompoundGroup|dict: The result, either again as a CompoundGroup structure or as a layered dict structure.
        """
        if result_as_dict:
            res_dict = {
                k: v.map_over_trajectories(map_func, result_as_dict, result_var_name)
                for k, v in self.children.items()
            }
            # res_dict["_group_info"] = self.get_group_info()
            return res_dict
        else:
            return type(self)(
                self.get_compound_info(),
                {
                    k: v.map_over_trajectories(
                        map_func, result_as_dict, result_var_name
                    )
                    for k, v in self.children.items()
                },
            )  # type: ignore
