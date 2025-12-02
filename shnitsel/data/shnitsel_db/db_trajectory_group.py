from dataclasses import dataclass
import logging
from typing import Any, Callable, Dict, List, Mapping, Self, TypeVar

from shnitsel.data.helpers import dataclass_from_dict

from .helpers import traj_list_to_child_mapping
from shnitsel.data.trajectory_format import Trajectory
from .db_trajectory_data import TrajectoryData
import xarray as xr

from .datatree_level import _datatree_level_attribute_key

T = TypeVar("T")


@dataclass
class GroupInfo:
    """Class to hold auxiliaryt info of a group of Trajectories in ShnitselDB"""

    group_name: str
    group_attributes: Dict[str, Any] | None = None


class TrajectoryGroup(xr.DataTree):
    """DataTree node to keep track of a group of trajectories where properties defining the group can be set"""

    def __init__(
        self,
        group_info: GroupInfo | None = None,
        children: Mapping[str, TrajectoryData | Self] | None = None,
    ):
        super().__init__(
            None, children, group_info.group_name if group_info is not None else None
        )

        self.attrs[_datatree_level_attribute_key] = "TrajectoryGroup"
        if group_info is not None:
            self.attrs["group_info"] = group_info.__dict__

    def is_level(self, target_level: str) -> bool:
        """Check whether we are at a certain level

        Args:
            target_level (str): Desired level to check for

        Returns:
            bool: True if this level satisfies the requirements
        """
        return target_level == "TrajectoryGroup"

    def get_group_info(self) -> GroupInfo:
        """Reconstruct the Group info object from settings stored in this node's attributes

        Returns:
            GroupInfo: The group information
        """
        if "group_info" in self.attrs:
            return dataclass_from_dict(GroupInfo, self.attrs["group_info"])
        else:
            return GroupInfo(
                group_name=self.name if self.name is not None else "unknown",
                group_attributes={},
            )

    def collect_trajectories(self) -> List[TrajectoryData]:
        """Function to retrieve all trajectories in this subtree

        Returns:
            List[TrajectoryData]: List of all nodes with TrajectoryData type
        """
        res = []

        for x in self.children.values():
            res += x.collect_trajectories()

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
            result_as_dict (bool, optional): Whether to return the result as a dict or as a TrajectoryGroup structure. Defaults to False which yields a CompoundGroup.
            result_var_name (str,optional): The name of the result variable to be assigned in either the result dataset or in the result dict.

        Returns:
            TrajectoryGroup|dict: The result, either again as a TrajectoryGroup structure or as a layered dict structure.
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
                self.get_group_info(),
                {
                    k: v.map_over_trajectories(
                        map_func, result_as_dict, result_var_name
                    )
                    for k, v in self.children.items()
                },
            )  # type: ignore

    def merge_with(self, other: Self) -> Self:
        """Function to merge two TrajectoryGroups into one.

        Called when merging two database states. Will fail if group_info differs between compounds to avoid loss of information.

        Args:
            other (CompoundGroup): The other CompoundGroup to be merged

        Raises:
            ValueError: Raised if the compound_info differs.

        Returns:
            CompoundGroup: A CompoundGroup object holding the entire merged subtree
        """
        own_info = self.get_group_info()
        other_info = other.get_group_info()
        if other_info != own_info:
            message = f"Cannot merge groups with conflicting group information: {other_info} vs. {own_info}"
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

        for k, v in self.children.items():
            key_str = str(k)
            if key_str in self.children:
                continue
            else:
                res_children.append(v.copy())

        return type(self)(own_info, traj_list_to_child_mapping(res_children))  # type: ignore
