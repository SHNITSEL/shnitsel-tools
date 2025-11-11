from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Mapping, Self, TypeVar

from shnitsel.data.trajectory_format import Trajectory

from .db_trajectory_data import TrajectoryData
import xarray as xr

from .datatree_level import _datatree_level_attribute_key

T = TypeVar("T")


@dataclass
class GroupInfo:
    """Class to hold auxiliaryt info of a group of Trajectories in ShnitselDB"""

    group_name: str
    grouped_settings: Dict[str, Any] | None = None


class TrajectoryGroup(xr.DataTree):
    """DataTree node to keep track of a group of trajectories where properties defining the group can be set"""

    def __init__(
        self,
        group_info: GroupInfo | None,
        children: Mapping[str, TrajectoryData | Self] | None = None,
    ):
        super().__init__(None, children, group_info.group_name)

        self.attrs[_datatree_level_attribute_key] = "TrajectoryGroup"
        if group_info is not None:
            self.attrs["group_info"] = group_info.grouped_settings

    def get_group_info(self) -> GroupInfo:
        """Reconstruct the Group info object

        Returns:
            GroupInfo: The group information
        """
        if "group_info" in self.attrs:
            return GroupInfo(
                group_name=self.name, grouped_settings=self.attrs["group_info"]
            )
        else:
            return GroupInfo(group_name=self.name, grouped_settings={})

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
                k: v.map_over_trajectories(map_func, result_as_dict)
                for k, v in self.children.items()
            }
            # res_dict["_group_info"] = self.get_group_info()
            return res_dict
        else:
            return type(self)(
                self.get_group_info(),
                {
                    k: v.map_over_trajectories(map_func, result_as_dict)
                    for k, v in self.children.items()
                },
            )  # type: ignore
