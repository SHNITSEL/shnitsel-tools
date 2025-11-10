from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Self

from .db_trajectory_data import TrajectoryData
import xarray as xr

from .datatree_level import _datatree_level_attribute_key


@dataclass
class GroupInfo:
    """Class to hold auxiliaryt info of a group of Trajectories in ShnitselDB"""

    grouped_settings: Dict[str, Any]


class TrajectoryGroup(xr.DataTree):
    """DataTree node to keep track of a group of trajectories where properties defining the group can be set"""

    def __init__(
        self,
        group_info: GroupInfo | None,
        children: Mapping[str, TrajectoryData | Self] | None = None,
        name: str | None = None,
    ):
        super().__init__(None, children, name)

        self.attrs[_datatree_level_attribute_key] = "TrajectoryGroup"
        if group_info is not None:
            self.attrs["group_info"] = group_info

    def collect_trajectories(self) -> List[TrajectoryData]:
        """Function to retrieve all trajectories in this subtree

        Returns:
            List[TrajectoryData]: List of all nodes with TrajectoryData type
        """
        res = []

        for x in self.children.values():
            res += x.collect_trajectories()

        return res
