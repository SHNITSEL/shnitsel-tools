from typing import List, Self
import xarray as xr

from shnitsel.data.trajectory_format import Trajectory
from .datatree_level import _datatree_level_attribute_key

class TrajectoryData(xr.DataTree):
    """DataTree node to keep track of a single trajectory entry"""

    def __init__(
        self,
        dataset: xr.Dataset | Trajectory | None = None,
        name: str | None = None,
    ):
        super().__init__(dataset=dataset, children=None, name=name)

        self.attrs[_datatree_level_attribute_key] = "TrajectoryData"

    def collect_trajectories(self) -> List[Self]:
        """Function to retrieve all trajectories in this subtree

        Returns:
            List[TrajectoryData]: List of all nodes with TrajectoryData type
        """
        return [self.copy()]
