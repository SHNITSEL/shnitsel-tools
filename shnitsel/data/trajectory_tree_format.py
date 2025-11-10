from dataclasses import dataclass
from typing import Any, Dict, Mapping
import xarray as xr

from shnitsel.data.trajectory_format import Trajectory


@dataclass
class GroupInfo:
    grouped_settings: Dict[str, Any]


@dataclass
class CompoundInfo:
    compound_name: str = "unknown"
    smile_repr: str | None = None


class TrajectoryData(xr.DataTree):
    """DataTree node to keep track of a single trajectory entry"""

    def __init__(
        self,
        dataset: xr.Dataset | Trajectory | None = None,
        name: str | None = None,
    ):
        super().__init__(dataset=dataset, children=None, name=name)

        self.attrs["DataTree_Level"] = "TrajectoryData"


class TrajectoryGroup(xr.DataTree):
    """DataTree node to keep track of a group of trajectories where properties defining the group can be set"""

    def __init__(
        self,
        group_info: GroupInfo,
        children: Mapping[str, TrajectoryData] | None = None,
        name: str | None = None,
    ):
        super().__init__(None, children, name)

        self.attrs["DataTree_Level"] = "TrajectoryGroup"
        self.attrs["group_info"] = group_info


class CompoundGroup(xr.DataTree):
    """DataTree node to keep track of all data associated with a common compound within the datatree"""

    def __init__(
        self,
        coumpound_info: CompoundInfo,
        children: Mapping[
            str,
            TrajectoryGroup | TrajectoryData,
        ]
        | None = None,
    ):
        super().__init__(None, children, coumpound_info.compound_name)
        self.attrs["DataTree_Level"] = "CompoundGroup"
        self.attrs["compound_info"] = coumpound_info


class ShnitselDBRoot(xr.DataTree):
    """DataTree root node to keep track of Shnitsel data in a predefined structural hierarchy"""

    def __init__(self, compounds: Mapping[str, CompoundGroup] | None = None):
        super().__init__(dataset=None, children=compounds, name="ROOT")
        self.attrs["DataTree_Level"] = "ShnitselDBRoot"
