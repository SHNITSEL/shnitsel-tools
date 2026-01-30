from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Literal, Self

from ..trajectory_grouping_params import TrajectoryGroupingMetadata


from .data_series import DataSeries
from .frames import Frames
import xarray as xr


@dataclass
class MetaInformation:
    """Meta information for trajectory setup"""

    input_format: Literal["sharc", "newtonx", "ase", "pyrai2md"] | None = None
    input_type: Literal['static', 'dynamic'] | None = None
    input_format_version: str | None = None
    theory_basis_set: str | None = None
    est_level: str | None = None


@dataclass
class Trajectory(DataSeries):
    # from .frames import Frames
    # from .per_state import PerState
    # from .inter_state import InterState

    def __init__(self, ds: xr.Dataset):
        assert 'frame' not in ds.dims, (
            "Dataset has `frame` dimension and cannot be considered a Trajectory"
        )
        assert "time" in ds.dims or 'time' in ds.coords, (
            "Dataset is missing `time` dimension and cannot be considered a Trajectory"
        )
        assert "atom" in ds.dims, (
            "Dataset is missing `atom` dimension and cannot be considered a Trajectory"
        )
        assert "state" in ds.dims, (
            "Dataset is missing `state` dimension and cannot be considered a Trajectory"
        )
        super().__init__(ds)
        self._is_multi_trajectory = False

    @cached_property
    def as_frames(self) -> "Frames":
        """Convert this trajectory to a frames version of this trajectory, where the leading dimension
        is `frame` instead of `time`.

        Returns
        --------
            Frames: The resulting frames instance with a stacked dimension `frame` and a new coordinate `active_trajectory` along the `frame` dimension
        """

        frame_ds = self.dataset.expand_dims(atrajectory=[self.trajectory_id]).stack(
            frame=["atrajectory", "time"]
        )
        return Frames(frame_ds)

    @cached_property
    def as_trajectory(self) -> Self:
        """Convert this trajectory to a trajectory.

        Returns
        --------
            Self: The same object that is already a trajectory
        """

        return self

    @property
    def is_multi_trajectory(self) -> bool:
        """Flag whether this is a multi-trajectory container.

        Overwritten by child classes that combine multiple trajectories into one object
        """
        return False

    @property
    def trajectory_input_path(self) -> str | None:
        """Input path from which the trajectory was loaded"""
        trajectory_input_path = self._param_from_vars_or_attrs('trajectory_input_path')
        return trajectory_input_path

    @property
    def leading_dim(self) -> str:
        """The leading dimension along which consistent configurations are indexed.
        Usually `time` or `frame`."""
        return "time"

    @classmethod
    def get_type_marker(cls) -> str:
        return "shnitsel::Trajectory"