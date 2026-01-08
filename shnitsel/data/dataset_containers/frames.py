from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import xarray as xr
from .trajectory import Trajectory


@dataclass
class Frames(Trajectory):
    # TODO: FIXME: This should not be a subclass of Trajectory. It should have similar accessors, but probably from a different base class.
    _is_multi_trajectory: bool = False

    def __init__(self, ds: xr.Dataset):
        assert (
            "time" not in ds.dims
        ), "Dataset has `time` dimension and cannot be considered a set of Frames"
        assert (
            "frame" in ds.dims
        ), "Dataset is missing `frame` dimension and cannot be considered a set of Frames"
        assert (
            "atom" in ds.dims
        ), "Dataset is missing `atom` dimension and cannot be considered a set of Frames"
        assert (
            "state" in ds.dims
        ), "Dataset is missing `state` dimension and cannot be considered a set of Frames"
        super().__init__(ds)

        # TODO: FIXME: This should be harmonized across all creation and use points. Make the frame-component `active_trajectory` and the per-trajectory property `trajectory`
        if "trajectory" in ds.dims:
            # Check if we have a dimension to select properties of different trajectories.
            self._is_multi_trajectory = True

    @property
    def is_multi_trajectory(self) -> bool:
        return self._is_multi_trajectory


@dataclass
class MultiTrajectoryFrames(Frames):
    _frame_sources: Sequence[Frames] | None = None

    def __init__(self, framesets: Sequence[Frames]):
        # TODO: FIXME: Concatenate frames into one single big frameset.
        self._frame_sources = framesets
        is_multi_trajectory = False
        if len(framesets) > 1:
            is_multi_trajectory = True
        elif len(framesets) == 1 and framesets[0].is_multi_trajectory:
            is_multi_trajectory = True

        # TODO: FIXME: Make sure that concatenation would work. Convert variables to same unit, etc.

        # Build the concatenated trajectory. May trigger exceptions
        combined_dataset = xr.concat(
            [frames.dataset for frames in framesets],
            dim='frame',
        )

        super().__init__(combined_dataset)
        self._is_multi_trajectory = is_multi_trajectory
