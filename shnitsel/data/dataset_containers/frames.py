from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
import xarray as xr
from .trajectory import Trajectory


@dataclass
class Frames(Trajectory):
    # TODO: FIXME: This should not be a subclass of Trajectory. It should have similar accessors, but probably from a different base class.
    _is_multi_trajectory: bool = False

    def __init__(self, ds: xr.Dataset):
        assert "time" not in ds.dims, (
            "Dataset has `time` dimension and cannot be considered a set of Frames"
        )
        assert "frame" in ds.dims, (
            "Dataset is missing `frame` dimension and cannot be considered a set of Frames"
        )
        assert "atom" in ds.dims, (
            "Dataset is missing `atom` dimension and cannot be considered a set of Frames"
        )
        assert "state" in ds.dims, (
            "Dataset is missing `state` dimension and cannot be considered a set of Frames"
        )
        super().__init__(ds)

        # TODO: FIXME: This should be harmonized across all creation and use points. Make the frame-component `active_trajectory` and the per-trajectory property `trajectory`
        if "trajectory" in ds.dims:
            # Check if we have a dimension to select properties of different trajectories.
            self._is_multi_trajectory = True

    @property
    def is_multi_trajectory(self) -> bool:
        return self._is_multi_trajectory


class MultiTrajectoryFrames(Frames):
    def __init__(self, framesets: Iterable[Frames]):
        # TODO: FIXME: Concatenate frames into one single big frameset.
        ...
        super().__init__(framesets)
        self._is_multi_trajectory = True
