from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import xarray as xr
from .data_series import DataSeries


@dataclass
class Frames(DataSeries):
    # TODO: FIXME: This should not be a subclass of Trajectory. It should have similar accessors, but probably from a different base class.
    _is_multi_trajectory: bool = False

    def __init__(self, ds: xr.Dataset):
        assert "time" not in ds.dims, (
            "Dataset has `time` dimension and cannot be considered a set of Frames"
        )
        assert "frame" in ds.dims or 'frame' in ds.coords, (
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
        if (
            "trajectory" in ds.dims
            and ds.sizes["trajectory"] > 1
            or "atrajectory" in ds.coords
            and len(set(ds.coords["atrajectory"].values)) > 1
        ):
            # Check if we have a dimension to select properties of different trajectories.
            self._is_multi_trajectory = True

    @property
    def leading_dim(self) -> str:
        """The leading dimension along which consistent configurations are indexed.
        Usually `time` or `frame`."""
        return "frame"

    @property
    def trajid(self) -> int | str | None:
        """Id of the trajectory. If assigned it is expected to be unique across the same input
        but may clash with other trajectory ids if multiple separate imports are combined
        or indepdendent simulation data is combined."""
        trajid = super().trajid
        if trajid is None:
            # Try and get the own trajectory id from the active trajectory
            tmp_atraj = self._param_from_vars_or_attrs('atrajectory')
            if tmp_atraj is not None:
                trajids = set(
                    tmp_atraj.values
                    if isinstance(tmp_atraj, xr.DataArray)
                    else tmp_atraj
                )
                if len(trajids) == 1:
                    return trajids.pop()
        return trajid

    @property
    def atrajectory(self) -> xr.DataArray | None:
        """Ids of the active trajectory in this frameset if present"""
        return self.coords['atrajectory']

    @property
    def active_trajectory(self) -> xr.DataArray | None:
        """Ids of the active trajectory in this frameset if present"""
        return self.atrajectory

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
