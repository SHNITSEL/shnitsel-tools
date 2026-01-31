from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
from typing import Self, TYPE_CHECKING
import xarray as xr

from .data_series import DataSeries

if TYPE_CHECKING:
    from .trajectory import Trajectory


@dataclass
class Frames(DataSeries):
    _is_multi_trajectory: bool = False

    def __new__(cls, ds: xr.Dataset):
        from .multi_stacked import MultiSeriesStacked

        if cls == Frames:
            if (
                "trajectory" in ds.dims
                and ds.sizes["trajectory"] > 1
                or "atrajectory" in ds.coords
                and len(set(ds.coords["atrajectory"].values)) > 1
            ):
                return object.__new__(MultiSeriesStacked)
        return object.__new__(cls)

    def __init__(self, ds: xr.Dataset):
        assert "time" not in ds.dims, (
            "Dataset has `time` dimension and cannot be considered a set of Frames"
        )
        if "frame" not in ds.dims and "frame" in ds.coords:
            curr_fdim = ds.coords["frame"].dims[0]
            ds = ds.swap_dims({curr_fdim: "frame"})
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
        if (
            "trajectory" in ds.dims
            and ds.sizes["trajectory"] > 1
            or "atrajectory" in ds.coords
            and len(set(ds.coords["atrajectory"].values)) > 1
        ):
            # Check if we have a dimension to select properties of different trajectories.
            self._is_multi_trajectory = True

    @cached_property
    def as_frames(self) -> Self:
        """Idempotent conversion to Frame instance"""
        return self

    @cached_property
    def as_trajectory(self) -> "Trajectory":
        """Attempt to convert this dataset into a Trajectory instance.

        Drops the `atrajectory` and `trajectory` dimensions of the Frameset and replaces the
        `frame` dimension with a `time` dimension before conversion.

        Returns
        -------
        Trajectory
            The converted dataset underlying this Frameset.
        """
        from .trajectory import Trajectory

        assert not self.is_multi_trajectory, (
            "Cannot convert multi-trajectory frames to trajectory format"
        )
        assert 'time' in self.coords, (
            "Frameset is lacking `time` coordinate for conversion to trajectory"
        )

        ds = self.dataset
        trajid = self.trajectory_id

        if (
            "frame" in ds.indexes
            and ds.indexes["frame"].names is not len(ds.indexes["frame"].names) > 1
        ):
            processed_ds = ds.unstack('frame')
        else:
            processed_ds = ds

        processed_ds = (
            processed_ds.drop_vars('atrajectory', errors='ignore')
            .swap_dims({'frame': 'time'})
            .isel(trajectory=0, missing_dims="ignore")
        )

        if trajid is not None:
            trajid_coord = xr.DataArray(
                trajid,
                dims=[],
                attrs={'long_name': "The id of the indexed/current trajectory"},
            )
            processed_ds = processed_ds.assign_coords(trajectory=trajid_coord)

        return Trajectory(processed_ds)

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
            trajid = self._param_from_vars_or_attrs('atrajectory')

        if trajid is not None and not isinstance(trajid, (int, str)):
            trajids = set(
                trajid.values[:].tolist()
                if isinstance(trajid, xr.DataArray)
                else trajid
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

    @classmethod
    def get_type_marker(cls) -> str:
        return "shnitsel::Frames"
