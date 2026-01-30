from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Self, Sequence

from shnitsel.data.dataset_containers.data_series import DataSeries
from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.dataset_containers.frames import Frames
import xarray as xr


if TYPE_CHECKING:
    from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked


@dataclass
class MultiSeriesLayered(MultiSeriesDataset):
    """A version of the multi-series dataset where the data is indexed along a new `trajectory` dimension.
    Missing data across trajectories is padded with `np.nan` values and can thus lead to typing issues.
    """

    _stacked_repr_cached: "MultiSeriesStacked | None" = None

    def __init__(
        self, framesets: xr.Dataset | Sequence[Frames | Trajectory | xr.Dataset]
    ):
        if isinstance(framesets, xr.Dataset):
            assert 'trajectory' in framesets.dims, (
                "Layered dataset must have `trajectory` dimension"
            )
            assert 'atrajectory' not in framesets.coords, (
                "Layered dataset must not have `atrajectory` coordinate"
            )
            assert 'time' in framesets.dims or 'frame' in framesets.dims, (
                "Layered dataset must have `time` or 'frame' dimension"
            )
            if 'frame' in framesets.coords and 'frame' in framesets.indexes:
                assert 'trajectory' not in framesets.indexes['frame'].names, (
                    "Layered dataset must not have `frames` dimension as a multi-index on `trajectory` dimension"
                )

            MultiSeriesDataset.__init__(self, framesets)
            self._is_multi_trajectory = framesets.sizes['trajectory'] > 1

        else:
            from shnitsel.data.traj_combiner_methods import layer_trajs

            # TODO: FIXME: Layer frames into one single big frameset.
            self._frame_sources = framesets
            is_multi_trajectory = False
            if len(framesets) > 1:
                is_multi_trajectory = True
            elif len(framesets) == 1 and framesets[0].is_multi_trajectory:
                is_multi_trajectory = True

            # TODO: FIXME: Make sure that concatenation would work. Convert variables to same unit, etc.

            # Build the concatenated trajectory. May trigger exceptions
            combined_dataset = layer_trajs(framesets)

            DataSeries.__init__(self, combined_dataset)
            MultiSeriesDataset.__init__(self, framesets, combined_dataset)
            self._is_multi_trajectory = is_multi_trajectory

    @property
    def grouping_dimension(self) -> str:
        return 'trajectory'

    @cached_property
    def as_stacked(self) -> "MultiSeriesStacked":
        """Get a stacked representation of the layered datasets in this object

        Returns
        -------
        MultiSeriesStacked
            The converted (or extracted from cache) stacked version of this multi-data dataset.
        """
        from .multi_stacked import MultiSeriesStacked
        from shnitsel.data.dataset_containers import wrap_dataset

        if self._stacked_repr_cached is not None and isinstance(
            self._stacked_repr_cached, MultiSeriesStacked
        ):
            if self._stacked_repr_cached._layered_repr_cached is not self:
                self._stacked_repr_cached._layered_repr_cached = self

            return self._stacked_repr_cached

        if self._basis_data is not None:
            tmp_res = MultiSeriesStacked(self._basis_data)
        else:
            # TODO: FIXME: Fix the typing issue with `wrap_dataset`
            ds: xr.Dataset = self.dataset
            datasets: Sequence[Frames | Trajectory] = [
                wrap_dataset(
                    (tmp_ds := ds.sel(trajectory=id))
                    .isel(time=slice(None, tmp_ds.max_ts.item() + 1))
                    .drop_dims('trajectory', errors="ignore"),
                    expected_types=Trajectory | Frames,
                )
                for id in ds.coords['trajectory'].values
            ]

            tmp_res = MultiSeriesStacked(datasets)

        # Set self as cached result of the inverse conversion
        tmp_res._layered_repr_cached = self
        self._stacked_repr_cached = tmp_res
        return tmp_res

    @cached_property
    def as_layered(self) -> Self:
        return self

    @classmethod
    def get_type_marker(cls) -> str:
        return "shnitsel::MultiSeriesLayered"
