from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Self, Sequence
from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.dataset_containers.frames import Frames
import xarray as xr


if TYPE_CHECKING:
    from .multi_layered import MultiSeriesLayered


@dataclass
class MultiSeriesStacked(Frames, MultiSeriesDataset):
    """A version of the multi-series dataset where the data is indexed along a sahred `frame` (Multi-index) dimension.
    There is no padding necessary to make the trajectories the same length.
    """

    _layered_repr_cached: "MultiSeriesLayered | None" = None

    def __init__(
        self, framesets: Sequence[Frames | Trajectory | xr.Dataset] | xr.Dataset
    ):
        if isinstance(framesets, xr.Dataset):
            assert 'frame' in framesets.coords, (
                "Stacked dataset must have `frame` dimension"
            )
            assert 'atrajectory' in framesets.coords, (
                "Stacked dataset must have `atrajectory` coordinate"
            )
            assert 'trajectory' in framesets.dims, (
                "Stacked dataset must have `trajectory` dimension"
            )

            MultiSeriesDataset.__init__(self, framesets)
            self._is_multi_trajectory = framesets.sizes['trajectory'] > 1
        else:
            from shnitsel.data.traj_combiner_methods import concat_trajs

            # TODO: FIXME: Stack frames into one single big frameset.
            is_multi_trajectory = False
            if len(framesets) > 1:
                is_multi_trajectory = True
            elif len(framesets) == 1 and framesets[0].is_multi_trajectory:
                is_multi_trajectory = True

            # TODO: FIXME: Make sure that concatenation would work. Convert variables to same unit, etc.

            # Build the concatenated trajectory. May trigger exceptions
            combined_dataset = concat_trajs(framesets)

            Frames.__init__(self, combined_dataset)
            MultiSeriesDataset.__init__(self, framesets, combined_dataset)
            self._is_multi_trajectory = is_multi_trajectory

    @property
    def grouping_dimension(self) -> str:
        return 'atrajectory'

    @cached_property
    def as_layered(self) -> "MultiSeriesLayered":
        """Get a layered representation of the stacked datasets in this object

        Returns
        -------
        MultiSeriesLayered
            The converted (or extracted from cache) layered version of this multi-data dataset.
        """
        from .multi_layered import MultiSeriesLayered
        from shnitsel.data.dataset_containers import wrap_dataset

        if self._layered_repr_cached is not None and isinstance(
            self._layered_repr_cached, MultiSeriesLayered
        ):
            if self._layered_repr_cached._stacked_repr_cached is not self:
                self._layered_repr_cached._stacked_repr_cached = self

            return self._layered_repr_cached

        if self._basis_data is not None:
            tmp_res = MultiSeriesLayered(self._basis_data)
        else:
            ds: xr.Dataset = self.dataset
            datasets: Sequence[Frames | Trajectory] = [
                wrap_dataset(
                    ds.sel(trajectory=id, atrajectory=id)
                    .drop_dims(['trajectory', 'atrajectory'], errors="ignore")
                    .drop_vars('atrajectory'),
                    expected_types=Trajectory | Frames,
                )
                for id in ds.coords['trajectory'].values
            ]

            tmp_res = MultiSeriesLayered(datasets)
        # Set self as cached result of the inverse conversion
        tmp_res._stacked_repr_cached = self
        self._layered_repr_cached = tmp_res
        return tmp_res

    @cached_property
    def as_stacked(self) -> Self:
        return self

    @classmethod
    def get_type_marker(cls) -> str:
        return "shnitsel::MultiSeriesStacked"
