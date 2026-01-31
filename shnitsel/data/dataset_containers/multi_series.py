from dataclasses import dataclass
from functools import cached_property
from shnitsel.data.dataset_containers.data_series import DataSeries
import xarray as xr
from typing import Sequence, TYPE_CHECKING

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.trajectory_grouping_params import TrajectoryGroupingMetadata

if TYPE_CHECKING:
    from .multi_layered import MultiSeriesLayered
    from .multi_stacked import MultiSeriesStacked


@dataclass
class MultiSeriesDataset(DataSeries):
    """Class to serve as the basis for Layered and Stacked multi-dataseries datasets.

    Is itself a DataSeries, but with different, more specific semantics than a generic DataSeries.

    """

    _basis_data: Sequence[Frames | Trajectory | xr.Dataset] | None = None

    def __init__(
        self,
        basis: xr.Dataset | Sequence[Frames | Trajectory | xr.Dataset],
        combined: xr.Dataset | None = None,
    ):
        combined_data: xr.Dataset
        if isinstance(basis, xr.Dataset):
            assert (
                'trajectory' in basis.dims
                or 'atrajectory' in basis.coords
                and 'frame' in basis.dims
            ), (
                "Dataset used as basis for a MultiSeries Dataset, needs either a `trajectory` or `frame` dimension with the latter accompanied by a `atrajectory` coordinate to select the individual trajectory."
            )
            self._basis_data = None
            combined_data = basis
        elif isinstance(basis, Sequence):
            self._basis_data = basis
            assert combined is not None, (
                "If a multi-dataseries basis is passed, the child class or caller needs to implement the combination logic and provide the `combined` argument."
            )
            combined_data = combined
        else:
            raise ValueError(
                f"Unsupported basis type for MultiSeriesDataset: {type(basis)}"
            )

        super().__init__(combined_data)

    @property
    def grouping_dimension(self) -> str:
        raise NotImplemented

    @cached_property
    def as_stacked(self) -> "MultiSeriesStacked":
        raise NotImplementedError(
            "Subclasses are required to implement the `as_stacked` property of Multi-Series datasets"
        )

    @cached_property
    def as_layered(self) -> "MultiSeriesLayered":
        raise NotImplementedError(
            "Subclasses are required to implement the `as_layered` property of Multi-Series datasets"
        )

    def get_grouping_metadata(self) -> TrajectoryGroupingMetadata:
        def to_single_val(entries):
            if isinstance(entries, xr.DataArray):
                if entries.size > 1:
                    return entries.values.flatten()[0]
                return entries.item()

            return entries

        res = TrajectoryGroupingMetadata(
            delta_t_in_fs=to_single_val(self.delta_t),
            input_format_name=to_single_val(self.input_format),
            input_format_version=to_single_val(self.input_format_version),
            est_level=to_single_val(self.est_level),
            theory_basis_set=to_single_val(self.theory_basis_set),
            charge_in_e=to_single_val(self.charge),
            # TODO: FIXME: We should differentiate by all state attributes.
            num_states=len(self.state_ids),
        )
        return res

    @classmethod
    def get_type_marker(cls) -> str:
        return "shnitsel::MultiSeriesDataset"
