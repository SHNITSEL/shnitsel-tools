from types import UnionType
from typing import get_args, overload, TypeVar

from shnitsel.data.dataset_containers.multi_layered import MultiSeriesLayered
from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked
from .data_series import DataSeries
from .shared import ShnitselDataset
from .trajectory import Trajectory
from .frames import Frames
from .inter_state import InterState
from .per_state import PerState
import xarray as xr

__all__ = [
    "Trajectory",
    "Frames",
    "InterState",
    "PerState",
    "DataSeries",
    "ShnitselDataset",
    "wrap_dataset",
]

ConvertedType = TypeVar(
    "ConvertedType",
    bound=(ShnitselDataset | UnionType),
)


@overload
def wrap_dataset(
    ds: xr.Dataset | ShnitselDataset,
    expected_types: type[ConvertedType],
) -> ConvertedType: ...


@overload
def wrap_dataset(
    ds: xr.Dataset | Trajectory | Frames | DataSeries | ShnitselDataset,
    expected_types: None = None,
) -> ShnitselDataset | xr.Dataset: ...


def wrap_dataset(
    ds: xr.Dataset | ShnitselDataset,
    expected_types: type[ConvertedType] | UnionType | None = None,
) -> ConvertedType | ShnitselDataset | xr.Dataset:
    """Helper function to wrap a generic xarray dataset in a wrapper container

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to wrap or an already wrapped dataset that may not need conversion.
    expected_types: type[ConvertedType] | UnionType, optional
        Can be used to limit which wrapped format would be acceptable as a result.
        If set, an assertion error will be triggered if the `ds` parameter could not be wrapped in the appropriate
        type.

    Returns
    -------
    ConvertedType | ShnitselDataset | xr.Dataset
        The wrapped dataset or the original dataset if no conversion was possible
    """
    if expected_types is not None:
        accepted_types = (
            get_args(expected_types)
            if isinstance(expected_types, UnionType)
            else [expected_types]
        )
        if isinstance(ds, expected_types):
            return ds
    else:
        if isinstance(ds, ShnitselDataset):
            return ds
        accepted_types = [
            MultiSeriesLayered,
            MultiSeriesStacked,
            Trajectory,
            Frames,
            DataSeries,
            ShnitselDataset,
        ]

    raw_ds: xr.Dataset
    if isinstance(ds, ShnitselDataset):
        raw_ds = ds.dataset
    else:
        raw_ds = ds

    if MultiSeriesLayered in accepted_types or MultiSeriesDataset in accepted_types:
        try:
            return MultiSeriesLayered(raw_ds)
        except:
            pass
    if MultiSeriesStacked in accepted_types or MultiSeriesDataset in accepted_types:
        try:
            return MultiSeriesStacked(raw_ds)
        except:
            pass
    if (
        Trajectory in accepted_types
        or DataSeries in accepted_types
        or ShnitselDataset in accepted_types
    ):
        try:
            return Trajectory(raw_ds)
        except:
            pass
    if (
        Frames in accepted_types
        or DataSeries in accepted_types
        or ShnitselDataset in accepted_types
    ):
        try:
            return Frames(raw_ds)
        except:
            pass

    if DataSeries in accepted_types or ShnitselDataset in accepted_types:
        try:
            return DataSeries(raw_ds)
        except:
            pass

    if ShnitselDataset in accepted_types:
        try:
            return ShnitselDataset(raw_ds)
        except:
            pass
    if expected_types is not None:
        raise AssertionError(
            f"Could not convert input dataset to expected types {expected_types}"
        )

    return ds
