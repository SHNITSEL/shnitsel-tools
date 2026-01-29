from types import UnionType
from typing import get_args, overload, TypeVar, TYPE_CHECKING

if TYPE_CHECKING:
    from shnitsel.data.tree.node import TreeNode

from .multi_layered import MultiSeriesLayered
from .multi_series import MultiSeriesDataset
from .multi_stacked import MultiSeriesStacked
from .data_series import DataSeries
from .shared import ShnitselDataset
from .trajectory import Trajectory
from .frames import Frames
from .inter_state import InterState
from .per_state import PerState
import xarray as xr

__all__ = [
    "MultiSeriesDataset",
    "MultiSeriesLayered",
    "MultiSeriesStacked",
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
    ds: xr.Dataset | Trajectory | Frames | DataSeries | ShnitselDataset,
    expected_types: type[ConvertedType],
) -> ConvertedType: ...


@overload
def wrap_dataset(
    ds: xr.Dataset | Trajectory | Frames | DataSeries | ShnitselDataset,
    expected_types: None = None,
) -> ShnitselDataset | xr.Dataset: ...


def wrap_dataset(
    ds: xr.Dataset | ShnitselDataset,
    expected_types: type[ConvertedType] | None = None,
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

    Notes
    -----
    This function can also be called with a tree structure as input and will automatically map itself over the leaves.
    This is only meant for internal Shnitsel tools use and may be removed at some point.

    """
    from shnitsel.data.tree import TreeNode

    # This is a special case we handle internally but do not wish to advertise
    if isinstance(ds, TreeNode):
        return ds.map_data(wrap_dataset, expected_types=expected_types)
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
            f"Could not convert input dataset to expected types {expected_types}.\n Input type was {type(raw_ds)}"
        )

    return ds
