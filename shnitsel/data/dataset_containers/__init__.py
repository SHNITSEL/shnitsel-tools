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


def wrap_dataset(
    ds: xr.Dataset | Trajectory | Frames | DataSeries | ShnitselDataset,
) -> Trajectory | Frames | DataSeries | ShnitselDataset | xr.Dataset:
    """Helper function to wrap a generic xarray dataset in a wrapper container

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to wrap

    Returns
    -------
    Trajectory | Frames | DataSeries | ShnitselDataset | xr.Dataset
        The wrapped dataset or the original dataset if no conversion was possible
    """
    if isinstance(ds, (Trajectory, Frames, DataSeries, ShnitselDataset)):
        return ds

    try:
        return Trajectory(ds)
    except:
        try:
            return Frames(ds)
        except:
            try:
                return DataSeries(ds)
            except:
                try:
                    return ShnitselDataset(ds)
                except:
                    return ds
