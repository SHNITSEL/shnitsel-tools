from dataclasses import dataclass
from .shared import ShnitselDataset
import xarray as xr


@dataclass
class Trajectory(ShnitselDataset):
    def __init__(self, ds: xr.Dataset):
        assert (
            "time" in ds.dims
        ), "Dataset is missing `time` dimension and cannot be considered a Trajectory"
        assert (
            "atom" in ds.dims
        ), "Dataset is missing `atom` dimension and cannot be considered a Trajectory"
        assert (
            "state" in ds.dims
        ), "Dataset is missing `state` dimension and cannot be considered a Trajectory"
        super().__init__(ds)
