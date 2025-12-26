from dataclasses import dataclass
import xarray as xr
from .shared import ShnitselDataset


@dataclass
class Frames(ShnitselDataset):
    def __init__(self, ds: xr.Dataset):
        assert (
            "time" not in ds.dims
        ), "Dataset has `time` dimension and cannot be considered a set of Frames"
        assert (
            "frame" in ds.dims
        ), "Dataset is missing `frame` dimension and cannot be considered a set of Frames"
        assert (
            "atom" in ds.dims
        ), "Dataset is missing `atom` dimension and cannot be considered a set of Frames"
        assert (
            "state" in ds.dims
        ), "Dataset is missing `state` dimension and cannot be considered a set of Frames"
        super().__init__(ds)


@dataclass
class Frames(ShnitselDataset):
    def __init__(self, ds: xr.Dataset):
        assert (
            "time" not in ds.dims
        ), "Dataset has `time` dimension and cannot be considered a set of Frames"
        assert (
            "frame" not in ds.dims
        ), "Dataset is missing `frame` dimension and cannot be considered a set of Frames"
        assert (
            "atom" in ds.dims
        ), "Dataset is missing `atom` dimension and cannot be considered a set of Frames"
        assert (
            "state" in ds.dims
        ), "Dataset is missing `state` dimension and cannot be considered a set of Frames"
        super().__init__(ds)
