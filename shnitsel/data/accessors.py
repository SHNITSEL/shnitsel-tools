from dataclasses import dataclass
import xarray as xr


@dataclass
class ShnitselDataset:
    _raw_dataset: xr.Dataset

    def __init__(self, ds: xr.Dataset):
        self._raw_dataset = ds

    def dataset(self) -> xr.Dataset:
        return self._raw_dataset

    # TODO: Forward all unmet requests to dataset.


@dataclass
class Trajectory(ShnitselDataset):
    def __init__(self, ds: xr.Dataset):
        assert 'time' in ds.dims, (
            'Dataset is missing `time` dimension and cannot be considered a Trajectory'
        )
        assert 'atom' in ds.dims, (
            'Dataset is missing `atom` dimension and cannot be considered a Trajectory'
        )
        assert 'state' in ds.dims, (
            'Dataset is missing `state` dimension and cannot be considered a Trajectory'
        )
        super().__init__(ds)


@dataclass
class Frames(ShnitselDataset):
    def __init__(self, ds: xr.Dataset):
        assert 'time' not in ds.dims, (
            'Dataset has `time` dimension and cannot be considered a set of Frames'
        )
        assert 'frame' in ds.dims, (
            'Dataset is missing `frame` dimension and cannot be considered a set of Frames'
        )
        assert 'atom' in ds.dims, (
            'Dataset is missing `atom` dimension and cannot be considered a set of Frames'
        )
        assert 'state' in ds.dims, (
            'Dataset is missing `state` dimension and cannot be considered a set of Frames'
        )
        super().__init__(ds)


@dataclass
class InterState(ShnitselDataset):
    def __init__(self, ds: xr.Dataset):
        assert 'state' in ds.dims, (
            'Dataset is missing `state` dimension and cannot be considered an InterState set of variables.'
        )
        # TODO: FIXME: Calculate inter-state variables and cache in original dataset
        super().__init__(ds)


@dataclass
class PerState(ShnitselDataset):
    def __init__(self, ds: xr.Dataset):
        assert 'state' in ds.dims, (
            'Dataset is missing `state` dimension and cannot be considered an PerState set of variables.'
        )
        # TODO: FIXME: Calculate per-state variables and cache in original dataset
        super().__init__(ds)
