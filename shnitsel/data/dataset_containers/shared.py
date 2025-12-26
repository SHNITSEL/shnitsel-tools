from dataclasses import dataclass
import xarray as xr


@dataclass
class ShnitselDataset:
    _raw_dataset: xr.Dataset

    def __init__(self, ds: xr.Dataset):
        self._raw_dataset = ds

    @property
    def dataset(self) -> xr.Dataset:
        return self._raw_dataset

    # TODO: Forward all unmet requests to dataset.


@dataclass
class ShnitselDerivedDataset:
    _base_dataset: xr.Dataset | None
    _raw_dataset: xr.Dataset

    def __init__(self, base_ds: xr.Dataset | None, derived_ds: xr.Dataset):
        self._raw_dataset = derived_ds
        self._base_dataset = base_ds

    @property
    def dataset(self) -> xr.Dataset:
        return self._raw_dataset

    @property
    def base(self) -> xr.Dataset | None:
        return self.base

    # TODO: Forward all unmet requests to dataset.
