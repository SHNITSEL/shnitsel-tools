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

    @property
    def leading_dimension(self) -> str:
        if "frame" in self.dataset.dims:
            return "frame"
        elif "time" in self.dataset.dims:
            return "time"
        else:
            raise ValueError(
                "Unknown leading dimension of the contained dataset. The Dataset may have been misconstructed or loaded from malformed data."
            )

    @property
    def state_names(self):
        pass

    @property
    def state_types(self):
        pass

    @property
    def state_magnetic_number(self):
        pass

    @property
    def state_degeneracy_group(self):
        pass

    @property
    def state_charges(self):
        pass

    @property
    def active_state(self):
        pass

    @property
    def state_diagonal(self):
        pass

    @property
    def atom_names(self):
        pass

    @property
    def atom_numbers(self):
        pass

    # TODO: Forward all unmet requests to dataset.


@dataclass
class ShnitselDerivedDataset(ShnitselDataset):
    _base_dataset: xr.Dataset | None

    def __init__(self, base_ds: xr.Dataset | None, derived_ds: xr.Dataset):
        super().__init__(derived_ds)
        self._base_dataset = base_ds

    @property
    def base(self) -> xr.Dataset | None:
        return self.base

    # TODO: Forward all unmet requests to dataset.
