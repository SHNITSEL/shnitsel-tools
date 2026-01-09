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
    def state_ids(self):
        if "state" not in self.dataset.coords:
            raise KeyError("No coordinate `state` provided for the trajectory")
        return self.dataset.coords["state"]

    @property
    def state_names(self):
        if "state_names" not in self.dataset.coords:
            raise KeyError("No coordinate `state_names` provided for the trajectory")
        return self.dataset.coords["state_names"]

    @property
    def state_types(self):
        if "state_types" not in self.dataset.coords:
            raise KeyError("No coordinate `state_types` provided for the trajectory")
        return self.dataset.coords["state_types"]

    @property
    def state_magnetic_number(self):
        if "state_magnetic_number" not in self.dataset.coords:
            raise KeyError(
                "No coordinate `state_magnetic_number` provided for the trajectory"
            )
        return self.dataset.coords["state_magnetic_number"]

    @property
    def state_degeneracy_group(self):
        if "state_degeneracy_group" not in self.dataset.coords:
            raise KeyError(
                "No coordinate `state_degeneracy_group` provided for the trajectory"
            )
        return self.dataset.coords["state_degeneracy_group"]

    @property
    def state_charges(self):
        if "state_charges" not in self.dataset.coords:
            raise KeyError("No coordinate `state_charges` provided for the trajectory")
        return self.dataset.coords["state_charges"]

    @property
    def active_state(self):
        if "astate" not in self.dataset.coords:
            if "astate" not in self.dataset:
                raise KeyError(
                    "No coordinate `astate` holding the active state id provided for the trajectory"
                )
            return self.dataset['astate']
        return self.dataset.coords["astate"]

    @property
    def state_diagonal(self):
        if "sdiag" not in self.dataset.coords:
            raise KeyError(
                "No coordinate `sdiag` holding the active state id provided for the trajectory"
            )
        return self.dataset.coords["sdiag"]

    @property
    def atom_names(self):
        if "atom_names" not in self.dataset.coords:
            raise KeyError("No coordinate `atom_names` provided for the trajectory")
        return self.dataset.coords["atom_names"]

    @property
    def atom_numbers(self):
        if "atom_numbers" not in self.dataset.coords:
            raise KeyError("No coordinate `atom_numbers` provided for the trajectory")
        return self.dataset.coords["atom_numbers"]

    # TODO: Forward all unmet requests to dataset.

    @property
    def dims(self):
        return self.dataset.dims

    @property
    def coords(self):
        return self.dataset.coords

    @property
    def sizes(self):
        return self.dataset.sizes

    @property
    def data_vars(self):
        return self.dataset.data_vars

    def has_variable(self, name: str) -> bool:
        return name in self.data_vars

    def has_dimension(self, name: str) -> bool:
        return name in self.dims

    def has_coordinate(self, name: str) -> bool:
        return name in self.coords


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
