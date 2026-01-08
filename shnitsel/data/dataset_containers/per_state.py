from dataclasses import dataclass
from typing import Literal

from .trajectory import Trajectory
from .shared import ShnitselDerivedDataset
from .frames import Frames

import xarray as xr


@dataclass
class PerState(ShnitselDerivedDataset):
    _original_frames: Frames | Trajectory | None

    def __init__(
        self,
        frames: Frames | Trajectory | None = None,
        /,
        direct_perstate_data: xr.Dataset | None = None,
    ):
        from shnitsel.analyze.stats import get_per_state
        # TODO: FIXME: Calculate per-state variables and cache in original dataset

        base_ds = None
        if frames is not None:
            assert (
                "state" in frames.dataset.dims
            ), "Dataset is missing `state` dimension and cannot be considered an PerState set of variables."
            # TODO: FIXME: Calculate per-state variables and cache in original dataset
            base_ds = frames.dataset

        if direct_perstate_data is not None:
            per_state_props = direct_perstate_data
        elif base_ds is not None:
            per_state_props = get_per_state(base_ds)
        else:
            per_state_props = xr.Dataset()

        self._original_frames = frames

        super().__init__(base_ds, per_state_props)

    @property
    def energy(self) -> xr.DataArray:
        if "energy" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `energy` to encode per-state energy in trajectory data"
            )
        return self.dataset.data_vars["energy"]

    @property
    def dipole_permanent(self) -> xr.DataArray:
        if "dip_perm" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `dip_perm` to encode per-state permanent dipole moments in trajectory data"
            )
        return self.dataset.data_vars["dip_perm"]

    @property
    def dipole_permanent_norm(self) -> xr.DataArray:
        if "dip_perm_norm" not in self.dataset.data_vars:
            if 'dip_perm' not in self.dataset.data_vars:
                raise KeyError(
                    "No variable `dip_perm_norm` to encode per-state permanent dipole moments in trajectory data"
                )
            from shnitsel.analyze.generic import keep_norming

            self.dataset["dip_perm_norm"] = keep_norming(self.dataset["dip_perm"])
        return self.dataset.data_vars["dip_perm_norm"]

    @property
    def forces(self) -> xr.DataArray:
        if "forces" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `forces` to encode per-state forces moments in trajectory data"
            )
        return self.dataset.data_vars["forces"]

    @property
    def forces_norm(self) -> xr.DataArray:
        if "forces_norm" not in self.dataset.data_vars:
            if 'forces' not in self.dataset.data_vars:
                raise KeyError(
                    "No variable `forces` to encode per-state forces in trajectory data"
                )
            from shnitsel.analyze.generic import keep_norming

            self.dataset["forces_norm"] = keep_norming(self.dataset["forces"])
        return self.dataset.data_vars["forces_norm"]

    @property
    def forces_format(self) -> bool | Literal["all", "active_only"] | None:
        return self.dataset.attrs[]
