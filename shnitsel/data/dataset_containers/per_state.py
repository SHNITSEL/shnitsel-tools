from dataclasses import dataclass
from typing import Literal

from shnitsel.analyze.stats import get_per_state
from .trajectory import Trajectory
from .shared import ShnitselDerivedDataset
from .frames import Frames

import xarray as xr


@dataclass
class PerState(ShnitselDerivedDataset):
    _original_frames: Frames | Trajectory

    def __init__(self, frames: Frames | Trajectory):
        assert "state" in frames.dims, (
            "Dataset is missing `state` dimension and cannot be considered an PerState set of variables."
        )
        # TODO: FIXME: Calculate per-state variables and cache in original dataset

        self._original_frames = frames
        per_state_props = get_per_state(frames.dataset)
        super().__init__(frames.dataset, per_state_props)

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
    def forces(self) -> xr.DataArray:
        if "forces" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `forces` to encode per-state forces moments in trajectory data"
            )
        return self.dataset.data_vars["forces"]

    @property
    def forces_format(self) -> bool | Literal["all", "active_only"] | None:
        return self._original_frames.forces_format
