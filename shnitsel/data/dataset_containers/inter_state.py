from dataclasses import dataclass

from shnitsel.analyze.stats import get_inter_state
from .shared import ShnitselDerivedDataset
from .frames import Frames
import xarray as xr


@dataclass
class InterState(ShnitselDerivedDataset):
    _original_frames: Frames

    def __init__(self, frames: Frames):
        assert "state" in frames.dataset.dims, (
            "Dataset is missing `state` dimension and cannot be considered an PerState set of variables."
        )
        # TODO: FIXME: Calculate per-state variables and cache in original dataset

        self._original_frames = frames
        per_state_props = get_inter_state(frames.dataset)
        super().__init__(frames.dataset, per_state_props)

    @property
    def delta_energy(self) -> xr.DataArray:
        return self.energy_interstate

    @property
    def energy_interstate(self) -> xr.DataArray:
        if "energy_interstate" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `energy` to calculate interstate delta energy in trajectory data"
            )
        return self.dataset.data_vars["energy_interstate"]

    @property
    def dipole_transition(self) -> xr.DataArray:
        if "dip_trans" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `dip_trans` to encode interstate transition dipole moments in trajectory data"
            )
        return self.dataset.data_vars["dip_trans"]

    @property
    def nacs(self) -> xr.DataArray:
        if "nacs" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `nacs` to encode non-adiabatic couplings in trajectory data"
            )
        return self.dataset.data_vars["nacs"]

    @property
    def socs(self) -> xr.DataArray:
        if "socs" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `socs` to encode spin-orbit-couplings in trajectory data"
            )
        return self.dataset.data_vars["socs"]
