from dataclasses import dataclass

from shnitsel.analyze.generic import keep_norming
from shnitsel.analyze.stats import get_inter_state
from .shared import ShnitselDerivedDataset
from .frames import Frames
from .trajectory import Trajectory
import xarray as xr


@dataclass
class InterState(ShnitselDerivedDataset):
    _original_frames: Frames | Trajectory | None

    def __init__(
        self,
        frames: Frames | Trajectory | None = None,
        /,
        direct_interstate_data: xr.Dataset | None = None,
    ):
        base_ds = None
        if frames is not None:
            assert (
                "state" in frames.dataset.dims
            ), "Dataset is missing `state` dimension and cannot be considered an InterState set of variables."
            # TODO: FIXME: Calculate per-state variables and cache in original dataset
            base_ds = frames.dataset

        if direct_interstate_data is not None:
            inter_state_props = direct_interstate_data
        elif base_ds is not None:
            inter_state_props = get_inter_state(base_ds)
        else:
            inter_state_props = xr.Dataset()

        self._original_frames = frames

        super().__init__(base_ds, inter_state_props)

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
    def dipole_transition_norm(self) -> xr.DataArray:
        if "dipole_trans_norm" not in self.dataset.data_vars:
            if not self.has_variable('dip_trans'):
                raise KeyError(
                    "No variable `dipole_trans_norm` to encode the norm of interstate transition dipole moments in trajectory data. Also no raw `dipole_trans` to calculate it from. "
                )
            self.dataset['dip_trans_norm'] = keep_norming(self.dipole_transition)
        return self.dataset.data_vars["dip_trans_norm"]

    @property
    def nacs(self) -> xr.DataArray:
        if "nacs" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `nacs` to encode non-adiabatic couplings in trajectory data"
            )
        return self.dataset.data_vars["nacs"]

    @property
    def nacs_norm(self) -> xr.DataArray:
        if "nacs_norm" not in self.dataset.data_vars:
            if not self.has_variable('nacs'):
                raise KeyError(
                    "No variable `nacs_norm` to encode the norm of non-adiabatic-couplings in trajectory data. Also no raw `nacs` to calculate it from. "
                )
            self.dataset['nacs_norm'] = keep_norming(self.nacs)
        return self.dataset.data_vars["nacs_norm"]

    @property
    def socs(self) -> xr.DataArray:
        if "socs" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `socs` to encode spin-orbit-couplings in trajectory data"
            )
        return self.dataset.data_vars["socs"]

    @property
    def socs_norm(self) -> xr.DataArray:
        if "socs_norm" not in self.dataset.data_vars:
            if not self.has_variable('socs'):
                raise KeyError(
                    "No variable `socs_norm` to encode the norm of spin-orbit-couplings in trajectory data. Also no raw `socs` to calculate it from. "
                )
            self.dataset['socs_norm'] = keep_norming(self.socs)
        return self.dataset.data_vars["socs_norm"]

    @property
    def fosc(self) -> xr.DataArray:
        if "fosc" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `fosc` to encode the strength of the oscillator in trajectory data"
            )
        return self.dataset.data_vars["fosc"]
