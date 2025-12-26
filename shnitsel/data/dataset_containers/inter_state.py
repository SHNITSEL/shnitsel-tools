from dataclasses import dataclass

from shnitsel.analyze.stats import get_inter_state
from .shared import ShnitselDerivedDataset
from .frames import Frames
import xarray as xr


@dataclass
class InterState(ShnitselDerivedDataset):
    _original_frames: Frames

    def __init__(self, frames: Frames):
        assert (
            "state" in frames.dataset.dims
        ), "Dataset is missing `state` dimension and cannot be considered an PerState set of variables."
        # TODO: FIXME: Calculate per-state variables and cache in original dataset

        self._original_frames = frames
        per_state_props = get_inter_state(frames.dataset)
        super().__init__(frames.dataset, per_state_props)

    @property
    def delta_energy(self):
        pass

    @property
    def dipole_transition(self):
        pass

    @property
    def nacs(self):
        pass

    @property
    def socs(self):
        pass
