from dataclasses import dataclass

from shnitsel.analyze.stats import get_per_state
from .shared import ShnitselDerivedDataset
from .frames import Frames
import xarray as xr


@dataclass
class PerState(ShnitselDerivedDataset):
    _original_frames: Frames

    def __init__(self, frames: Frames):
        assert (
            "state" in frames.dataset.dims
        ), "Dataset is missing `state` dimension and cannot be considered an PerState set of variables."
        # TODO: FIXME: Calculate per-state variables and cache in original dataset

        self._original_frames = frames
        per_state_props = get_per_state(frames.dataset)
        super().__init__(frames.dataset, per_state_props)
