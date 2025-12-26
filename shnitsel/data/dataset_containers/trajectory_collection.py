from dataclasses import dataclass
from typing import Sequence

from .trajectory import Trajectory


@dataclass
class TrajectoryCollection:
    _trajectories: Sequence[Trajectory]

    def __init__(self, trajectories: Sequence[Trajectory]):
        self._trajectories = trajectories
