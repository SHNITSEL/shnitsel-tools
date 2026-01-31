from dataclasses import dataclass
from typing import Iterable, Sequence

from .trajectory import Trajectory


@dataclass
class TrajectoryCollection:
    _trajectories: Sequence[Trajectory]

    def __init__(self, trajectories: Sequence[Trajectory]):
        self._trajectories = trajectories

    @property
    def trajectories(self) -> Iterable[Trajectory]:
        return self._trajectories
