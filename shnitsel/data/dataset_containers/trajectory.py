from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Literal, Self

from ..trajectory_grouping_params import TrajectoryGroupingMetadata


from .data_series import DataSeries
from .frames import Frames
import xarray as xr


@dataclass
class MetaInformation:
    """Meta information for trajectory setup"""

    input_format: Literal["sharc", "newtonx", "ase", "pyrai2md"] | None = None
    input_type: Literal['static', 'dynamic'] | None = None
    input_format_version: str | None = None
    theory_basis_set: str | None = None
    est_level: str | None = None


@dataclass
class Trajectory(DataSeries):
    # from .frames import Frames
    # from .per_state import PerState
    # from .inter_state import InterState

    def __init__(self, ds: xr.Dataset):
        assert 'frame' not in ds.dims, (
            "Dataset has `frame` dimension and cannot be considered a Trajectory"
        )
        assert "time" in ds.dims or 'time' in ds.coords, (
            "Dataset is missing `time` dimension and cannot be considered a Trajectory"
        )
        assert "atom" in ds.dims, (
            "Dataset is missing `atom` dimension and cannot be considered a Trajectory"
        )
        assert "state" in ds.dims, (
            "Dataset is missing `state` dimension and cannot be considered a Trajectory"
        )
        super().__init__(ds)

    @cached_property
    def as_frames(self) -> "Frames":
        """Convert this trajectory to a frames version of this trajectory, where the leading dimension
        is `frame` instead of `time`.

        Returns
        --------
            Frames: The resulting frames instance with a stacked dimension `frame` and a new coordinate `active_trajectory` along the `frame` dimension
        """

        frame_ds = self.dataset.expand_dims(atrajectory=[self.trajectory_id]).stack(
            frame=["atrajectory", "time"]
        )
        return Frames(frame_ds)

    @cached_property
    def as_trajectory(self) -> Self:
        """Idempotent conversion to Trajectory instance"""
        return self

    @property
    def leading_dim(self) -> str:
        """The leading dimension along which consistent configurations are indexed.
        Usually `time` or `frame`."""
        return "time"

    @property
    def t_max(self) -> float:
        """Maximum time up to which the simulation could have run if not interrupted.

        It may actually have run to this time."""
        t_max: float | str | None = self._param_from_vars_or_attrs('t_max')
        if t_max is None:
            t_max = -1
        if isinstance(t_max, str):
            t_max = float(t_max)
        if not isinstance(t_max, float):
            t_max = -1
        return t_max

    @property
    def delta_t(self) -> float:
        """The simulation timestep usually in the same units as `time`"""
        delta_t = self._param_from_vars_or_attrs('delta_t')
        if delta_t is None:
            delta_t = -1
        if isinstance(delta_t, str):
            delta_t = float(delta_t)
        if not isinstance(delta_t, float):
            delta_t = -1
        return delta_t

    @property
    def trajectory_id(self) -> int | str | None:
        """An alias for `trajid` with a more telling name"""
        return self.trajid

    @property
    def trajid(self) -> int | str | None:
        """Id of the trajectory. If assigned it is expected to be unique across the same input
        but may clash with other trajectory ids if multiple separate imports are combined
        or indepdendent simulation data is combined."""
        trajid = self._param_from_vars_or_attrs('trajid')
        if trajid is None:
            trajid = self._param_from_vars_or_attrs('id')
        if trajid is None:
            trajid = self._param_from_vars_or_attrs('trajectory_id')
        return trajid

    @property
    def max_timestep(self) -> int:
        """Alias for `max_ts` with a more telling name"""
        return self.max_ts

    @property
    def max_ts(self) -> int:
        """The maximum time step to which the simulation progressed before termination."""
        max_ts = self._param_from_vars_or_attrs('max_ts')
        if max_ts is None:
            return self.sizes[self.leading_dimension]
        return max_ts

    @property
    def completed(self) -> bool:
        """A flag whether the imported Trajectory had successfully completed."""
        completed = self._param_from_vars_or_attrs('completed')
        if completed is None:
            return False
        return completed

    @property
    def input_format(
        self,
    ) -> Literal["sharc", "newtonx", "ase", "pyrai2md", "unknown"] | str:
        """Name of the simulation software or input file type from which the data was originally imported."""
        input_format = self._param_from_vars_or_attrs('input_format')
        if input_format is None:
            return "unknown"
        return input_format

    @property
    def input_type(self) -> Literal["static", "dynamic", "unknown"]:
        """Whether the data in this trajectory is static (independently optimized) or continuous
        time-resolved data or whether the type is not known"""
        input_type = self._param_from_vars_or_attrs('input_type')
        if input_type is None:
            return "unknown"
        return input_type

    @property
    def input_format_version(self) -> str:
        """The version of the simulation software used to create this trajectory"""
        input_format_version = self._param_from_vars_or_attrs('input_format_version')
        if input_format_version is None:
            return "unknown"
        return input_format_version

    @property
    def num_singlets(self) -> int:
        """Number of singlet states in the system"""
        num_singlets = self._param_from_vars_or_attrs('num_singlets')
        if num_singlets is None:
            num_singlets = self._param_from_vars_or_attrs('nsinglets')
        if num_singlets is None:
            return 0
        return num_singlets

    @property
    def num_doublets(self) -> int:
        """Number of doublet states in the system"""
        num_doublets = self._param_from_vars_or_attrs('num_doublets')
        if num_doublets is None:
            num_doublets = self._param_from_vars_or_attrs('ndoublets')
        if num_doublets is None:
            return 0
        return num_doublets

    @property
    def num_triplets(self) -> int:
        """Number of triplet states in the system"""
        num_triplets = self._param_from_vars_or_attrs('num_triplets')
        if num_triplets is None:
            num_triplets = self._param_from_vars_or_attrs('ntriplets')
        if num_triplets is None:
            return 0
        return num_triplets

    @property
    def forces_format(self) -> bool | Literal["all", "active_only"] | None:
        """The `forces` format in the trajectory.

        Options are a binary flag to signify whether there are forces or not.
        If the flag is True, the forces still might not be available for all states but only for the active state.
        If `'all'` is the format, then there will be forces for all states.
        If the mode is `'active_only'` there will definitely only be forces for the active state in the trajectory.
        If The mode is `None`, more specific manual analysis may be required."""
        has_forces = self._param_from_vars_or_attrs('has_forces')
        if has_forces is None:
            has_forces = self._param_from_vars_or_attrs('forces_format')
        return has_forces

    @property
    def is_multi_trajectory(self) -> bool:
        """Flag whether this is a multi-trajectory container.

        Overwritten by child classes that combine multiple trajectories into one object
        """
        return False

    @property
    def trajectory_input_path(self) -> str | None:
        """Input path from which the trajectory was loaded"""
        trajectory_input_path = self._param_from_vars_or_attrs('trajectory_input_path')
        return trajectory_input_path

    @property
    def theory_basis_set(self) -> str | None:
        """The theory basis set identifier for the underlying simulation"""
        theory_basis_set = self._param_from_vars_or_attrs('theory_basis_set')
        return theory_basis_set

    @property
    def est_level(self) -> str | None:
        """The electronic structure theory level used during the simulation."""
        est_level = self._param_from_vars_or_attrs('est_level')
        return est_level

    @property
    def misc_input_settings(self) -> dict | None:
        """A dictionary of miscalleneous input settings read from trajectory output

        Arbitrary mapping from file names to settings within those files.
        """
        # To keep track of input settings we do not explicitly use anywhere else.
        misc_input_settings = self._param_from_vars_or_attrs('misc_input_settings')
        return misc_input_settings
