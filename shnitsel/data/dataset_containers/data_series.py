from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Literal, Self

from ..trajectory_grouping_params import TrajectoryGroupingMetadata


from .shared import ShnitselDataset
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
class DataSeries(ShnitselDataset):
    # from .frames import Frames
    # from .per_state import PerState
    # from .inter_state import InterState

    def __init__(self, ds: xr.Dataset):
        assert 'state' in ds.dims
        assert 'atom' in ds.dims
        super().__init__(ds)

    @cached_property
    def per_state(self) -> "PerState":
        """Convert this trajectory to a PerState object only allowing access to the per-state data encoded in this entity

        Returns
        --------
            PerState: The wrapper for the per-state properties
        """
        from .per_state import PerState

        return PerState(self)

    @cached_property
    def inter_state(self) -> "InterState":
        """Convert this trajectory to an InterState object only allowing access to the inter-state data encoded in this entity.

        Will calculate some interstate properties like state-to-state energy differences.

        Returns
        --------
            InterState: The wrapper for the inter-state properties
        """
        from .inter_state import InterState

        return InterState(self)

    @property
    def leading_dim(self) -> str:
        """The leading dimension along which consistent configurations are indexed.
        Usually `time` or `frame`."""
        return "frame"

    @property
    def positions(self):
        """The atom position data stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "atXYZ" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `atXYZ` to encode positions provided for the trajectory"
            )
        return self.dataset.data_vars["atXYZ"]

    @property
    def atXYZ(self):
        """The positional data for atoms stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "atXYZ" not in self.dataset.data_vars:
            raise KeyError("No variable `atXYZ` provided for the trajectory")
        return self.dataset.data_vars["atXYZ"]

    @property
    def energy(self):
        """The energy information stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "energy" not in self.dataset.data_vars:
            raise KeyError("No variable `energy` provided for the trajectory")
        return self.dataset.data_vars["energy"]

    @property
    def forces(self):
        """The force data stored in this dataset if accessible.
        Note that depending on `forces_format`, there may only be data for the active state or
        for some of the states.

        Will throw a `KeyError` if no data is accessible."""
        if "forces" not in self.dataset.data_vars:
            raise KeyError("No variable `forces` provided for the trajectory")
        return self.dataset.data_vars["forces"]

    @property
    def nacs(self):
        """The non adiabatic coupling data stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "nacs" not in self.dataset.data_vars:
            raise KeyError("No variable `nacs` provided for the trajectory")
        return self.dataset.data_vars["nacs"]

    @property
    def socs(self):
        """The spin orbit coupling data stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "socs" not in self.dataset.data_vars:
            raise KeyError("No variable `socs` provided for the trajectory")
        return self.dataset.data_vars["socs"]

    @property
    def dipole_permanent(self):
        """The permanent dipole data stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "dip_perm" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `dip_perm` containing permanent dipole moments provided for the trajectory"
            )
        return self.dataset.data_vars["dip_perm"]

    @property
    def dipole_transition(self):
        """The transition dipole data stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "dip_trans" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `dip_trans` containing transitional dipole moments provided for the trajectory"
            )
        return self.dataset.data_vars["dip_trans"]

    @property
    def e_kin(self):
        """The kinetic energy information stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "e_kin" not in self.dataset.data_vars:
            raise KeyError("No variable `e_kin` provided for the trajectory")
        return self.dataset.data_vars["e_kin"]

    @property
    def velocities(self):
        """The velocity information stored in this dataset if accessible.

        Will throw a `KeyError` if no data is accessible."""
        if "velocities" not in self.dataset.data_vars:
            raise KeyError("No variable `velocities` provided for the trajectory")
        return self.dataset.data_vars["velocities"]

    @property
    def charge(self) -> float:
        """The charge of the molecule if set on the trajectory data.
        Loaded from `charge` attribute (or variable) or `state_charges` coordinate
        if provided.

        If no information is found, 0 is returned."""
        charge = self._param_from_vars_or_attrs('charge')
        if charge is None:
            if 'state_charges' in self._raw_dataset:
                return float(self._raw_dataset['state_charges'][0].item())
            return 0
        return float(charge)

    @charge.setter
    def charge(self, value):
        """Setter for the charge of a trajectory

        Parameters
        ----------
        float : float
            The charge in units of elementary charges.
        """
        # TODO: FIXME: We probably do not want to support setters here!
        self._raw_dataset.attrs['charge'] = value
        self._raw_dataset['state_charges'].values[:] = value

    def set_charge(self, value: float | xr.DataArray) -> Self:
        """Method to set the charge on a dataset, clear conflicting positions
        of charge info on the dataset and return a new instance of the wrapped dataset.


        Parameters
        ----------
        value : float | xr.DataArray
            Either a single value (optionally wrapped in a DataArray already) to indicate
            the charge of the full molecule in all states (will be set to coordinate `charge`) or a DataArray that represents
            state-dependent charges (which will be set to `state_charges`)

        Returns
        -------
        Self
            The updated object as a copy.

        Raises
        ------
        ValueError
            If an unsupported `value` was provided.
        """
        new_attrs = dict(self._raw_dataset.attrs)
        if 'charge' in new_attrs:
            del new_attrs['charge']

        tmp_ds = (
            self._raw_dataset.drop_attrs(deep=False)
            .assign_attrs(new_attrs)
            .drop_vars('charge', errors='ignore')
        )
        if isinstance(value, (float | int)):
            # Create a dummy, zero-dim
            charge_da = xr.DataArray(
                float(value),
                dims=[],
                attrs={'units': 'e', 'unitdim': 'charge'},
                name='charge',
            )

            new_ds = tmp_ds.assign_coords(charge=charge_da)
        elif isinstance(value, xr.DataArray):
            if len(value.sizes) == 0:
                # a scalar variable
                new_ds = tmp_ds.assign_coords(charge=value)
            else:
                # This is a state_charge contender:
                new_ds = tmp_ds.assign_coords(state_charges=value)
        else:
            raise ValueError("Unknown type for charge value: %s" % type(value))

        return type(self)(new_ds)

    def _param_from_vars_or_attrs(self, key: str) -> Any | None:
        """Helper function to extract information either from a data var or from
        a coordinate or from the attributes of the dataset

        Parameters
        ----------
        key : str
            The key under which we expect to find the data

        Returns
        -------
        Any|None
            the value associated with the key that has been found
        """
        if key in self._raw_dataset.data_vars:
            return self._raw_dataset.data_vars[key]
        if key in self._raw_dataset.coords:
            return self._raw_dataset.coords[key]
        elif key in self._raw_dataset.attrs:
            return self._raw_dataset.attrs.get(key, None)

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

    @property
    def attrs(self) -> dict:
        """A dictionary of the attributes set on this Trajectory.

        Arbitrary mapping from attribute keys (str) to attribute values.
        """
        return self.dataset.attrs

    def get_grouping_metadata(self) -> TrajectoryGroupingMetadata:
        return TrajectoryGroupingMetadata(
            delta_t_in_fs=self.delta_t,
            input_format_name=self.input_format,
            input_format_version=self.input_format_version,
            est_level=self.est_level,
            theory_basis_set=self.theory_basis_set,
            charge_in_e=self.charge,
            # TODO: FIXME: We should differentiate by all state attributes.
            num_states=len(self.state_ids),
        )
