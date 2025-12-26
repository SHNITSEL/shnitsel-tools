from dataclasses import dataclass
from functools import cached_property
from typing import Literal


from .shared import ShnitselDataset
import xarray as xr


@dataclass
class Trajectory(ShnitselDataset):
    from .frames import Frames

    def __init__(self, ds: xr.Dataset):
        assert (
            "time" in ds.dims
        ), "Dataset is missing `time` dimension and cannot be considered a Trajectory"
        assert (
            "atom" in ds.dims
        ), "Dataset is missing `atom` dimension and cannot be considered a Trajectory"
        assert (
            "state" in ds.dims
        ), "Dataset is missing `state` dimension and cannot be considered a Trajectory"
        super().__init__(ds)

    @cached_property
    def as_frames(self) -> "Frames":
        """Convert this trajectory to a frames version of this trajectory, where the leading dimension
        is `frame` instead of `time`.

        Returns
        --------
            Frames: The resulting frames instance with a stacked dimension `frame` and a new coordinate `active_trajectory` along the `frame` dimension
        """
        from .frames import Frames

        frame_ds = self.dataset.expand_dims(
            active_trajectory=[self.dataset.attrs["trajid"]]
        ).stack(frame=["active_trajectory", "time"])
        return Frames(frame_ds)

    @property
    def leading_dim(self) -> str:
        return "time"

    @property
    def positions(self):
        if "atXYZ" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `atXYZ` to encode positions provided for the trajectory"
            )
        return self.dataset.data_vars["atXYZ"]

    @property
    def atXYZ(self):
        if "atXYZ" not in self.dataset.data_vars:
            raise KeyError("No variable `atXYZ` provided for the trajectory")
        return self.dataset.data_vars["atXYZ"]

    @property
    def energy(self):
        if "energy" not in self.dataset.data_vars:
            raise KeyError("No variable `energy` provided for the trajectory")
        return self.dataset.data_vars["energy"]

    @property
    def forces(self):
        if "forces" not in self.dataset.data_vars:
            raise KeyError("No variable `forces` provided for the trajectory")
        return self.dataset.data_vars["forces"]

    @property
    def nacs(self):
        if "nacs" not in self.dataset.data_vars:
            raise KeyError("No variable `nacs` provided for the trajectory")
        return self.dataset.data_vars["nacs"]

    @property
    def socs(self):
        if "socs" not in self.dataset.data_vars:
            raise KeyError("No variable `socs` provided for the trajectory")
        return self.dataset.data_vars["socs"]

    @property
    def dipole_permanent(self):
        if "dip_perm" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `dip_perm` containing permanent dipole moments provided for the trajectory"
            )
        return self.dataset.data_vars["dip_perm"]

    @property
    def dipole_transition(self):
        if "dip_trans" not in self.dataset.data_vars:
            raise KeyError(
                "No variable `dip_trans` containing transitional dipole moments provided for the trajectory"
            )
        return self.dataset.data_vars["dip_trans"]

    @property
    def e_kin(self):
        if "e_kin" not in self.dataset.data_vars:
            raise KeyError("No variable `e_kin` provided for the trajectory")
        return self.dataset.data_vars["e_kin"]

    @property
    def velocities(self):
        if "velocities" not in self.dataset.data_vars:
            raise KeyError("No variable `velocities` provided for the trajectory")
        return self.dataset.data_vars["velocities"]

    @property
    def charge(self) -> float:
        pass

    @property
    def t_max(self) -> float:
        pass

    @property
    def delta_t(self) -> float:
        pass

    @property
    def trajectory_id(self) -> int:
        pass

    @property
    def max_timestep(self) -> int:
        pass

    @property
    def completed(self) -> bool:
        pass

    @property
    def input_format(self) -> Literal["sharc", "newtonx", "ase", "pyrai2md"]:
        pass

    @property
    def input_type(self) -> Literal["static", "dynamic", "unknown"]:
        pass

    @property
    def input_format_version(self) -> str:
        pass

    @property
    def num_singlets(self) -> int:
        pass

    @property
    def num_doublets(self) -> int:
        pass

    @property
    def num_triplets(self) -> int:
        pass

    @property
    def forces_format(self) -> bool | Literal["all", "active_only"] | None:
        pass

    @property
    def is_multi_trajectory(self) -> bool:
        return False

    @property
    def trajectory_input_path(self) -> str | None:
        pass

    @property
    def theory_basis_set(self) -> str | None:
        pass

    @property
    def est_level(self) -> str | None:
        pass

    @property
    def misc_input_settings(self) -> dict | None:
        # To keep track of input settings we do not explicitly use anywhere else.
        pass

    @property
    def attrs(self) -> dict:
        return self.dataset.attrs
