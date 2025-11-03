import logging
import os
from typing import Any, List

from shnitsel.data.TrajectoryFormat import Trajectory
import xarray as xr

from shnitsel.io import read
from shnitsel.io.sharc.parse_trajectory import read_traj
from shnitsel.io.sharc.parse_initial_conditions import dir_of_iconds
from shnitsel.units import standard_shnitsel_units
from shnitsel.units.definitions import unit_dimensions


class TestSHARC:
    """Class to test functionality related to SHARC initial conditions and trajectories"""

    def is_permitted_traj_result(self, obj: Any) -> True:
        logging.debug(f"Obj has type: {type(obj)}")
        return (
            isinstance(obj, xr.Dataset)
            or isinstance(obj, Trajectory)
            or (
                isinstance(obj, list)
                and (len(obj) == 0 or isinstance(obj[0], Trajectory))
            )
        )

    def has_required_properties(self, traj: List[Trajectory] | Trajectory) -> bool:
        res = True
        if isinstance(traj, list):
            for i_traj in traj:
                res = res and self.has_required_properties(i_traj)
        else:

            print(traj.variables.keys())
            print(traj["atXYZ"].attrs.keys())

            check_prop_units = {
                "atXYZ": {"unit": standard_shnitsel_units[unit_dimensions.length]},
                "forces": {"unit": standard_shnitsel_units[unit_dimensions.force]},
                "energy": {"unit": standard_shnitsel_units[unit_dimensions.energy]},
                "dip_trans": {"unit": standard_shnitsel_units[unit_dimensions.dipole]},
                "dip_perm": {"unit": standard_shnitsel_units[unit_dimensions.dipole]},
                "states": {},
            }

            for prop in check_prop_units:
                assert (
                    prop in traj.variables.keys()
                ), f"Property {prop} is missing in resulting trajectory"
                if "unit" in check_prop_units[prop]:
                    assert (
                        "units" in traj["atXYZ"].attrs
                    ), f"Property {prop} has no unit set in resulting trajectory"
                    required_unit = check_prop_units[prop]["unit"]
                    actual_unit = traj["atXYZ"].attrs["units"]
                    assert (
                        actual_unit == required_unit
                    ), f"Property {prop} has unit {actual_unit} instead of required unit {required_unit}"

            for var in traj:
                if "unitdim" in traj[var].attrs:
                    assert "units" in traj[var].attrs
                    unit_dim = traj[var].attrs["unitdim"]
                    assert (
                        standard_shnitsel_units[unit_dim] == "1"
                        or traj[var].attrs["units"] == standard_shnitsel_units[unit_dim]
                    )

            return True

    def test_read_dirs_of_iconds(self):
        path = os.path.join("tutorials", "test_data", "sharc", "iconds_butene")
        iconds = dir_of_iconds(path)

        assert self.is_permitted_traj_result(iconds)
        assert self.has_required_properties(iconds)

    def test_sharc_process_and_output_iconds(self):
        # parse icond data
        iconds_butene = dir_of_iconds(path="tutorials/test_data/sharc/iconds_butene")
        # save the parsed data to h5netcdf file
        savepath = os.path.join(
            os.getcwd(), "tutorials", "test_data", "sharc", "iconds_butene.hdf5"
        )
        iconds_butene.reset_index("statecomb").to_netcdf(savepath, engine="h5netcdf")

    def test_read_sharc_trajs_direct(self):
        # parse trajectory data from SHARC output files
        traj_frames_butene = read_traj(
            "tutorials/test_data/sharc/traj_butene/TRAJ_00002"
        )

        assert self.is_permitted_traj_result(traj_frames_butene)
        assert self.has_required_properties(traj_frames_butene)

    def test_read_sharc_glob_match(self):
        # Read a set of trajectories from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene/", kind="sharc"
        )
        assert self.is_permitted_traj_result(traj_frames_butene)
        assert self.has_required_properties(traj_frames_butene)

    def test_read_sharc_wrapper_direct(self):
        # Directly read one single trajectory from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene/TRAJ_00002", kind="sharc"
        )

        assert self.is_permitted_traj_result(traj_frames_butene)
        assert self.has_required_properties(traj_frames_butene)
