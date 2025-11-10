import os

from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.io.shnitsel.write import write_shnitsel_file
from shnitsel.test_support.trajectory_verification import verify_trajectory_format
import xarray as xr

import shnitsel as sh
from shnitsel.io import read
from shnitsel.io.newtonx import parse_newtonx
from shnitsel.io.newtonx.format_reader import NewtonXFormatReader


class TestNewtonX:
    """Class to test functionality related to NewtonX file format"""

    asserted_properties_in_trajectory = [
        "energy",
        "forces",
        "atXYZ",
        "state_types",
        "nacs",
        "astate",
        "atNums",
        "atNames",
    ]

    def test_nx_direct_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = NewtonXFormatReader().read_trajectory(
            "tutorials/test_data/newtonx/test_R02/TRAJ1/"
        )
        assert verify_trajectory_format(
            traj_frames_chd, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from NewtonX trajectory does not satisfy the Shnitsel standard format"

    def test_nx_FormatDetection_R02(self):
        # check if the directory has appropriate files
        format_info = NewtonXFormatReader().check_path_for_format_info(
            "tutorials/test_data/newtonx/test_R02/TRAJ1/"
        )
        assert format_info is not None, "Did not detect any format information"
        assert (
            format_info.format_name == "newtonx"
        ), "Did not correctly detect NewtonX format"

    def test_nx_FormatReader_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = NewtonXFormatReader().read_trajectory(
            "tutorials/test_data/newtonx/test_R02/TRAJ1/"
        )
        assert verify_trajectory_format(
            traj_frames_chd, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from NewtonX trajectory does not satisfy the Shnitsel standard format"

    def test_nx_general_wrapper_with_kind_R02_multiple(self):
        # parse multiple trajectory datasets from Newton-X output directories
        traj_frames_chd = read("tutorials/test_data/newtonx/test_R02/", kind="newtonx")
        assert verify_trajectory_format(
            traj_frames_chd, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from NewtonX trajectory does not satisfy the Shnitsel standard format"

    def test_nx_general_wrapper_with_kind_R02_single(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = read(
            "tutorials/test_data/newtonx/test_R02/TRAJ1", kind="newtonx"
        )
        assert verify_trajectory_format(
            traj_frames_chd, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from NewtonX trajectory does not satisfy the Shnitsel standard format"

    def test_nx_detected_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = read("tutorials/test_data/newtonx/test_R02/", kind=None)
        assert verify_trajectory_format(
            traj_frames_chd, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from NewtonX trajectory does not satisfy the Shnitsel standard format"

    def test_nx_detected_R02_db_tree(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = read(
            "tutorials/test_data/newtonx/test_R02/", kind=None, concat_method="db"
        )
        assert isinstance(
            traj_frames_chd, ShnitselDB
        ), "Failed to build ShnitselDB structure from trajectories"

        print(traj_frames_chd)
        write_shnitsel_file(traj_frames_chd, "tutorials/test_data/newtonx/test_R02.nc/")
        loaded_db = read("tutorials/test_data/newtonx/test_R02.nc/")
        assert isinstance(
            loaded_db, ShnitselDB
        ), "Failed to read ShnitselDB structure from stored file"
