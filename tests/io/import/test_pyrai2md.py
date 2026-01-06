import os
import shnitsel as sh
from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.data.trajectory_format import Trajectory
from shnitsel.io.pyrai2md import parse_pyrai2md
from shnitsel.io import read
from shnitsel.io.pyrai2md.format_reader import PyrAI2mdFormatReader
from shnitsel.test_support.trajectory_verification import verify_trajectory_format


class TestPyrai2mdFunctionality:
    input_path = "tutorials/test_data/pyrai2md/traj_I01/"

    """Class to test functionality related to the pyrai2md file format"""

    asserted_properties_in_trajectory = [
        "energy",
        "forces",
        "atXYZ",
        "state_types",
        "state_names",
        "atNums",
        "atNames",
    ]

    def test_read_pyrai2md_trajs_direct(self):
        # parse trajectory data directly from pyrai2md output
        traj = PyrAI2mdFormatReader().read_trajectory(self.input_path + "traj1")
        assert traj is not None, (
            f"Failed to load PyRAI2md trajectory from {self.input_path + 'traj1'}"
        )

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_pyrai2md_FormatDetection_R02(self):
        # check if the directory has appropriate files
        format_info = PyrAI2mdFormatReader().check_path_for_format_info(
            self.input_path + "traj1"
        )
        assert format_info is not None, "Did not detect any format information"
        assert (
            format_info.format_name == "pyrai2md"
        ), "Did not correctly detect PyrAI2md format"

    def test_parse_wrapper_single_format(self):
        # Test if read wrapper works if kind specified
        traj = read(self.input_path + "traj1", kind="pyrai2md")
        assert traj is not None, (
            f"Failed to load PyRAI2md trajectory from {self.input_path + 'traj1'}"
        )
        assert isinstance(traj, Trajectory), (
            f"Loaded PyRAI2md trajectory from {self.input_path + 'traj1'} has wrong format: {type(traj)}"
        )

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_parse_wrapper_single_guess_kind(self):
        # Test if read wrapper works if kind specified
        traj = read(self.input_path + "traj1")
        assert traj is not None, (
            f"Failed to load PyRAI2md trajectory from {self.input_path + 'traj1'}"
        )
        assert isinstance(traj, Trajectory), (
            f"Loaded PyRAI2md trajectory from {self.input_path + 'traj1'} has wrong format: {type(traj)}"
        )

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_parse_wrapper_multi_trajectory(self):
        # Test if read wrapper works if kind specified
        traj = read(self.input_path, kind="pyrai2md")
        assert (
            traj is not None
        ), f"Failed to load PyRAI2md trajectory db from {self.input_path}"
        assert isinstance(
            traj, ShnitselDB
        ), f"Loaded PyRAI2md trajectory from {self.input_path} has wrong format: {type(traj)}"

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_parse_wrapper_multi_trajectory_guess_kind(self):
        # Test if read wrapper works if kind not specified
        traj = read(self.input_path)
        assert (
            traj is not None
        ), f"Failed to load PyRAI2md trajectory db from {self.input_path}"
        assert isinstance(
            traj, ShnitselDB
        ), f"Loaded PyRAI2md trajectory from {self.input_path} has wrong format: {type(traj)}"

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"
