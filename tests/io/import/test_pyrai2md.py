import os
import shnitsel as sh
from shnitsel.io.pyrai2md import parse_pyrai2md
from shnitsel.io import read
from shnitsel.io.pyrai2md.format_reader import PyrAI2mdFormatReader
from shnitsel.test_support.TrajectoryVerification import verify_trajectory_format


class TestPyrai2mdFunctionality:
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
        traj = PyrAI2mdFormatReader().read_trajectory("tutorials/test_data/pyrai2md/traj")

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_pyrai2md_FormatDetection_R02(self):
        # check if the directory has appropriate files
        format_info = PyrAI2mdFormatReader().check_path_for_format_info(
            "tutorials/test_data/pyrai2md/traj"
        )
        assert format_info is not None, "Did not detect any format information"
        assert (
            format_info.format_name == "pyrai2md"
        ), "Did not correctly detect PyrAI2md format"

    def test_parse_wrapper_single_format(self):
        # Test if read wrapper works if kind specified
        traj = read("tutorials/test_data/pyrai2md/traj", kind="pyrai2md")

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_parse_wrapper_single_guess_kind(self):
        # Test if read wrapper works if kind specified
        traj = read("tutorials/test_data/pyrai2md/traj")

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_parse_wrapper_trajectory(self):
        # Test if read wrapper works if kind specified
        traj = read("tutorials/test_data/pyrai2md/", kind="pyrai2md")

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"

    def test_parse_wrapper_trajectory_guess_kind(self):
        # Test if read wrapper works if kind specified
        traj = read("tutorials/test_data/pyrai2md/traj")

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from PyrAI2md trajectory does not satisfy the Shnitsel standard format"
