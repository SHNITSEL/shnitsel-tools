import logging
import os


from shnitsel.io import read
from shnitsel.io.helpers import LoadingParameters
from shnitsel.io.sharc.format_reader import SHARCFormatReader
from shnitsel.io.sharc.parse_trajectory import read_traj
from shnitsel.io.sharc.parse_initial_conditions import read_iconds_individual
from shnitsel.io import write_shnitsel_file
from shnitsel.test_support.TrajectoryVerification import verify_trajectory_format


class TestSHARC:
    """Class to test functionality related to SHARC initial conditions and trajectories"""

    asserted_properties_in_trajectory = [
        "energy",
        "forces",
        "atXYZ",
        "state_types",
        "state_names",
        "nacs",
        "atNums",
        "atNames",
    ]

    def test_read_single_iconds_folder_direct(self):
        # Parse iconds data with direct method call
        path = os.path.join(
            "tutorials", "test_data", "sharc", "iconds_butene", "ICOND_00000"
        )
        iconds = SHARCFormatReader().read_trajectory(path, loading_parameters=LoadingParameters())

        assert verify_trajectory_format(
            iconds, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

    def test_read_multi_iconds_folder_wrapper_detect(self):
        # parse icond data with wrapper
        iconds_butene = read("tutorials/test_data/sharc/iconds_butene")
        assert verify_trajectory_format(
            iconds_butene, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

    def test_read_multi_iconds_folder_wrapper_specific(self):
        # Parse iconds data with wrapper but specified format
        iconds_butene = read("tutorials/test_data/sharc/iconds_butene", kind="sharc")
        assert verify_trajectory_format(
            iconds_butene, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

    def test_sharc_process_and_output_iconds_directory(self):
        # parse icond data
        iconds_butene = read("tutorials/test_data/sharc/iconds_butene")
        assert verify_trajectory_format(
            iconds_butene, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

        # save the parsed data to h5netcdf file
        savepath = os.path.join(
            os.getcwd(), "tutorials", "test_data", "sharc", "iconds_butene.hdf5"
        )
        assert iconds_butene is not None
        if not isinstance(iconds_butene, list):
            write_shnitsel_file(iconds_butene, savepath)

    def test_read_sharc_trajs_direct(self):
        # parse trajectory data from SHARC output files
        traj_frames_butene = SHARCFormatReader().read_trajectory(
            "tutorials/test_data/sharc/traj_butene/TRAJ_00002"
        )
        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_glob_match(self):
        # Read a set of trajectories from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene/", kind="sharc"
        )
        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_wrapper_direct(self):
        # Directly read one single trajectory from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene/TRAJ_00002", kind="sharc"
        )

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"
