import logging
import os


from shnitsel.io import read
from shnitsel.io.sharc.parse_trajectory import read_traj
from shnitsel.io.sharc.parse_initial_conditions import dir_of_iconds
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

    def test_read_dirs_of_iconds(self):
        path = os.path.join("tutorials", "test_data", "sharc", "iconds_butene")
        iconds = dir_of_iconds(path)

        assert verify_trajectory_format(
            iconds, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

    def test_sharc_process_and_output_iconds(self):
        # parse icond data
        iconds_butene = dir_of_iconds(path="tutorials/test_data/sharc/iconds_butene")
        assert verify_trajectory_format(
            iconds_butene, self.asserted_properties_in_trajectory
        ), f"Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

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
