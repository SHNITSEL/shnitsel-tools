import os


from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.data.trajectory_format import Trajectory
from shnitsel.io import read
from shnitsel.io.helpers import LoadingParameters
from shnitsel.io.sharc.format_reader import SHARCFormatReader
from shnitsel.io.sharc.parse_trajectory import read_traj
from shnitsel.io import write_shnitsel_file
from shnitsel.test_support.trajectory_verification import verify_trajectory_format


class TestSHARCTrajectories:
    """Class to test functionality related to SHARC initial conditions and trajectories"""

    asserted_properties_in_trajectory = [
        "energy",
        "forces",
        "atXYZ",
        "state_types",
        "state_names",
        "atNums",
        "atNames",
    ]

    asserted_properties_in_gradless_trajectory = [
        "energy",
        "atXYZ",
        "state_types",
        "state_names",
        "atNums",
        "atNames",
    ]

    asserted_properties_in_v4_trajectory = [
        "energy",
        "atXYZ",
        "state_types",
        "state_names",
        "state_charges",
        "atNums",
        "atNames",
    ]

    def test_read_sharc_trajs_direct_v2_1(self):
        # parse trajectory data from SHARC output files
        traj_frames_butene = SHARCFormatReader().read_trajectory(
            "tutorials/test_data/sharc/traj_butene_v2.1/TRAJ_00002"
        )
        assert isinstance(traj_frames_butene, Trajectory)

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory + ["nacs"]
        ), "Resulting trajectories from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_glob_match_v2_1(self):
        # Read a set of trajectories from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene_v2.1/",
            kind="sharc",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory + ["nacs"]
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_detect_type_v2_1(self):
        # Read a set of trajectories from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene_v2.1/", concat_method='db'
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory + ["nacs"]
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_wrapper_direct_v2_1(self):
        # Directly read one single trajectory from a directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_butene_v2.1/",
            kind="sharc",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory + ["nacs"]
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

        write_shnitsel_file(
            traj_frames_butene, "tutorials/test_data/sharc/traj_butene.nc"
        )

    def test_read_sharc_wrapper_detect_v3_0(self):
        # Read a bundle of trajectories from a v3.0 directory using type detection
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v3.0/",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory + ["nacs"]
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_wrapper_direct_v3_0(self):
        # Read trajectory bundle from a v3.0 directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v3.0/",
            kind="sharc",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_trajectory + ["nacs"]
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

        write_shnitsel_file(
            traj_frames_butene, "tutorials/test_data/sharc/traj_I01_v3.nc"
        )

    def test_read_sharc_wrapper_detect_v3_0_triplets(self):
        # Read a bundle of trajectories from a v3.0 directory with triplet states using type detection
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v3.0_triplets/",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_gradless_trajectory
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_wrapper_direct_v3_0_triplets(self):
        # Read trajectory bundle from a v3.0 directory
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v3.0_triplets/",
            kind="sharc",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_gradless_trajectory
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

        write_shnitsel_file(
            traj_frames_butene, "tutorials/test_data/sharc/traj_I01_v3_triplets.nc"
        )

    def test_read_sharc_wrapper_detect_v3_0_triplets_nacs_socs(self):
        # Read a bundle of trajectories from a v3.0 directory using type detection containing triplet data with nacs and socs
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs/",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene,
            self.asserted_properties_in_trajectory + ["nacs", "socs"],
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

    def test_read_sharc_wrapper_direct_v3_0_triplets_nacs_socs(self):
        # Read trajectory bundle from a v3.0 directory containing triplet data with nacs and socs
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs/",
            kind="sharc",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene,
            self.asserted_properties_in_trajectory + ["nacs", "socs"],
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

        write_shnitsel_file(
            traj_frames_butene,
            "tutorials/test_data/sharc/traj_I01_v3_triplets_nacs_socs.nc",
        )

    def test_read_sharc_wrapper_detect_v4_0(self):
        # Read a bundle of trajectories from a v4.0 directory using type detection
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v4.0/",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene,
            self.asserted_properties_in_v4_trajectory,
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

        # Check charge
        assert traj_frames_butene.variables["state_charges"][0] == 1

    def test_read_sharc_wrapper_direct_v4_0(self):
        # Read trajectory bundle from a v4.0 directory. Needs to detect the charges
        traj_frames_butene = read(
            "tutorials/test_data/sharc/traj_I01_v4.0/",
            kind="sharc",
            concat_method='db',
        )
        assert isinstance(
            traj_frames_butene, ShnitselDB
        ), f"Read result was not of requested type `db`. Was {type(traj_frames_butene)} instead."

        assert verify_trajectory_format(
            traj_frames_butene, self.asserted_properties_in_v4_trajectory
        ), "Resulting trajectory from SHARC trajectory does not satisfy the Shnitsel standard format"

        # Check charge
        assert traj_frames_butene.variables["state_charges"][0] == 1

        write_shnitsel_file(
            traj_frames_butene, "tutorials/test_data/sharc/traj_I01_v4.nc"
        )
