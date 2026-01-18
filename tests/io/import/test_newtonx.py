import pathlib

import pytest

from shnitsel.data.dataset_containers import Frames, Trajectory
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.io.shnitsel.write import write_shnitsel_file
from shnitsel.test_support.trajectory_verification import verify_trajectory_format

from shnitsel.io import read
from shnitsel.io.newtonx.format_reader import NewtonXFormatReader


input_pyrazene = "tutorials/test_data/newtonx/test_pyrazene_v2.6/"
input_R02_a = "tutorials/test_data/newtonx/test_R02_a_v2.2/"
input_R02_b = "tutorials/test_data/newtonx/test_R02_b/"

I01_additional_properties = ["forces", "astate", "e_kin", "velocities"]
R02_a_additional_properties = ["forces", "astate", "nacs"]
R02_b_additional_properties = []


class TestNewtonX:
    """Class to test functionality related to NewtonX file format"""

    asserted_properties_in_trajectory = [
        "energy",
        "atXYZ",
        "state_types",
        "atNums",
        "atNames",
    ]

    @pytest.mark.parametrize(
        'path, add_props',
        [
            (input_pyrazene + "TRAJ1/", I01_additional_properties),
            (input_R02_a + "TRAJ1/", R02_a_additional_properties),
            (input_R02_b + "TRAJ49/", R02_b_additional_properties),
        ],
    )
    def test_nx_direct_R02_single_direct(self, path, add_props):
        # parse trajectory data from Newton-X output files
        traj = NewtonXFormatReader().read_data(path)
        assert traj is not None, f"Failed to load NewtonX trajectory from {path}"
        assert isinstance(traj, Trajectory), (
            f"Loaded NewtonX trajectory from {path} has wrong format: {type(traj)}"
        )

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory + add_props
        ), (
            f"Resulting trajectory from NewtonX trajectory at {path} does not satisfy the Shnitsel standard format"
        )

    @pytest.mark.parametrize(
        'path, add_props',
        [
            (input_pyrazene + "TRAJ1/", I01_additional_properties),
            (input_R02_a + "TRAJ1/", R02_a_additional_properties),
            (input_R02_b + "TRAJ49/", R02_b_additional_properties),
        ],
    )
    def test_nx_FormatDetection_R02(self, path, add_props):
        # check if the directory has appropriate files
        format_info = NewtonXFormatReader().check_path_for_format_info(path)
        assert format_info is not None, (
            f"Did not detect any format information for {path}"
        )
        assert format_info.format_name == "newtonx", (
            f"Did not correctly detect NewtonX format for {path}"
        )

    @pytest.mark.parametrize(
        'path, add_props',
        [
            (input_pyrazene, I01_additional_properties),
            (input_R02_a, R02_a_additional_properties),
            (input_R02_b, R02_b_additional_properties),
        ],
    )
    def test_nx_general_wrapper_with_kind_R02_multiple(self, path, add_props):
        # parse multiple trajectory datasets from Newton-X output directories
        traj = read(path, kind="newtonx")
        assert traj is not None, f"Failed to load NewtonX trajectory db from {path}"
        assert isinstance(traj, ShnitselDB), (
            f"Loaded NewtonX db from {path} has wrong format: {type(traj)}"
        )

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory + add_props
        ), (
            f"Resulting trajectory db from NewtonX trajectory from {path} does not satisfy the Shnitsel standard format"
        )

    @pytest.mark.parametrize(
        'path, add_props',
        [
            (input_pyrazene + "TRAJ1/", I01_additional_properties),
            (input_R02_a + "TRAJ1/", R02_a_additional_properties),
            (input_R02_b + "TRAJ49/", R02_b_additional_properties),
        ],
    )
    def test_nx_general_wrapper_with_kind_R02_single(self, path, add_props):
        # parse trajectory data from Newton-X output files
        traj = read(path, kind="newtonx")
        assert traj is not None, f"Failed to load NewtonX trajectory from {path}"
        assert isinstance(traj, Trajectory), (
            f"Loaded NewtonX trajectory from {path} has wrong format: {type(traj)}"
        )

        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory + add_props
        ), (
            f"Resulting trajectory from NewtonX trajectory from {path} does not satisfy the Shnitsel standard format"
        )

    @pytest.mark.parametrize(
        'path, add_props',
        [
            (input_pyrazene + "TRAJ1/", I01_additional_properties),
            (input_R02_a + "TRAJ1/", R02_a_additional_properties),
            (input_R02_b + "TRAJ49/", R02_b_additional_properties),
        ],
    )
    def test_nx_detected_R02_single(self, path, add_props):
        # parse trajectory data from Newton-X output files
        traj = read(path, kind=None)
        assert traj is not None, f"Failed to load NewtonX trajectory from {path}"
        assert isinstance(traj, Trajectory), (
            f"Loaded NewtonX trajectory from {path} has wrong format: {type(traj)}"
        )
        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory + add_props
        ), (
            f"Resulting trajectory from NewtonX trajectory from {path} does not satisfy the Shnitsel standard format"
        )

    @pytest.mark.parametrize(
        'path, add_props',
        [
            (input_pyrazene, I01_additional_properties),
            (input_R02_a, R02_a_additional_properties),
            (input_R02_b, R02_b_additional_properties),
        ],
    )
    def test_nx_detected_R02_multi(self, path, add_props):
        # parse trajectory data from Newton-X output files
        traj = read(path, kind=None)
        assert traj is not None, f"Failed to load NewtonX trajectory from {path}"
        assert isinstance(traj, ShnitselDB), (
            f"Loaded NewtonX trajectory from {path} has wrong format: {type(traj)}"
        )
        assert verify_trajectory_format(
            traj, self.asserted_properties_in_trajectory + add_props
        ), (
            f"Resulting trajectory from NewtonX trajectory db from {path} does not satisfy the Shnitsel standard format"
        )

    @pytest.mark.parametrize(
        'path, add_props',
        [
            (input_pyrazene, I01_additional_properties),
            (input_R02_a, R02_a_additional_properties),
            (input_R02_b, R02_b_additional_properties),
        ],
    )
    def test_nx_detected_R02_multi_round_trip(self, path, add_props):
        traj = read(path, kind=None)
        assert traj is not None, f"Failed to load NewtonX trajectory from {path}"
        assert isinstance(traj, ShnitselDB), (
            f"Loaded NewtonX trajectory from {path} has wrong format: {type(traj)}"
        )
        path_obj = pathlib.Path(path)
        db_path = path_obj.parent / (path_obj.name + ".nc")
        write_shnitsel_file(traj, db_path)
        loaded_db = read(db_path)
        assert isinstance(loaded_db, ShnitselDB), (
            "Failed to read ShnitselDB structure from stored file"
        )
