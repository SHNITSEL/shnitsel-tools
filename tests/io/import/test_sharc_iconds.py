import os


from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.data.trajectory_format import Trajectory
from shnitsel.io import read
from shnitsel.io.helpers import LoadingParameters
from shnitsel.io.sharc.format_reader import SHARCFormatReader
from shnitsel.io import write_shnitsel_file
from shnitsel.test_support.trajectory_verification import verify_trajectory_format


class TestSHARC:
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

    def test_read_single_iconds_folder_direct(self):
        # Parse iconds data with direct method call
        path = "tutorials/test_data/sharc/iconds_butene/ICOND_00000"
        iconds = SHARCFormatReader().read_trajectory(
            path, loading_parameters=LoadingParameters()
        )
        assert isinstance(
            iconds, Trajectory
        ), "Wrong return format. Was expecting `Trajectory` result for single trajectory folder"

        assert verify_trajectory_format(
            iconds, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

    def test_read_multi_iconds_folder_wrapper_detect(self):
        # parse icond data with wrapper
        iconds_butene = read(
            "tutorials/test_data/sharc/iconds_butene", concat_method='db'
        )

        assert isinstance(iconds_butene, ShnitselDB)

        assert verify_trajectory_format(
            iconds_butene, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

    def test_read_multi_iconds_folder_wrapper_specific(self):
        # Parse iconds data with wrapper but specified format
        iconds_butene = read("tutorials/test_data/sharc/iconds_butene", kind="sharc")
        assert isinstance(iconds_butene, ShnitselDB)

        assert verify_trajectory_format(
            iconds_butene, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

    def test_sharc_process_and_output_iconds_directory(self):
        # parse icond data
        iconds_butene = read("tutorials/test_data/sharc/iconds_butene")
        assert isinstance(iconds_butene, ShnitselDB)

        assert verify_trajectory_format(
            iconds_butene, self.asserted_properties_in_trajectory
        ), "Resulting trajectory from initial conditions does not satisfy the Shnitsel standard format"

        # save the parsed data to h5netcdf file
        savepath = os.path.join(
            os.getcwd(), "tutorials", "test_data", "sharc", "iconds_butene.hdf5"
        )
        assert iconds_butene is not None
        if not isinstance(iconds_butene, list):
            write_shnitsel_file(iconds_butene, savepath)
