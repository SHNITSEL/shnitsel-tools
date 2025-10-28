import os

import xarray as xr

import shnitsel as sh


class TestSHARC:
    """Class to test functionality related to SHARC initial conditions and trajectories"""

    def test_read_dir_of_iconds(self):
        path = os.path.join("tutorials", "test_data", "sharc", "iconds_butene")
        iconds = sh.parse.sharc_icond.dirs_of_iconds(path)

        assert isinstance(iconds, xr.Dataset)

    def test_sharc_process_and_output_iconds(self):
        # parse icond data
        iconds_butene = sh.parse.sharc_icond.dirs_of_iconds(
            path="tutorials/test_data/sharc/iconds_butene"
        )
        # save the parsed data to h5netcdf file
        savepath = os.path.join(os.getcwd(), "tutorials", "test_data", "sharc", "iconds_butene.hdf5")
        iconds_butene.reset_index("statecomb").to_netcdf(savepath, engine="h5netcdf")

    def test_read_sharc_trajs(self):
        # parse trajectory data from SHARC output files
        traj_frames_butene = sh.parse.read_trajs(
            "tutorials/test_data/sharc/traj_butene", kind="sharc"
        )
