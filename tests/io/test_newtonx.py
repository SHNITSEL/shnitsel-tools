import os

import xarray as xr

import shnitsel as sh


class TestNewtonX:
    """Class to test functionality related to NewtonX file format"""

    def test_nx():
        # parse trajectory data from Newton-X output files
        traj_frames_chd = sh.parse.read_trajs('test_data/newtonx/', kind='nx')

