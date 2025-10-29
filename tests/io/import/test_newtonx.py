import os

import xarray as xr

import shnitsel as sh
from shnitsel.io import read
from shnitsel.io.newtonx import parse_newtonx


class TestNewtonX:
    """Class to test functionality related to NewtonX file format"""

    def test_nx_direct_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = parse_newtonx(
            'tutorials/test_data/newtonx/test_trajectory_R02/TRAJ1/')

    def test_nx_general_wrapper_with_kind_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = read(
            'tutorials/test_data/newtonx/test_trajectory_R02/', kind='newtonx')

    def test_nx_detected_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = read(
            'tutorials/test_data/newtonx/test_trajectory_R02/', kind=None)
