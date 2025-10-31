import os

import xarray as xr

import shnitsel as sh
from shnitsel.io import read
from shnitsel.io.newtonx import parse_newtonx
from shnitsel.io.newtonx.format_reader import NewtonXFormatReader


class TestNewtonX:
    """Class to test functionality related to NewtonX file format"""

    def test_nx_direct_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = parse_newtonx(
            "tutorials/test_data/newtonx/test_trajectory_R02/TRAJ1/"
        )

    def test_nx_FormatDetection_R02(self):
        # check if the directory has appropriate files
        format_info = NewtonXFormatReader().check_path_for_format_info(
            "tutorials/test_data/newtonx/test_trajectory_R02/TRAJ1/"
        )

    def test_nx_FormatReader_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = NewtonXFormatReader().read_from_path(
            "tutorials/test_data/newtonx/test_trajectory_R02/TRAJ1/"
        )

    def test_nx_general_wrapper_with_kind_R02_multiple(self):
        # parse multiple trajectory datasets from Newton-X output directories
        traj_frames_chd = read(
            "tutorials/test_data/newtonx/test_trajectory_R02/", kind="newtonx"
        )

    def test_nx_general_wrapper_with_kind_R02_single(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = read(
            "tutorials/test_data/newtonx/test_trajectory_R02/TRAJ1", kind="newtonx"
        )

    def test_nx_detected_R02(self):
        # parse trajectory data from Newton-X output files
        traj_frames_chd = read(
            "tutorials/test_data/newtonx/test_trajectory_R02/", kind=None
        )
