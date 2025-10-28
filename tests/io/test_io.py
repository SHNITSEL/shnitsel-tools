import os

import xarray as xr

import shnitsel as sh


class TestGeneralIO:
    """Class to test general frame I/O functionality"""

    def test_open(self):
        sh.open_frames('tutorials/test_data/butene.nc')
