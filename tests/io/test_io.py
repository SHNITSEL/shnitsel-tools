import os

import xarray as xr

import shnitsel as sh


class TestGeneralIO:
    """Class to test general frame I/O functionality"""

    def test_open():
        sh.open_frames('test_data/butene.nc')
