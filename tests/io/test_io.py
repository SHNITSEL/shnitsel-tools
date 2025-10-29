import os

import xarray as xr

from shnitsel.io import read_shnitsel_file


class TestGeneralIO:
    """Class to test general frame I/O functionality"""

    def test_open(self):
        read_shnitsel_file('tutorials/test_data/fixtures/butene/data.nc')
