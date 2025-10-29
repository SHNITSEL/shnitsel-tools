import os

import xarray as xr

from shnitsel.io import read_shnitsel_file


class TestGeneralIO:
    """Class to test general frame I/O functionality"""

    def test_open_static(self):
        read_shnitsel_file(
            'tutorials/test_data/shnitsel/fixtures/butene_static/data.nc')

    def test_open_grid(self):
        read_shnitsel_file(
            'tutorials/test_data/shnitsel/fixtures/butene_grid/data.nc')

    def test_open_dynamic(self):
        read_shnitsel_file(
            'tutorials/test_data/shnitsel/fixtures/butene_dynamic/data.nc')
