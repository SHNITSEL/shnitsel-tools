import os

import xarray as xr

from shnitsel.io.shnitsel.format_reader import ShnitselFormatReader


class TestShnitselIO:
    """Class to test general frame I/O functionality"""

    def test_open_static(self):
        ShnitselFormatReader().read_data(
            'tutorials/test_data/shnitsel/fixtures/butene_static/data.nc'
        )

    def test_open_grid(self):
        ShnitselFormatReader().read_data(
            'tutorials/test_data/shnitsel/fixtures/butene_grid/data.nc'
        )

    def test_open_dynamic(self):
        ShnitselFormatReader().read_data(
            'tutorials/test_data/shnitsel/fixtures/butene_dynamic/data.nc'
        )
