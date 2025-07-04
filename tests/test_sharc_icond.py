import os

import xarray as xr

import shnitsel as sh


def test_dir_of_iconds():
    path = os.path.join('tutorials', 'test_data', 'sharc', 'iconds_butene')
    iconds = sh.parse.sharc_icond.dirs_of_iconds(path)

    assert isinstance(iconds, xr.Dataset)