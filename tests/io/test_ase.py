import os

import pytest
from xarray.testing import assert_equal
from shnitsel import read_ase
import shnitsel.xarray

class TestASEFunctionality:
    """Class to test all functions of the shnitsel tools related to ASE (Atomic simulation Environment) Databases"""

    @pytest.mark.parametrize(
        'path,kind',
        [
            ('tutorials/test_data/old_CH2NH2.db', 'schnet'),
            ('tutorials/test_data/tobias_cis_new.db', 'spainn'),
            ('tutorials/test_data/old_CH2NH2.db', None),
            ('tutorials/test_data/tobias_cis_new.db', None),
        ],
    )
    def test_ase_round_trip(path, kind):
        tmp_path = '/tmp/test_round_trip.db'
        try:
            os.remove(tmp_path)
        except FileNotFoundError:
            pass

        frames1 = read_ase(path, kind=kind)
        frames1.sh.write_ase(tmp_path, kind=kind)
        frames2 = read_ase(tmp_path, kind=kind)
        assert_equal(frames1, frames2)


if __name__ == '__main__':
    TestASEFunctionality.test_ase_round_trip()