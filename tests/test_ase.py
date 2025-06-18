import os

import pytest
from xarray.testing import assert_equal
from shnitsel.dynamic.ase import read_ase_db, write_ase_db

@pytest.mark.parametrize(
    'path,kind',
    [
        ('/nc/SHNITSEL-data/old_CH2NH2.db', 'schnet'),
        ('/nc/SHNITSEL-data/tobias_cis_new.db', 'spainn'),
    ],
)
def test_ase_round_trip(path, kind):
    tmp_path = '/tmp/test_round_trip.db'
    try:
        os.remove(tmp_path)
    except FileNotFoundError:
        pass

    frames1 = read_ase_db(path, kind=kind)
    write_ase_db(frames1, tmp_path, kind=kind)
    frames2 = read_ase_db(tmp_path, kind=kind)
    assert_equal(frames1, frames2)


if __name__ == '__main__':
    test_ase_round_trip()