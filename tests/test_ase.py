import os

from xarray.testing import assert_equal
from shnitsel.dynamic.ase import read_ase_db, write_ase_db


def test_ase_round_trip():
    tmp_path = '/tmp/test_round_trip.db'
    try:
        os.remove(tmp_path)
    except FileNotFoundError:
        pass

    frames1 = read_ase_db('/nc/SHNITSEL-data/old_CH2NH2.db')
    write_ase_db(frames1, tmp_path)
    frames2 = read_ase_db(tmp_path)
    assert_equal(frames1, frames2)


if __name__ == '__main__':
    test_ase_round_trip()