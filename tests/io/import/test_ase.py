import os

import pytest
from xarray.testing import assert_equal
from shnitsel.io.ase import read_ase
import shnitsel.xarray


@pytest.fixture()
def FauxBulkDataDB():
    # Supposed to generate an ASE database with the wrong format so we can test if our code successfully catches wrong formats
    import pathlib

    bulk_db_path = pathlib.Path("tutorials/test_data/ase/bulk.db")
    if not bulk_db_path.exists() or not bulk_db_path.is_file():

        def build_bulk():
            from ase.build import bulk
            from ase.calculators.emt import EMT
            from ase.db import connect
            from ase.eos import calculate_eos

            db = connect(bulk_db_path.resolve().as_posix())
            for symb in ['C']:
                atoms = bulk(symb, 'diamond')
                atoms.calc = EMT()
                atoms.get_potential_energy()
                eos = calculate_eos(atoms)
                v, e, B = eos.fit()  # find minimum
                # Do one more calculation at the minimu and write to database:
                atoms.cell *= (v / atoms.get_volume()) ** (1 / 3)
                atoms.get_potential_energy()
                db.write(atoms, bm=B)

        build_bulk()
        print("Created bulk at: ", bulk_db_path.resolve().as_posix())

    return bulk_db_path.resolve().as_posix()


class TestASEFunctionality:
    """Class to test all functions of the shnitsel tools related to ASE (Atomic simulation Environment) Databases"""

    @pytest.mark.parametrize(
        'path, kind',
        [
            ('tutorials/test_data/ase/schnarc_ch2nh2+.db', 'schnet'),
            ('tutorials/test_data/ase/spainn_ch2nh2+.db', 'spainn'),
            ('tutorials/test_data/ase/schnarc_ch2nh2+.db', None),
            ('tutorials/test_data/ase/spainn_ch2nh2+.db', None),
        ],
    )
    def test_ase_round_trip(self, path, kind):
        tmp_path = '/tmp/test_round_trip.db'
        try:
            os.remove(tmp_path)
        except FileNotFoundError:
            pass

        frames1 = read_ase(path, kind=kind)
        frames1.st.write_ase(tmp_path, kind=kind)
        frames2 = read_ase(tmp_path, kind=kind)
        assert_equal(frames1, frames2)

    def test_generic_ase(self, FauxBulkDataDB):
        # Just test if we can open it. Should not fail
        frames1 = read_ase(FauxBulkDataDB, kind=None)

    def test_invalid_kinds_raises_valueerror(self, FauxBulkDataDB):
        with pytest.raises(ValueError) as excinfo:
            # Should fail because of invalid kind of DB
            frames1 = read_ase(FauxBulkDataDB, kind="harlow")
        assert "'kind' should be one of 'schnet' or 'spainn'" in str(excinfo.value)

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError) as excinfo:
            # Should fail because of file not existing
            read_ase("./nowhere.db", kind=None)

    def test_invalid_formats(self, FauxBulkDataDB):
        with pytest.raises(ValueError) as excinfo:
            # Should fail because DB is not in correct format
            read_ase(FauxBulkDataDB, kind="spainn")
        assert "No rows with the appropriate format" in str(excinfo.value)

        with pytest.raises(ValueError) as excinfo:
            # Should fail because DB is not in correct format
            read_ase(FauxBulkDataDB, kind="schnet")
        assert "No rows with the appropriate format" in str(excinfo.value)


if __name__ == '__main__':
    test = TestASEFunctionality()
    bulk_db_path = FauxBulkDataDB()
    # test.test_ase_round_trip()
    test.test_generic_ase(bulk_db_path)
    test.test_invalid_kinds_raises_valueerror(bulk_db_path)
    test.test_invalid_formats(bulk_db_path)
