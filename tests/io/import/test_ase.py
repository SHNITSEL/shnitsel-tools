import os
import pathlib

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
        'path, db_format',
        [
            ('tutorials/test_data/ase/schnarc_ch2nh2+.db', 'schnet'),
            ('tutorials/test_data/ase/spainn_ch2nh2+.db', 'spainn'),
            # ('tutorials/test_data/ase/schnarc_ch2nh2+.db', None),
            # ('tutorials/test_data/ase/spainn_ch2nh2+.db', None),
        ],
    )
    def test_ase_round_trip(self, path, db_format):
        tmp_path = 'tutorials/test_data/ase/tmp_test_round_trip.db'
        try:
            os.remove(tmp_path)
        except FileNotFoundError:
            pass

        path = pathlib.Path(path)
        tmp_path = pathlib.Path(tmp_path)

        frames1 = read_ase(path, db_format=db_format)
        frames1.st.write_ase_db(tmp_path, db_format=db_format)
        frames2 = read_ase(tmp_path, db_format=db_format)
        assert_equal(frames1, frames2)

    def test_generic_ase(self, FauxBulkDataDB):
        # Just test if we can open it. Should not fail
        with pytest.raises(ValueError) as excinfo:
            frames1 = read_ase(pathlib.Path(FauxBulkDataDB), db_format=None)

    def test_invalid_db_format_raises_valueerror(self, FauxBulkDataDB):
        with pytest.raises(ValueError) as excinfo:
            # Should fail because of invalid db_format of DB
            frames1 = read_ase(pathlib.Path(FauxBulkDataDB), db_format="harlow")
        assert "'db_format' should be one of 'schnet' or 'spainn'" in str(excinfo.value)

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError) as excinfo:
            # Should fail because of file not existing
            read_ase(pathlib.Path("./nowhere.db"), db_format=None)

    def test_invalid_formats(self, FauxBulkDataDB):
        with pytest.raises(ValueError) as excinfo:
            # Should fail because DB is not in correct format
            read_ase(
                pathlib.Path("tutorials/test_data/ase/schnarc_ch2nh2+.db"),
                db_format="spainn",
            )
        # assert "No rows with the appropriate format" in str(excinfo.value)

        with pytest.raises(ValueError) as excinfo:
            # Should fail because DB is not in correct format
            read_ase(
                pathlib.Path("tutorials/test_data/ase/spainn_ch2nh2+.db"),
                db_format="schnet",
            )
        # assert "No rows with the appropriate format" in str(excinfo.value)


if __name__ == '__main__':
    test = TestASEFunctionality()
    bulk_db_path = FauxBulkDataDB()
    # test.test_ase_round_trip()
    test.test_generic_ase(bulk_db_path)
    test.test_invalid_db_format_raises_valueerror(bulk_db_path)
    test.test_invalid_formats(bulk_db_path)
