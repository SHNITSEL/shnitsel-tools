import shnitsel as sh
import shnitsel.xarray


def test_da_accessors():
    # frames = sh.dynamic.ase.read_ase_db('/nc/SHNITSEL-data/old_CH2NH2.db')
    frames = sh.open_frames('/nc/SHNITSEL_databases/dynamic/A01_ethene_dynamic.nc')
    xyz = frames.energy.sh.sudi()
    assert isinstance(xyz, str)


if __name__ == '__main__':
    test_da_accessors()