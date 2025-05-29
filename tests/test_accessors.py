import shnitsel as sh
import shnitsel.xarray
import xarray as xr


def test_da_accessors():
    # frames = sh.dynamic.ase.read_ase_db('/nc/SHNITSEL-data/old_CH2NH2.db')
    frames = sh.open_frames('/nc/SHNITSEL_databases/dynamic/A01_ethene_dynamic.nc')
    e_step = frames.energy.sh.sudi()
    assert isinstance(e_step, xr.DataArray)
    xyz = frames.atXYZ.isel(frame=0).sh.to_xyz()
    assert isinstance(xyz, str)
    dihedrals = frames.atXYZ.sh.dihedral(0, 1, 2, 3)
    assert isinstance(dihedrals, xr.DataArray)


if __name__ == '__main__':
    test_da_accessors()