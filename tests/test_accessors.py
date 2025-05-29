import shnitsel as sh
import shnitsel.xarray
import xarray as xr


def test_da_accessors():
    # frames = sh.dynamic.ase.read_ase_db('/nc/SHNITSEL-data/old_CH2NH2.db')
    frames = sh.open_frames('/nc/SHNITSEL_databases/dynamic/A01_ethene_dynamic.nc')
    e_step = frames.energy.sh.sudi()
    assert isinstance(e_step, xr.DataArray)
    xyz = frames.atXYZ.isel(frame=0).squeeze().sh.to_xyz()
    assert isinstance(xyz, str)
    dihedrals = frames.atXYZ.sh.dihedral(0, 1, 2, 3)
    assert isinstance(dihedrals, xr.DataArray)

    atom_methods = {'pairwise_dists_pca', 'dihedral', 'angle', 'distance'}
    assert atom_methods <= set(dir(frames.atXYZ.sh))
    astate_methods = {'hop_indices', 'trajs_with_hops', 'get_hop_types'}
    assert astate_methods <= set(dir(frames.astate.sh))
    # astatesh = frames.astate.sh


if __name__ == '__main__':
    test_da_accessors()