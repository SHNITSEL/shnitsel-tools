import shnitsel as sh
import shnitsel.xarray
import xarray as xr
import pytest


@pytest.fixture
def traj_butene():
    frames = sh.parse.read_trajs('tutorials/test_data/sharc/traj_butene', kind='sharc')
    return frames


@pytest.fixture
def iconds_butene():
    iconds = sh.parse.sharc_icond.dirs_of_iconds(
        'tutorials/test_data/sharc/iconds_butene'
    )
    return iconds


def test_da_accessors(traj_butene):
    frames = sh.postprocess.ts_to_time(traj_butene)
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


def test_ds_accessors(traj_butene, iconds_butene):
    assert 'ts_to_time' in dir(traj_butene.sh)
    frames = traj_butene.sh.ts_to_time()
    # TODO Add more methods -- by parametrizing.
    assert 'iconds_to_frames' in dir(iconds_butene.sh)