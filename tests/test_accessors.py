import shnitsel.xarray
import xarray as xr
import pytest
from shnitsel.io import read
from shnitsel.data.tree import tree_to_frames


class TestAccessors:
    """Class to test all functions of the shnitsel accessors provided to DataArray and DataSet objects"""

    @pytest.fixture
    def traj_butene(self):
        db = read('tutorials/test_data/sharc/traj_butene', kind='sharc')
        assert db is not None
        frames = tree_to_frames(db)
        return frames

    @pytest.fixture
    def iconds_butene(self):
        iconds = read('tutorials/test_data/sharc/iconds_butene')
        return iconds

    def test_da_accessors(self, traj_butene):
        frames = traj_butene
        e_step = frames.energy.st.mdiff()
        assert isinstance(e_step, xr.DataArray)
        xyz = frames.atXYZ.isel(frame=0).squeeze().st.to_xyz()
        assert isinstance(xyz, str)
        dihedrals = frames.atXYZ.st.dihedral(0, 1, 2, 3)
        assert isinstance(dihedrals, xr.DataArray)

        atom_methods = {'pairwise_dists_pca', 'dihedral', 'angle', 'distance'}
        assert atom_methods <= set(dir(frames.atXYZ.st))
        astate_methods = {'hop_indices', 'trajs_with_hops', 'get_hop_types'}
        assert astate_methods <= set(dir(frames.astate.st))
        # astatesh = frames.astate.sh

    def test_ds_accessors(self, traj_butene, iconds_butene):
        assert 'ts_to_time' in dir(traj_butene.st)
        frames = traj_butene.st.ts_to_time()
        # TODO Add more methods -- by parametrizing.
        assert 'iconds_to_frames' in dir(iconds_butene.st)
