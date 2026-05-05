from shnitsel.data.dataset_containers.multi_layered import MultiSeriesLayered
from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked
from shnitsel.data.tree.node import TreeNode
import shnitsel.xarray
import xarray as xr
import pytest


class TestAccessors:
    """Class to test all functions of the shnitsel accessors provided to DataArray and Dataset objects"""

    @pytest.fixture(
        params=[
            ('./tutorials/test_data/shnitsel/traj_I02.nc', 1),
        ]
    )
    def ds_stacked(self, request) -> MultiSeriesStacked:
        from shnitsel.io import read

        path, charge = request.param
        db = read(path).set_charge(charge)
        assert isinstance(db, TreeNode)
        res = db.as_stacked
        return res

    @pytest.fixture(
        params=[
            ('./tutorials/test_data/shnitsel/traj_I02.nc', 1),
        ]
    )
    def ds_layered(self, request) -> MultiSeriesLayered:
        from shnitsel.io import read

        path, charge = request.param
        db = read(path).set_charge(charge)
        assert isinstance(db, TreeNode)
        res = db.as_layered
        return res

    # TODO: Subtests is a feature of pytest 9, which we do not list as requirement
    def test_da_accessors(
        self, ds_stacked: MultiSeriesStacked, ds_layered: MultiSeriesLayered, subtests
    ):
        kws = {
            'norm': dict(),
            'subtract_combinations': dict(),
            'keep_norming': dict(),
            'pwdists': dict(),
            'calc_confidence_interval': dict(),
            'time_grouped_confidence_interval': dict(),
            'to_xyz': dict(),
            'traj_to_xyz': dict(),
            'to_mol': dict(),
            'smiles_map': dict(),
            'default_mol': dict(),
            'pairwise_dists_pca': dict(),
            'convert_energy': dict(to='eV'),
            'convert_force': dict(),
            'convert_dipole': dict(),
            'convert_length': dict(),
            'convert_time': dict(),
            'convert_nacs': dict(),
            'convert_socs': dict(),
            'mdiff': dict(),
            'flatten_levels': dict(),
            'expand_midx': dict(),
            'assign_levels': dict(),
            'mgroupby': dict(),
            'msel': dict(),
            'sel_trajs': dict(),
            'sel_trajids': dict(),
            'true_upto': dict(),
            'dihedral': dict(),
            'angle': dict(),
            'distance': dict(),
            'get_bond_lengths': dict(),
            'get_bond_angles': dict(),
            'get_bond_torsions': dict(),
            'get_pyramids': dict(),
            'get_bats': dict(),
            'kabsch': dict(),
            'FrameSelector': dict(),
            'TrajSelector': dict(),
            'frame3D': dict(),
            'frames3Dgrid': dict(),
            'traj3D': dict(),
            'trajs3Dgrid': dict(),
            'traj_vmd': dict(),
            'pca': dict(),
            'lda': dict(),
            'pls': dict(),
        }
        for multi_key, multi_set in {
            "stacked": ds_stacked,
            "layered": ds_layered,
        }.items():
            for var_name, da in multi_set.data_vars.items():
                for method_name in da.st.suitable:
                    with subtests.test(
                        f"ds({multi_key})['{var_name}'].st.{method_name}",
                        var_name=var_name,
                        method_name=method_name,
                    ):
                        assert hasattr(da.st, method_name)
                        # TODO: FIXME: The generic invocation on all variables leads to issues with some accessors due to missing dimensions etc.
                        # if method_name in kws:
                        #     getattr(da.st, method_name)(**kws[method_name])

    # TODO: Subtests is a feature of pytest 9, which we do not list as requirement
    # TODO: FIXME: Multi-index selectors/helpers are broken
    # @pytest.mark.xfail
    def test_ds_accessors(
        self, ds_stacked: MultiSeriesStacked, ds_layered: MultiSeriesLayered, subtests
    ):
        kws = {
            'pca_and_hops': dict(center_mean=False),
            'flatten_levels': dict(idx_name='frame', levels=['time']),
            'expand_midx': dict(midx_name='frame', level_name='compound', value=''),
            # 'assign_levels': dict(),
            'mgroupby': dict(levels=['atrajectory', 'time']),
            # 'msel': dict(),
            'sel_trajs': dict(trajids_or_mask=[1]),
            'stack_trajs': dict(),
            'unstack_trajs': dict(),
            'energy_filtranda': dict(),
            'sanity_check': dict(),
            'omit': dict(),
            'truncate': dict(),
            'transect': dict(cutoff_time=100),
            'pls_ds': dict(xname='energy', yname='astate'),
            'FrameSelector': dict(data_var='energy', dim='frame'),
            'TrajSelector': dict(),
        }
        blacklist = [
            'msel',  # unclear how to use on generic Dataset
            'write_shnitsel_file',  # avoid side-effects
            'assign_levels',  # acceptable values depend on Dataset contents
            'FrameSelector',  # test Dataset may lack 2d data_var with size 2 dimension
            'TrajSelector',  # cf. FrameSelector
            'truncate', # Cannot be invoked without filtranda set, which we do not have by default
            'transect', # Cannot be invoked without filtranda set, which we do not have by default
        ]
        variant_blacklist = {
            "stacked": [],
            "layered": ["flatten_levels", "expand_midx", "mgroupby"],
        }
        for multi_key, multi_set in {
            "stacked": ds_stacked,
            "layered": ds_layered,
        }.items():
            for method_name in multi_set.st.suitable:
                with subtests.test(
                    f"ds({multi_key}).st.{method_name}", method_name=method_name
                ):
                    assert hasattr(multi_set.st, method_name)
                    if (
                        method_name in kws
                        and method_name not in blacklist
                        and method_name not in variant_blacklist[multi_key]
                    ):
                        print(method_name)
                        getattr(multi_set.st, method_name)(**kws[method_name])
