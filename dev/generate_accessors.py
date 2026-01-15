import logging

from libgenerate import generate_class_code


def main():
    try:
        import shnitsel as st
        from shnitsel import bridges
        from shnitsel import clean
        from shnitsel.clean import common, filter_energy, filter_geo
        from shnitsel import units
        from shnitsel.analyze import (
            generic,
            spectra,
            stats,
            pca,
            lda,
            pls,
            populations,
            hops,
        )
        from shnitsel.data import multi_indices
        import shnitsel.data.helpers as data_helpers
        from shnitsel.geo import geocalc
        from shnitsel.geo.geocalc_ import angles, distances, dihedrals
        from shnitsel.io.ase.write import write_ase_db
        from shnitsel.vis.plot import p3mhelpers
        from shnitsel.vis.plot import select
        from shnitsel.vis import vmd

    except ImportError as e:
        logging.error(
            f"Import of module for generation of accessor classes failed: {e.msg} \n{repr(e)}. \n Please ensure all modules are available."
        )
        raise

    da_funcs = [
        # analyze
        generic.norm,
        generic.subtract_combinations,
        generic.keep_norming,
        generic.pwdists,
        stats.calc_confidence_interval,
        stats.time_grouped_confidence_interval,
        bridges.to_xyz,
        bridges.traj_to_xyz,
        bridges.to_mol,
        bridges.smiles_map,
        bridges.construct_default_mol,
        pca.pca,
        # postprocess converters
        units.convert_energy,
        units.convert_force,
        units.convert_dipole,
        units.convert_length,
        units.convert_time,
        units.convert_nacs,
        units.convert_socs,
        # midx
        multi_indices.mdiff,
        multi_indices.flatten_levels,
        multi_indices.expand_midx,
        multi_indices.assign_levels,
        multi_indices.mgroupby,
        multi_indices.msel,
        multi_indices.sel_trajs,
        multi_indices.sel_trajids,
        # clean
        clean.sanity_check,
        # geom
        dihedrals.dihedral,
        angles.angle,
        distances.distance,
        geocalc.get_distances,
        geocalc.get_angles,
        geocalc.get_dihedrals,
        geocalc.get_pyramidalization,
        geocalc.get_max_chromophor_BLA,
        geocalc.get_bats,
        geocalc.kabsch,
        # geomatch_exact.get_bats_matching,
        # select
        select.FrameSelector,
        select.TrajSelector,
        # p3mhelpers
        p3mhelpers.frame3D,
        p3mhelpers.frames3Dgrid,
        p3mhelpers.traj3D,
        p3mhelpers.trajs3Dgrid,
        # vmd
        vmd.traj_vmd,
        # ml
        pca.pca,
        lda.lda,
        pls.pls,
        # hops
        hops.hops_mask_from_active_state,
        hops.filter_data_at_hops,
        hops.focus_hops,
        hops.assign_hop_time,
    ]

    ds_funcs = [
        # postprocess
        pca.pca_and_hops,
        data_helpers.validate,
        spectra.get_fosc,
        spectra.apply_gauss_broadening,
        stats.get_per_state,
        stats.get_inter_state,
        populations.calc_pops,
        bridges.construct_default_mol,
        # xrhelpers
        multi_indices.flatten_levels,
        multi_indices.expand_midx,
        multi_indices.assign_levels,
        multi_indices.mgroupby,
        multi_indices.msel,
        multi_indices.sel_trajs,
        multi_indices.unstack_trajs,
        multi_indices.stack_trajs,
        st.io.shnitsel.write_shnitsel_file,
        # plot
        spectra.get_spectra,
        # filtration
        filter_energy.calculate_energy_filtranda,
        filter_energy.filter_by_energy,
        clean.sanity_check,
        filter_geo.calculate_bond_length_filtranda,
        filter_geo.filter_by_length,
        common.omit,
        common.truncate,
        common.transect,
        # ase
        write_ase_db,
        # ml
        pls.pls_ds,
        # hops
        hops.hops_mask_from_active_state,
        hops.filter_data_at_hops,
        hops.focus_hops,
        hops.assign_hop_time,
        # select
        select.FrameSelector,
        select.TrajSelector,
    ]

    code = generate_class_code(
        {
            "DataArrayAccessor(DAManualAccessor)": da_funcs,
            "DatasetAccessor(DSManualAccessor)": ds_funcs,
        },
        imports={
            'Union': 'typing',
            'Optional': 'typing',
            'List': 'typing',
            'Dict': 'typing',
            'Hashable': 'typing',
            'Sequence': 'typing',
            'Literal': 'typing',
            'Callable': 'typing',
            'PCAResult' : 'shnitsel.analyze.pca',
            'PopulationStatistics' : 'shnitsel.analyze.populations',
            'DataArrayGroupBy': 'xarray.core.groupby',
            'DatasetGroupBy': 'xarray.core.groupby',
            'needs': '._contracts',
            'DAManualAccessor': '._accessors',
            'DSManualAccessor': '._accessors',
            'DatasetOrArray': 'shnitsel.core.typedefs',
            'nan': 'numpy',
        },
        plain_imports={
            'shnitsel',
            'xarray as xr',
            'xarray',
            'collections',
            'numpy',
            'numpy.typing as npt',
            'numbers',
            'os',
            'pathlib',
            'typing',
            'sklearn',
            'rdkit',
            'os',
            'pathlib',
        },
    )
    with open("../shnitsel/_generated_accessors.py", "w") as f:
        print(code, file=f)


if __name__ == "__main__":
    main()
