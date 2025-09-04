import inspect
from textwrap import indent
from typing import get_type_hints


def generate_class_code(class_name: str, functions: list) -> str:
    """\
    Generate source code for a class with methods wrapping given functions, preserving
    signatures (including type hints when possible).

    Parameters
    ----------
        class_name
            Name of the class to generate.
        functions
            List of function objects.

    Returns
    -------
        A string with the generated Python source code.
    """

    def get_ann_str(annotation):
        if inspect.isclass(annotation):
            if annotation.__module__ != "builtins":
                imports[annotation.__name__] = annotation.__module__
            return f"{annotation.__name__}"
        elif annotation is not None:
            return repr(annotation)
        else:
            return ""

    lines = []
    imports = {}
    plain_imports = {'typing', 'xarray'}
    # optional imports
    for func in functions:
        imports[func.__name__] = func.__module__

    lines.append("")

    # class header
    lines.append(f"class {class_name}:")
    if not functions:
        lines.append(indent("pass", "    "))
    else:
        for func in functions:
            name = func.__name__
            module = func.__module__
            sig = inspect.signature(func)
            hints = get_type_hints(func)

            # Build signature string with type hints
            # and arguments for wrapped function
            params = []
            args = []
            for pname, param in sig.parameters.items():
                annotation = hints.get(pname)
                # print(f"# {annotation.__name__} :: {annotation.__module__}")

                if param.default is not inspect.Parameter.empty:
                    pdefault = f"={param.default!r}"
                    adefault = f"={pname}"
                else:
                    pdefault = ""
                    adefault = ""
                ann_str = get_ann_str(annotation)
                if ann_str:
                    ann_str = ": " + ann_str
                params.append(f"{pname}{ann_str}{pdefault}")
                args.append(f"{pname}{adefault}")

            ret_ann = hints.get("return")
            ret_str = get_ann_str(ret_ann)
            if ret_str:
                ret_str = " -> " + ret_str
            param_str = ", ".join(["self"] + params[1:])
            arg_str = ", ".join(["self._obj"] + args[1:])

            # Build method
            method = f"""
    def {name}({param_str}){ret_str}:
        \"\"\"Wrapper for :py:func:`{module}.{name}`.\"\"\"
        return {name}({arg_str})
            """
            lines.append(method)

    import_str = "\n".join(f"import {module}" for module in plain_imports)
    import_str += "\n"
    import_str += "\n".join(
        [f"from {module} import {name}" for name, module in imports.items()]
    )
    # lines.append(f"from {func.__module__} import {func.__name__}")

    return import_str + "\n\n" + "\n".join(lines)


if __name__ == "__main__":
    import shnitsel as st
    from shnitsel.core.plot import p3mhelpers

    da_funcs = [
        # postprocess
        st.postprocess.norm,
        st.postprocess.subtract_combinations,
        st.postprocess.pairwise_dists_pca,
        st.postprocess.sudi,
        st.postprocess.hop_indices,
        st.postprocess.relativize,
        st.postprocess.ts_to_time,
        st.postprocess.keep_norming,
        st.postprocess.calc_ci,
        st.postprocess.time_grouped_ci,
        st.postprocess.to_xyz,
        st.postprocess.traj_to_xyz,
        st.postprocess.dihedral,
        st.postprocess.angle,
        st.postprocess.distance,
        st.postprocess.trajs_with_hops,
        st.postprocess.get_hop_types,
        st.postprocess.to_mol,
        st.postprocess.default_mol,
        # xrhelpers
        st.xrhelpers.flatten_levels,
        st.xrhelpers.expand_midx,
        st.xrhelpers.assign_levels,
        st.xrhelpers.mgroupby,
        st.xrhelpers.msel,
        st.xrhelpers.sel_trajs,
        st.xrhelpers.sel_trajs,
        st.xrhelpers.sel_trajids,
        # filtre_unphysical
        st.core.filter_unphysical.smiles_map,
        # filtre
        st.core.filtre.last_time_where,
        # geom
        st.core.geom.get_bond_lengths,
        st.core.geom.get_bond_angles,
        st.core.geom.get_bond_torsions,
        st.core.geom.get_bats,
        st.core.geom.kabsch,
        # select
        p3mhelpers.frame3D,
        p3mhelpers.frames3Dgrid,
        p3mhelpers.traj3D,
        p3mhelpers.trajs3Dgrid,
        # ml
        st.core.ml.pca,
        st.core.ml.lda,
        st.core.ml.pls,
    ]

    ds_funcs = [
        # postprocess
        st.postprocess.pca_and_hops,
        st.postprocess.validate,
        st.postprocess.ts_to_time,
        st.postprocess.setup_frames,
        st.postprocess.setup_frames,
        st.postprocess.assign_fosc,
        st.postprocess.broaden_gauss,
        st.postprocess.get_per_state,
        st.postprocess.get_inter_state,
        st.postprocess.calc_pops,
        st.postprocess.find_hops,
        st.postprocess.default_mol,
        # xrhelpers
        st.xrhelpers.flatten_levels,
        st.xrhelpers.expand_midx,
        st.xrhelpers.assign_levels,
        st.xrhelpers.mgroupby,
        st.xrhelpers.msel,
        st.xrhelpers.save_frames,
        st.xrhelpers.sel_trajs,
        st.xrhelpers.unstack_trajs,
        st.xrhelpers.stack_trajs,
        # parse
        st.parse.sharc_icond.iconds_to_frames,
        # plot
        st.core.plot.spectra3d.spectra_all_times,
        # filtre
        st.core.filtre.energy_filtranda,
        st.core.filtre.get_cutoffs,
        st.core.filtre.truncate,
        # ase
        st.core.ase.write_ase,
        # ml
        st.core.ml.pls_ds,
    ]

    da_code = generate_class_code("GeneratedDAAcessor", da_funcs)
    ds_code = generate_class_code("GeneratedDSAcessor", ds_funcs)
    with open('../shnitsel/_generated_accessors.py', 'w') as f:
        print(da_code, ds_code, sep="\n", file=f)