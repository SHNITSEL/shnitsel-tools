import inspect


def generate_class_code(classes: dict[str, list[callable]]) -> str:
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
        """Convert annotation to string, handling complex types better."""
        if annotation is None or annotation is inspect._empty:
            return ""

        # Handle string annotations (forward references)
        if isinstance(annotation, str):
            return annotation

        # For simple types, use __name__ if available
        if hasattr(annotation, '__name__'):
            if (
                hasattr(annotation, '__module__')
                and annotation.__module__ != "builtins"
            ):
                module_name = annotation.__module__
                if module_name not in ['typing', 'builtins']:
                    imports[annotation.__name__] = module_name
            return annotation.__name__

        return repr(annotation)

    lines = []
    imports = {
        'Union': 'typing',
        'Optional': 'typing',
        'List': 'typing',
        'Dict': 'typing',
        'Hashable': 'typing',
        'Sequence': 'typing',
        'Literal': 'typing',
        'DataArrayGroupBy': 'xarray.core.groupby',
        'DatasetGroupBy': 'xarray.core.groupby',
        'needs': '._contracts',
        'DAManualAccessor': '._accessors',
        'DSManualAccessor': '._accessors',
    }
    plain_imports = {
        'xarray as xr',
        'xarray',
        'collections',
        'numpy',
        'numpy.typing as npt',
        'typing',
        'sklearn',
    }

    # Collect imports for all functions
    for class_name, functions in classes.items():
        for func in functions:
            imports[func.__name__] = func.__module__

        lines.append("")

        # class header
        lines.append(f"class {class_name}:")
        if not functions:
            lines.append("    pass")
        else:
            lines.append("    _methods = [")
            for func in functions:
                lines.append(f"        {func.__name__!r},")

            lines.append("    ]\n")

            for func in functions:
                name = func.__name__
                module = func.__module__
                sig = inspect.signature(func)

                hints = getattr(func, '__annotations__', {})

                # Build signature string with type hints
                # and arguments for wrapped function
                params = []
                args = []

                for pname, param in sig.parameters.items():
                    annotation = hints.get(pname, param.annotation)
                    ann_str = get_ann_str(annotation)
                    ann_str = f": {ann_str}" if ann_str else ""

                    # Handle different parameter kinds
                    if param.kind == inspect.Parameter.VAR_POSITIONAL:
                        # *args
                        params.append(f"*{pname}{ann_str}")
                        args.append(f"*{pname}")
                    elif param.kind == inspect.Parameter.VAR_KEYWORD:
                        # **kwargs
                        params.append(f"**{pname}{ann_str}")
                        args.append(f"**{pname}")
                    else:
                        # Regular parameters
                        if param.default is not inspect.Parameter.empty:
                            pdefault = f"={param.default!r}"
                            adefault = f"={pname}"
                        else:
                            pdefault = ""
                            adefault = ""
                        params.append(f"{pname}{ann_str}{pdefault}")
                        args.append(f"{pname}{adefault}")

                # Handle return annotation
                ret_ann = hints.get("return")
                ret_str = get_ann_str(ret_ann)
                ret_str = f" -> {ret_str}" if ret_str else ""

                # Build parameter and argument strings
                param_str = ", ".join(["self"] + params[1:])
                arg_str = ", ".join(["self._obj"] + args[1:])

                # Carry over needs decorator
                needs_kws = []
                if hasattr(func, '_needs'):
                    for field in func._needs._fields:
                        val = getattr(func._needs, field)
                        if val is None:
                            continue
                        elif isinstance(val, set):
                            # Generate sets with consistent order
                            reprs = [repr(x) for x in sorted(list(val))]
                            set_str = '{' + ', '.join(reprs) + '}'
                            needs_kws.append(f"{field}={set_str}")
                        else:
                            needs_kws.append(f"{field}={val}")
                    needs_str = ", ".join(needs_kws)
                    needs_str = f"""\
    @needs({needs_str})
"""
                else:
                    needs_str = ""

                # Build method
                method = f"""{needs_str}\
    def {name}({param_str}){ret_str}:
        \"\"\"Wrapper for :py:func:`{module}.{name}`.\"\"\"
        return {name}({arg_str})
"""
                lines.append(method)

    # Generate import statements
    import_lines = []
    for module in sorted(plain_imports):
        import_lines.append(f"import {module}")

    # Group imports by module
    module_imports = {}
    for name, module in imports.items():
        if module not in module_imports:
            module_imports[module] = []
        module_imports[module].append(name)

    for module, names in sorted(module_imports.items()):
        import_lines.append(f"from {module} import {', '.join(sorted(names))}")

    import_str = "\n".join(import_lines)

    return import_str + "\n\n" + "\n".join(lines)

def main():
    import shnitsel as st
    from shnitsel.core.plot import p3mhelpers
    from shnitsel.core.plot import select

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
        st.postprocess.smiles_map,
        st.postprocess.default_mol,
        # postprocess converters
        st.postprocess.convert_energy,
        st.postprocess.convert_forces,
        st.postprocess.convert_dipoles,
        st.postprocess.convert_length,
        # xrhelpers
        st.xrhelpers.flatten_levels,
        st.xrhelpers.expand_midx,
        st.xrhelpers.assign_levels,
        st.xrhelpers.mgroupby,
        st.xrhelpers.msel,
        st.xrhelpers.sel_trajs,
        st.xrhelpers.sel_trajids,
        # filtre
        st.core.filtre.last_time_where,
        # geom
        st.core.geom.get_bond_lengths,
        st.core.geom.get_bond_angles,
        st.core.geom.get_bond_torsions,
        st.core.geom.get_bats,
        st.core.geom.kabsch,
        # select
        select.FrameSelector,
        select.TrajSelector,
        # p3mhelpers
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
        st.postprocess.assign_fosc,
        st.postprocess.ds_broaden_gauss,
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

    code = generate_class_code(
        {
            "DataArrayAccessor(DAManualAccessor)": da_funcs,
            "DatasetAccessor(DSManualAccessor)": ds_funcs,
        }
    )
    with open('../shnitsel/_generated_accessors.py', 'w') as f:
        print(code, file=f)

if __name__ == "__main__":
    main()