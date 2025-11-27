import inspect
import logging
from typing import Callable, Dict, List



def generate_class_code(classes: Dict[str, List[Callable]]) -> str:
    """
    Generate source code for a class with methods wrapping given functions, preserving
    signatures (including type hints when possible).

    The class names are indicated by the keys in the provided dict, whereas the
    functions to be wrapped in the class methods are denoted by the list in
    the values of the provided dict.

    Parameters
    ----------
        classes dict[str, list[callable]]:
            A dict mapping the names of classes to be generated to the list
            of functions that should be wrapped in methods of that class.

    Returns
    -------
        A string with the generated Python source code.
    """

    # TODO: FIXME: We should account for the import of modules to fail if they are imported from optional packages.\
    # Maybe we should include a flag if the import failed for each individual module and throw an error only if an accessor method
    # using that import is called?

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
        'Callable': 'typing',
        'DataArrayGroupBy': 'xarray.core.groupby',
        'DatasetGroupBy': 'xarray.core.groupby',
        'needs': '._contracts',
        'DAManualAccessor': '._accessors',
        'DSManualAccessor': '._accessors',
        'DatasetOrArray': 'shnitsel.core.midx',
    }
    plain_imports = {
        'xarray as xr',
        'xarray',
        'collections',
        'numpy',
        'numpy.typing as npt',
        'os',
        'pathlib',
        'typing',
        'sklearn',
        'rdkit',
        'os',
        'pathlib'
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
    try:
        import shnitsel as st
        import shnitsel.units as units
        from shnitsel.io.ase.write import write_ase_db
        from shnitsel.core.plot import p3mhelpers
        from shnitsel.core.plot import select
        from shnitsel import bridges
        from shnitsel.core import (
            convenience,
            filtration,
            geom,
            midx,
            ml,
            generic,
            vmd,
            spectra,
            stats,
        )
    except ImportError as e:
        logging.error(
            f"Import of module for generation of accessor classes failed: {e.msg} \n{repr(e)}. \n Please ensure all modules are available."
        )

    da_funcs = [
        # postprocess
        generic.norm,
        generic.subtract_combinations,
        generic.keep_norming,
        stats.calc_ci,
        stats.time_grouped_ci,
        bridges.to_xyz,
        bridges.traj_to_xyz,
        bridges.to_mol,
        bridges.smiles_map,
        bridges.default_mol,
        convenience.pairwise_dists_pca,
        # postprocess converters
        units.convert_energy,
        units.convert_force,
        units.convert_dipole,
        units.convert_length,
        units.convert_time,
        units.convert_nacs,
        units.convert_socs,
        # midx
        midx.mdiff,
        midx.flatten_levels,
        midx.expand_midx,
        midx.assign_levels,
        midx.mgroupby,
        midx.msel,
        midx.sel_trajs,
        midx.sel_trajids,
        # filtration
        filtration.last_time_where,
        # geom
        geom.dihedral,
        geom.angle,
        geom.distance,
        geom.get_bond_lengths,
        geom.get_bond_angles,
        geom.get_bond_torsions,
        geom.get_pyramids,
        geom.get_bats,
        geom.kabsch,
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
        ml.pca,
        ml.lda,
        ml.pls,
    ]

    ds_funcs = [
        # postprocess
        convenience.pca_and_hops,
        convenience.validate,
        spectra.assign_fosc,
        spectra.ds_broaden_gauss,
        stats.get_per_state,
        stats.get_inter_state,
        st.core.populations.calc_pops,
        bridges.default_mol,
        # xrhelpers
        midx.flatten_levels,
        midx.expand_midx,
        midx.assign_levels,
        midx.mgroupby,
        midx.msel,
        midx.sel_trajs,
        midx.unstack_trajs,
        midx.stack_trajs,
        st.io.shnitsel.write_shnitsel_file,
        # parse
        st.io.sharc.parse_initial_conditions.iconds_to_frames,
        # plot
        st.core.spectra.spectra_all_times,
        # filtration
        filtration.energy_filtranda,
        filtration.get_cutoffs,
        filtration.truncate,
        # ase
        write_ase_db,
        # ml
        ml.pls_ds,
    ]

    code = generate_class_code(
        {
            "DataArrayAccessor(DAManualAccessor)": da_funcs,
            "DatasetAccessor(DSManualAccessor)": ds_funcs,
        }
    )
    with open("../shnitsel/_generated_accessors.py", "w") as f:
        print(code, file=f)


if __name__ == "__main__":
    main()
