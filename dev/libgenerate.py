import inspect
from typing import Callable, Dict, List


def generate_class_code(
    classes: Dict[str, List[Callable]], imports, plain_imports, name_overrides
) -> str:
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

    # Collect imports for all functions
    for class_name, functions in classes.items():
        for func in functions:
            imports[func.__name__] = func.__module__
            func.__name__ = name_overrides.get(func.__name__, func.__name__)

        lines.append("")

        # class header
        lines.append(f"class {class_name}:")
        if not functions:
            lines.append("    pass")
        else:
            lines.append("    _methods = [")
            lines.extend([f"        {func.__name__!r}," for func in functions])
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
    import_lines = [f"import {module}" for module in sorted(plain_imports)]

    # Group imports by module
    module_imports = {}
    for name, module in imports.items():
        if module not in module_imports:
            module_imports[module] = []
        module_imports[module].append(name)

    for module, names in sorted(module_imports.items()):
        import_lines.append(f"from {module} import {', '.join(sorted(names))}")

    import_str = "\n".join(import_lines)

    # Rename functions
    rename_lines = [
        f"{new_name} = {old_name}\n{new_name}.__name__ = {new_name!r}"
        for old_name, new_name in name_overrides.items()
    ]
    rename_str = '\n'.join(rename_lines)

    return import_str + "\n\n" + rename_str + "\n\n" + "\n".join(lines)
