from itertools import chain

import xarray as xr


class InconsistentAttributeError(ValueError):
    pass


class MissingValue:
    "Sentinel value for ``tree_to_frames``"

    pass


def tree_to_frames(tree, allow_inconsistent: set | None = None) -> xr.Dataset:
    """Transforms a DataTree into a single stacked Dataset

    Parameters
    ----------
    tree
        The DataTree to transform
    allow_inconsistent, optional
        A list specifying attributes that should *not* be checked
        for consistency, whereas they normally would be. By default None

    Returns
    -------
        A single Dataset with trajectories stacked along a dimension ``frame``;
        attributes required to be consistent across trajectories remain attributes;
        attributes permitted to vary across trajectories become coordinates;
        other Dataset-level attributes are ignored and omitted.
        Variable-level attributes are checked for consistency and propagated to the
        result.

    Raises
    ------
    InconsistentAttributeError
        If any of those attributes required to be unique across trajectories
        violate this condition, or if any of them are missing in all trajectories
        (in which case their value is consistent but invalid); this error can be suppressed
        by specifying the appropriate attribute names in the ``allow_inconsistent`` parameter.
        Note that suppression only works for Dataset-level attributes; inconsistency
        amongst Variable-level attributes always raises.

    Examples
    --------
    >>> frames = tree_to_frames(dt['/unknown'], allow_inconsistent={'delta_t'})
    """
    per_traj_dim_name = 'trajid_'
    exclude_attrs = {
        'DataTree_Level',
        'misc_input_settings',
        '__original_dataset',
        'trajid',
    }
    ensure_unique = {
        'input_type',
        'input_format_version',
        'delta_t',
        'num_singlets',
        'num_doublets',
        'num_triplets',
    }
    # vars_ensure_unique = {'long_name', 'unitdim', 'units', 'original_units'}
    if allow_inconsistent is not None:
        ensure_unique = ensure_unique.difference(allow_inconsistent)
        # vars_ensure_unique = vars_ensure_unique.difference(allow_inconsistent)
    datasets = []
    trajids = []
    coord_names = []
    unique_values = {k: {} for k in iter(ensure_unique)}
    var_attrs_unique_values = {}
    for node in tree.children.values():
        for k in node.attrs:
            if k not in exclude_attrs | ensure_unique:
                coord_names.append(k)

        for varname, var in chain(node.data_vars.items(), node.coords.items()):
            if varname not in var_attrs_unique_values:
                var_attrs_unique_values[varname] = {}
            for k in var.attrs:
                if k not in var_attrs_unique_values[varname]:
                    var_attrs_unique_values[varname][k] = {}

    coords = {k: ('trajid_', []) for k in coord_names}
    for i, node in enumerate(tree.children.values()):
        ds = (
            node.to_dataset()
            .expand_dims(trajid=[node.attrs['trajid']])
            .stack(frame=['trajid', 'time'])
            .drop_attrs()
        )
        datasets.append(ds)
        trajids.append(node.attrs['trajid'])
        for k in coords:
            coords[k][1].append(node.attrs.get(k, MissingValue))
        for k in iter(ensure_unique):
            v = node.attrs.get(k, MissingValue)
            if v not in unique_values[k]:
                unique_values[k][v] = []
            unique_values[k][v].append(node.attrs['trajid'])

        for var_name, var in chain(node.data_vars.items(), node.coords.items()):
            for var_attr_name in var_attrs_unique_values[var_name]:
                v = var.attrs.get(var_attr_name, MissingValue)
                if v not in var_attrs_unique_values[var_name][var_attr_name]:
                    var_attrs_unique_values[var_name][var_attr_name][v] = []
                var_attrs_unique_values[var_name][var_attr_name][v].append(
                    node.attrs['trajid']
                )

    attrs = {}
    messages = ""
    for k, vals in unique_values.items():
        if len(vals) != 1:
            messages += f"- There are {len(vals)} different values for {k}:\n"
            for val, ids in vals.items():
                messages += f"  - {k} = {val} in {len(ids)} trajectories, "
                if len(ids) < 20:
                    messages += "IDs: " + " ".join([str(x) for x in ids])
                else:
                    messages += "including IDs: " + " ".join([str(x) for x in ids[:20]])
                messages += "\n"
        elif list(vals)[0] is MissingValue:
            messages += f"- The attribute {k} is missing in all trajectories."
        else:
            attrs[k] = list(vals)[0]

    res = xr.concat(datasets, 'frame').assign_coords(
        trajid_=(per_traj_dim_name, trajids)
    )

    for var_name, var_attr_data in var_attrs_unique_values.items():
        for var_attr_name, vals in var_attr_data.items():
            if len(vals) != 1:
                messages += f"- There are {len(vals)} different values for {var_attr_name} in {var_name}:\n"
                for val, ids in vals.items():
                    messages += f"  - {k} = {val} in {len(ids)} trajectories, "
                    if len(ids) < 20:
                        messages += "IDs: " + " ".join([str(x) for x in ids])
                    else:
                        messages += "including IDs: " + " ".join(
                            [str(x) for x in ids[:20]]
                        )
                    messages += "\n"
            elif list(vals)[0] is MissingValue:
                messages += f"- The attribute {var_attr_name} in {var_name} is missing in all trajectories."
            else:
                if var_name in res.coords:
                    res.coords[var_name].attrs[var_attr_name] = list(vals)[0]
                else:
                    res.data_vars[var_name].attrs[var_attr_name] = list(vals)[0]

    if messages:
        raise InconsistentAttributeError("The following issues arose --\n" + messages)

    return res.assign_coords(coords).assign_attrs(attrs)
