from itertools import chain
from logging import info

import xarray as xr

MAX_IDS = 20

class InconsistentAttributeError(ValueError):
    pass


class MultipleCompoundsError(ValueError):
    pass


class MissingValue:
    """Sentinel value for ``tree_to_frames``."""


def _common_coords_attrs(tree, ensure_unique):
    exclude_attrs = {
        'DataTree_Level',
        'misc_input_settings',
        '__original_dataset',
        'trajid',
    }
    coord_names = []
    unique_values = {k: {} for k in iter(ensure_unique)}
    dvauv = {}
    for node in tree.children.values():
        coord_names.extend(
            [k for k in node.attrs if k not in exclude_attrs | ensure_unique]
        )

        for varname, var in chain(node.data_vars.items(), node.coords.items()):
            if varname not in dvauv:
                dvauv[varname] = {}
            for k in var.attrs:
                if k not in dvauv[varname]:
                    dvauv[varname][k] = {}

    return coord_names, unique_values, dvauv


def _collect_values(tree, ensure_unique):
    coord_names, unique_values, dvauv = _common_coords_attrs(tree, ensure_unique)

    datasets = []
    trajids = []
    coords = {k: ('trajid_', []) for k in coord_names}
    for node in tree.children.values():
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
            for var_attr_name in dvauv[var_name]:
                v = var.attrs.get(var_attr_name, MissingValue)
                if v not in dvauv[var_name][var_attr_name]:
                    dvauv[var_name][var_attr_name][v] = []
                dvauv[var_name][var_attr_name][v].append(node.attrs['trajid'])

    return datasets, trajids, coords, unique_values, dvauv

def _get_message(vals):
    message = ""
    for val, ids in vals.items():
        inc = " including" if len(ids) > MAX_IDS else ""
        ids_shown = " ".join([str(x) for x in ids[:MAX_IDS]])
        message += (
            f"  - value {val!r} in {len(ids)} trajectories,{inc} IDs: {ids_shown}\n"
        )
    return message

def _validate(unique_values, dvauv):
    attrs = {}
    messages = ""
    for k, vals in unique_values.items():
        if len(vals) != 1:
            messages += f"- There are {len(vals)} different values for {k}:\n"
            messages += _get_message(vals)
        elif next(iter(vals)) is MissingValue:
            messages += f"- The attribute {k} is missing in all trajectories."
        else:
            attrs[k] = next(iter(vals))

    dv_attrs = {}
    for var_name, var_attr_data in dvauv.items():
        for var_attr_name, vals in var_attr_data.items():
            if len(vals) != 1:
                messages += f"- There are {len(vals)} different values for {var_attr_name} in {var_name}:\n"
                messages += _get_message(vals)
            elif next(iter(vals)) is MissingValue:
                messages += f"- The attribute {var_attr_name} in {var_name} is missing in all trajectories."
            else:
                dv_attrs.setdefault(var_name, {})[var_attr_name] = next(iter(vals))

    if messages:
        raise InconsistentAttributeError("The following issues arose --\n" + messages)

    return attrs, dv_attrs


def _concat(tree, ensure_unique):
    per_traj_dim_name = 'trajid_'
    datasets, trajids, coords, unique_values, dvauv = _collect_values(
        tree, ensure_unique
    )

    res = xr.concat(datasets, 'frame').assign_coords(
        trajid_=(per_traj_dim_name, trajids)
    )

    attrs, dv_attrs = _validate(unique_values, dvauv)
    for k, v in dv_attrs.items():
        res[k].attrs.update(v)

    return res.assign_coords(coords).assign_attrs(attrs)


def tree_to_frames(tree, allow_inconsistent: set | None = None) -> xr.Dataset:
    """Transform a DataTree into a single stacked Dataset.

    Parameters
    ----------
    tree
        The py:class:`xarray.DataTree` to transform: a node which is either

            - at the CompoundGroup level, with TrajectoryData children
            - at the ShnitselDBRoot level, containing a single CompoundGroup

    allow_inconsistent, optional
        A list specifying attributes that should *not* be checked
        for consistency, whereas they normally would be.

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

    MultipleCompoundsError
        If ``tree`` is at the ShnitselDBRoot level and has multiple children.


    Examples
    --------
    >>> frames = tree_to_frames(dt['/unknown'], allow_inconsistent={'delta_t'})
    """
    # TODO: Is there a guarantee that the tree structure follows the hierarchy
    # ShnitselDBRoot -> CompoundGroup -> TrajectoryData?
    # Would it be better to check whether the children of `tree` have the 'trajid' attr set?
    if len(tree.children) == 0:
        raise ValueError("The root node of 'tree' has no children.")
    elif (
        tree.attrs.get('DataTree_Level') == 'ShnitselDBRoot'
        or 'trajid' not in next(iter(tree.children.values())).attrs
    ):
        # Assume we have been given a tree containing a CompoundGroup level
        compound_names = list(tree.children)
        if len(compound_names) == 1:
            target = compound_names[0]
            info(f"Converting the only CompoundGroup, named '{target}'")
            tree = tree.children[target]
        else:
            raise ValueError(
                "'tree' contains trajectories for multiple CompoundGroups, namely "
                f"{compound_names}. Please extract a single CompoundGroup first, "
                "e.g. instead of `tree_to_frames(dt)`, try "
                f"`tree_to_frames(dt['{compound_names[0]}'])`."
            )

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

    return _concat(tree, ensure_unique)
