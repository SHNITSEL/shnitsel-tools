from itertools import combinations, permutations
import json
import logging
import os
import pathlib
from typing import TYPE_CHECKING, Any, Literal

from ase.db import connect
import numpy as np
import pandas as pd
import xarray as xr
import random

from shnitsel.bridges import construct_default_mol, to_mol
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.compound import CompoundGroup, CompoundInfo
from shnitsel.data.tree.data_group import DataGroup
from shnitsel.data.tree.data_leaf import DataLeaf
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDBRoot
from shnitsel.io.shared.helpers import LoadingParameters
from shnitsel.io.shared.trajectory_finalization import normalize_dataset
from shnitsel.io.shared.trajectory_setup import (
    RequiredTrajectorySettings,
    assign_required_settings,
)
from shnitsel.io.shared.variable_flagging import mark_variable_assigned
from shnitsel.units.defaults import get_default_input_attributes
from shnitsel.data.atom_helpers import get_atom_number_from_symbol

dummy_leading_dim: str = "leading_dim_unknown"
multi_level_prefix: str = "_MultiIndex_levels_for_"

# FIXME: Add other wrapped formats and add tree support for ASE

SUPPORTED_COORD_KEYS = {
    "direction",
    "atNames",
    "atNums",
    "state_names",
    "state_types",
    "state_charges",
    "astate",
    "sdiag",
    "time",
    "atrajectory",
    "trajid",
    "from",
    "to",
    "full_statecomb_from",
    "full_statecomb_to",
    "trajectory",
    "trajid_",
    "charge",
}


def complete_shapes_guesses_from_variables(
    tmp_vars: dict[str, np.ndarray],
    meta_shape_vars: dict[str, list[str]],
    meta_shape_coords: dict[str, list[str]],
    db_format_meta: Literal['spainn', 'schnet'] | None,
    leading_dim_name: Literal['frame', 'time'] | str,
    leading_dim_target_guess: Literal['frame', 'time'] | str | None = None,
) -> tuple[
    dict[str, list[str]],
    dict[str, list[str]],
    Literal['spainn', 'schnet'],
    Literal['frame', 'time'] | str,
]:
    """Helper function to guess the format from variable shapes if no metadata could be found.

    _extended_summary_

    Parameters
    ---------
    tmp_vars: dict[str, np.ndarray]
        The temporary arrays of the variables read for this dataset
    meta_shape_vars: dict[str, list[str]]
        The shapes per frame of the variable arrays read for this dataset
    meta_shape_coords: dict[str, list[str]],
        The variable shape s
    db_format_meta: Literal['spainn', 'schnet'] | None,
    leading_dim_name: Literal['frame', 'time'] | str,
    leading_dim_target_guess: Literal['frame', 'time'] | str | None = None,

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    ValueError
        _description_
    ValueError
        _description_
    ValueError
        _description_
    """
    var_shape_dict: dict[str, list[str]] = dict(meta_shape_vars)
    coord_shape_dict: dict[str, list[str]] = dict(meta_shape_coords)
    leading_dim_target: Literal['frame', 'time'] | str = (
        leading_dim_target_guess or 'frame'
    )
    db_format_res: Literal['spainn', 'schnet'] | None = db_format_meta

    if 'energy' in tmp_vars:
        if len(tmp_vars['energy'].shape) == 3:
            db_format_res = 'spainn'
        elif len(tmp_vars['energy'].shape) == 2:
            db_format_res = 'schnet'

    schnet_shapes: dict[str, list[str]] = {
        'atXYZ': [leading_dim_name, 'atom', 'direction'],
        'velocities': [leading_dim_name, 'atom', 'direction'],
        'energy': [leading_dim_name, 'state'],
        'e_kin': [leading_dim_name],
        'forces': [leading_dim_name, 'state', 'atom', 'direction'],
        'nacs': [leading_dim_name, 'statecomb', 'atom', 'direction'],
        'smooth_nacs': [leading_dim_name, 'statecomb', 'atom', 'direction'],
        'socs': [leading_dim_name, 'full_statecomb'],
        'dipoles': [leading_dim_name, 'state_or_statecomb', 'direction'],
        "phases": [leading_dim_name, "state"],
    }
    spainn_shapes: dict[str, list[str]] = {
        'atXYZ': [leading_dim_name, 'atom', 'direction'],
        'velocities': [leading_dim_name, 'atom', 'direction'],
        # Note the extra dim, removed later
        'energy': [leading_dim_name, 'tmp', 'state'],
        'e_kin': [leading_dim_name],
        'forces': [leading_dim_name, 'atom', 'state', 'direction'],
        'nacs': [leading_dim_name, 'atom', 'statecomb', 'direction'],
        'smooth_nacs': [leading_dim_name, 'atom', 'statecomb', 'direction'],
        'socs': [leading_dim_name, 'full_statecomb'],
        'dipoles': [leading_dim_name, 'state_or_statecomb', 'direction'],
        "phases": [leading_dim_name, "state"],
    }

    if db_format_res is None:
        if "atXYZ" in tmp_vars:
            num_atoms = tmp_vars["atXYZ"].shape[-2]
        elif "velocity" in tmp_vars:
            num_atoms = tmp_vars["atXYZ"].shape[-2]
        elif "atNames" in tmp_vars:
            num_atoms = tmp_vars["atNames"].shape[-1]
        elif db_format_res:
            num_atoms = -1
        else:
            raise ValueError(
                "No indications for correct format available. Cannot parse ASE db."
            )

        # Try to find the correct format based on the atom coordinate and the number of atoms
        if num_atoms >= 0:
            is_schnet_compatible = True
            is_spainn_compatible = True
            for key, shapes in schnet_shapes.items():
                if key in tmp_vars:
                    try:
                        atom_index = shapes.index("atom")
                        if tmp_vars[key].shape[atom_index] != num_atoms:
                            is_schnet_compatible = False
                            break
                    except:
                        pass
            for key, shapes in spainn_shapes.items():
                if key in tmp_vars:
                    try:
                        atom_index = shapes.index("atom")
                        if tmp_vars[key].shape[atom_index] != num_atoms:
                            is_spainn_compatible = False
                            break
                    except:
                        pass

            db_format_res = (
                "schnet"
                if is_schnet_compatible
                else "spainn"
                if is_spainn_compatible
                else None
            )

    if db_format_res is None:
        raise ValueError("Could not discern db format from data input")
    if leading_dim_target is None:
        raise ValueError("Could not discern leading dimension from ase data input.")

    var_shape_dict = schnet_shapes if db_format_res == 'schnet' else spainn_shapes

    return var_shape_dict, coord_shape_dict, db_format_res, leading_dim_target


def shapes_from_metadata(
    db_meta: dict, db_format: Literal['spainn', 'schnet'] | None = None
) -> tuple[
    dict[str, list[str]],
    dict[str, list[str]],
    Literal['frame', 'time'] | str,
    Literal['spainn', 'schnet'] | None,
]:
    """Function to assign shapes based on the chosen db_format and potential information in the metadata of a database.

    If conflicting information on the db_format is provided and present in the database, en error will be raised.

    Parameters
    ----------
    db_meta : dict
        The metadata dict of an ASE database.
    db_format : Literal['spainn', 'schnet'] | None, optional
        The requested format of the database. Defaults to None.

    Returns
    -------
    dict[str, list[str]]
        Dict of data_var shapes
    dict[str, list[str]]
        Dict of coordinate shapes
    str
        The name of the leading dimension. Should be `frame` or `time`, but can be `leading_dim_unknown` if unknown

    Raises
    ------
    ValueError
        If a db_format of database was requested that conflicts with the format of the database.
    """

    if "__shnitsel_meta" in db_meta:
        db_meta["__shnitsel_meta"] = json.loads(db_meta["__shnitsel_meta"])
        shnitsel_meta = db_meta["__shnitsel_meta"]
    else:
        shnitsel_meta = {}

    if 'shnitsel_leading_dim' in shnitsel_meta:
        leading_dim_name = shnitsel_meta['shnitsel_leading_dim']
    else:
        leading_dim_name = dummy_leading_dim

    schnet_shapes: dict[str, list[str]] = {
        'atXYZ': [leading_dim_name, 'atom', 'direction'],
        'velocities': [leading_dim_name, 'atom', 'direction'],
        'energy': [leading_dim_name, 'state'],
        'e_kin': [leading_dim_name],
        'forces': [leading_dim_name, 'state', 'atom', 'direction'],
        'nacs': [leading_dim_name, 'statecomb', 'atom', 'direction'],
        'smooth_nacs': [leading_dim_name, 'statecomb', 'atom', 'direction'],
        'socs': [leading_dim_name, 'full_statecomb'],
        'dipoles': [leading_dim_name, 'state_or_statecomb', 'direction'],
        "phases": [leading_dim_name, "state"],
    }

    spainn_shapes: dict[str, list[str]] = {
        'atXYZ': [leading_dim_name, 'atom', 'direction'],
        'velocities': [leading_dim_name, 'atom', 'direction'],
        # Note the extra dim, removed later
        'energy': [leading_dim_name, 'tmp', 'state'],
        'e_kin': [leading_dim_name],
        'forces': [leading_dim_name, 'atom', 'state', 'direction'],
        'nacs': [leading_dim_name, 'atom', 'statecomb', 'direction'],
        'smooth_nacs': [leading_dim_name, 'atom', 'statecomb', 'direction'],
        'socs': [leading_dim_name, 'full_statecomb'],
        'dipoles': [leading_dim_name, 'state_or_statecomb', 'direction'],
        "phases": [leading_dim_name, "state"],
    }

    coord_shapes = {
        "direction": ["direction"],
        "atNames": ["atom"],
        "atNums": ["atom"],
        "state_names": ["state"],
        "state_types": ["state"],
        "state_charges": ["state"],
        "astate": [leading_dim_name],
        "sdiag": [leading_dim_name],
        "time": [leading_dim_name],
        "atrajectory": [leading_dim_name],  # This only exists for frame dimensions
        "trajid": [leading_dim_name],  # This only exists for frame dimensions
        "from": ["statecomb"],
        "to": ["statecomb"],
        "full_statecomb_from": ["full_statecomb"],
        "full_statecomb_to": ["full_statecomb"],
        "trajectory": [leading_dim_name],  # This only exists for frame dimensions
        "trajid_": [leading_dim_name],  # This only exists for frame dimensions
        "charge": [leading_dim_name],  # This only exists for frame dimensions
        # "statecomb": ["statecomb"],
        # "full_statecomb": ["full_statecomb"]
    }

    if "db_format" in shnitsel_meta:
        meta_format = shnitsel_meta["db_format"]
        if meta_format not in ["schnet", "spainn"]:
            raise ValueError(
                f"Database is of unsupported format: {meta_format}. Only `schnet` and `spainn` are supported."
            )

        if db_format is None:
            db_format = meta_format
            logging.info(
                "Automatically detected format: %s" % db_format
            )

        if meta_format != db_format:
            raise ValueError(
                f"Database is of format: {meta_format} instead of requested format {db_format}."
            )
    shapes: dict[str, list[str]]
    # Determine basis shapes based on the format
    if db_format == 'schnet':
        shapes = schnet_shapes
    elif db_format == 'spainn':
        shapes = spainn_shapes
    elif db_format is None:
        shapes = {}
        logging.warning(
            "Correct format could not be extracted from the database metadata. No dimension names assigned"
        )
    else:
        raise ValueError(
            f"'db_format' should be one of 'schnet' or 'spainn', not '{db_format}'."
        )

    # Read further shape data from the database
    if "var_meta" in shnitsel_meta:
        variable_metadata = shnitsel_meta["var_meta"]
        for varname, vardict in variable_metadata.items():
            if "dims" in vardict:
                shapes[varname] = vardict["dims"]

    if "coords" in shnitsel_meta:
        coord_metadata = shnitsel_meta["coords"]
        for coordname, coorddict in coord_metadata.items():
            if "dims" in coorddict:
                coord_shapes[coordname] = coorddict["dims"]

    return shapes, coord_shapes, leading_dim_name, db_format


def apply_dataset_meta_from_db_metadata(
    dataset: xr.Dataset,
    db_meta: dict,
    default_attrs: dict,
) -> xr.Dataset:
    """Apply attributes from db metadata and perform some validation checks on the result.

    Loads remaining missing coordinate variables from db metadata if available.
    Checks size of resulting dimensions if specified in db metadata.
    Further initializes the multi indices if specified in the metadata.

    Parameters
    ----------
    dataset : xr.Dataset
        Trajectory/Frames dataset parsed from ASE db
    db_meta : dict
        Metadata from the trajectory db file
    default_attrs : dict
        Attributes to apply to variables by default


    Returns
    -------
    xr.Dataset
        Dataset with attributes set from from db metadata and dimension sizes asserted
    """
    if "__shnitsel_meta" in db_meta:
        if isinstance(db_meta["__shnitsel_meta"], str):
            db_meta["__shnitsel_meta"] = json.loads(db_meta["__shnitsel_meta"])
        shnitsel_meta = db_meta["__shnitsel_meta"]
    else:
        shnitsel_meta = {}

    # Restore missing coordinates
    if "coords" in shnitsel_meta:
        coords_data = shnitsel_meta["coords"]
        for coordname, coorddict in coords_data.items():
            if coordname not in dataset.coords:
                # We do not want to attempt to overwrite an existing coordinate that may be different
                # due to a split during deserialization
                skip_index = False
                for dim in coorddict['dims']:
                    # Don't introduce new dimensions
                    if dim not in dataset.dims or dataset.sizes[dim] != len(
                        coorddict["values"]
                    ):
                        skip_index = True
                        break
                if skip_index:
                    continue

                dataset = dataset.assign_coords(
                    {
                        coordname: (
                            coorddict["dims"],
                            np.array(coorddict["values"]),
                        )
                    }
                )
                mark_variable_assigned(dataset[coordname])

    # Potentially reconstruct multiindex levels
    if (
        "_MultiIndex_levels_from_attrs" in shnitsel_meta
        and shnitsel_meta["_MultiIndex_levels_from_attrs"] == 1
    ):
        for k, v in shnitsel_meta["__multi_indices"].items():
            if str(k).startswith(multi_level_prefix):
                index_name = str(k)[len(multi_level_prefix) :]
                index_levels = v["level_names"]
                index_tuples = v["index_tuples"]
                index_tuples = [tuple(x) for x in index_tuples]

                index_len = len(index_tuples)

                skip_index = False

                if index_name in dataset.coords:
                    # and dataset.coords[index_name].size != len(index_tuples):
                    # We do not want to attempt to overwrite an existing coordinate that may have been split during deserialization
                    skip_index = True
                    continue

                if index_name in dataset.dims:
                    # If the multi-index dimension already exists, don't try and overwrite existing levels.
                    for level in index_levels:
                        # if level in dataset.coords:
                        #     print(
                        #         f"{level=} {dataset.coords[level].size=} != {index_len=}"
                        #     )
                        if (
                            level in dataset.coords
                            and dataset.coords[level].size != index_len
                        ):
                            skip_index = True
                            break

                if skip_index:
                    continue

                # print(f"{index_name=} not in {dataset.coords.keys()}")

                # Stack the existing dimensions instead of setting an xindex

                # tuples = list(
                #     zip(*[dataset.coords[level].values for level in index_levels])
                # )
                # print(index_name, ":\t", tuples)

                multi_coords = xr.Coordinates.from_pandas_multiindex(
                    pd.MultiIndex.from_tuples(
                        index_tuples,
                        names=index_levels,
                    ),
                    dim=index_name,
                )

                dataset = dataset.assign_coords(multi_coords)

                mark_variable_assigned(dataset[index_name])
                for level in index_levels:
                    mark_variable_assigned(dataset[level])

                # dataset =dataset.stack({index_name: index_levels})
                # dataset = dataset.set_xindex(index_levels)

    # Fill in missing frame/time coordinates
    if "frame" in dataset.dims:
        # Add dummy frame coordinate values treating all entries as single frames
        if "frame" not in dataset.coords:
            frame_vals = np.arange(0, dataset.sizes["frame"], 1)
            dataset = dataset.assign_coords(
                {
                    "frame": (["frame"], frame_vals, default_attrs.get("frame", None)),
                }
            )

        if "atrajectory" not in dataset.coords:
            # Set the atrajectory property to the same id across all frames.
            # The split should have happened in input parsing

            id_vals = np.full([dataset.sizes["frame"]], random.randint(0, 1000000))
            dataset = dataset.assign_coords(
                {
                    "atrajectory": (
                        ["frame"],
                        id_vals,
                        default_attrs.get("atrajectory", None),
                    ),
                }
            )
        if "time" not in dataset.coords:
            # Do not add dummy `time` coordinate.
            pass
    elif "time" in dataset.dims:
        # Fill in missing time coordinate with dummy values if no frame is set as dimension
        if "time" not in dataset:
            time_vals = np.arange(0, dataset.sizes["time"], 1) * (
                dataset.attrs.get(
                    "delta_t", dataset.delta_t if "delta_t" in dataset.coords else 1.0
                )
            )
            dataset = dataset.assign_coords(
                {
                    "time": (["time"], time_vals, default_attrs.get("time", None)),
                }
            )
    else:
        raise ValueError(
            f"Neither `frame` nor `time` dimension generated. Indicates that no data could be read. Available dimensions: `{dataset.sizes.keys()}`. Available coordinates: `{dataset.coords.keys()}`"
        )

    # Apply variable metadata where available
    if "var_meta" in shnitsel_meta:
        vars_dict = shnitsel_meta["var_meta"]
        for varname, vardict in vars_dict.items():
            if varname in dataset:
                if "attrs" in vardict:
                    var_attrs = vardict["attrs"]
                    if varname == "dipoles" and (
                        "dip_perm" in dataset or "dip_trans" in dataset
                    ):
                        # Dipoles should have been split back up and the names should be updated accordingly
                        if "dip_perm" in dataset:
                            dataset["dip_perm"].attrs.update(var_attrs)
                            if "dip_perm" in default_attrs:
                                dataset["dip_perm"]["long_name"] = default_attrs[
                                    "dip_perm"
                                ]["long_name"]
                        if "dip_trans" in dataset:
                            dataset["dip_trans"].attrs.update(var_attrs)
                            if "dip_trans" in default_attrs:
                                dataset["dip_trans"]["long_name"] = default_attrs[
                                    "dip_trans"
                                ]["long_name"]
                    else:
                        dataset[varname].attrs.update(var_attrs)

    if "_distance_unit" in db_meta:
        if "atXYZ" in dataset:  # and "units" not in dataset["atXYZ"].attrs:
            dataset["atXYZ"].attrs["units"] = db_meta["_distance_unit"]

    if "_property_unit_dict" in db_meta:
        unit_dict = db_meta["_property_unit_dict"]

        for varname, unit in unit_dict.items():
            if varname == "dipoles":
                if "dip_perm" in dataset and "unit" not in dataset["dip_perm"].attrs:
                    dataset["dip_perm"].attrs["unit"] = unit
                if "dip_trans" in dataset and "unit" not in dataset["dip_trans"].attrs:
                    dataset["dip_trans"].attrs["unit"] = unit
            else:
                if varname in dataset and "unit" not in dataset[varname].attrs:
                    dataset[varname].attrs["unit"] = unit

    # print(dataset["time"])
    # print(dataset["trajid"])

    delta_t = (
        float(dataset.attrs["delta_t"])
        if "delta_t" in dataset.attrs
        else float(dataset.delta_t)
        if "delta_t" in dataset.coords
        else None
    )
    if delta_t is None:
        # Try and extract from time info
        if "time" in dataset:
            # Sort times
            diff_t = list(set(dataset["time"].values))
            diff_t = sorted(diff_t)
            if len(diff_t) > 1:
                # If we have multiple times, calculate successive distances.
                dt_opts = np.diff(diff_t)
                # If the dt steps are within a tolerance, use that
                if abs(np.min(dt_opts) - np.min(dt_opts)) < 1e-2:
                    delta_t = float(np.mean(dt_opts))
                else:
                    # Make a best-effort guess
                    delta_t = float(np.min(dt_opts))
            else:
                delta_t = 0
        else:
            delta_t = -1

    num_singlets = (
        dataset.attrs["num_singlets"]
        if "num_singlets" in dataset.attrs
        else db_meta["n_singlets"]
        if "n_singlets" in db_meta
        else 0
    )
    num_doublets = (
        dataset.attrs["num_doublets"]
        if "num_doublets" in dataset.attrs
        else db_meta["n_doublets"]
        if "n_doublets" in db_meta
        else 0
    )
    num_triplets = (
        dataset.attrs["num_triplets"]
        if "num_triplets" in dataset.attrs
        else db_meta["n_triplets"]
        if "n_triplets" in db_meta
        else 0
    )

    # miscallaneous properties:
    extract_settings = RequiredTrajectorySettings(
        t_max=dataset.attrs["t_max"]
        if "t_max" in dataset.attrs
        else float(dataset.t_max)
        if "t_max" in dataset.coords
        else dataset.time.max()
        if "time" in dataset.coords
        else -1,
        delta_t=float(delta_t),
        max_ts=dataset.attrs["max_ts"]
        if "max_ts" in dataset.attrs
        else int(dataset.max_ts)
        if "max_ts" in dataset.coords
        else (
            dataset.sizes["time"]
            if "time" in dataset.sizes
            else dataset.sizes["frame"]
            if "frame" in dataset.sizes
            else 0
        ),
        completed=dataset.attrs["completed"] if "completed" in dataset.attrs else False,
        input_format=dataset.attrs["input_format"]
        if "input_format" in dataset.attrs
        else "ase",
        input_type=dataset.attrs["input_type"]
        if "input_type" in dataset.attrs
        else "unknown",
        input_format_version=dataset.attrs["input_format_version"]
        if "input_format_version" in dataset.attrs
        else "unknown",
        num_singlets=num_singlets,
        num_doublets=num_doublets,
        num_triplets=num_triplets,
    )

    dataset = assign_required_settings(dataset, extract_settings)

    # Fix derived coordinates if they are missing
    if "state" in dataset.dims:
        # Fix state coordinates if they are missing
        if "state_names" not in dataset or "state_types" not in dataset:
            if "states" in db_meta:
                state_name_data = np.array(str(db_meta["states"]).split(), dtype='U8')
                state_type_data = np.array(
                    [
                        1
                        if x.startswith("S")
                        else 2
                        if x.startswith("D")
                        else 3
                        if x.startswith("T")
                        else -1
                        for x in state_name_data
                    ]
                )

                dataset = dataset.assign_coords(
                    state_types=(
                        ["state"],
                        state_type_data,
                        default_attrs.get("state_types", None),
                    ),
                    state_names=(
                        ["state"],
                        state_name_data,
                        default_attrs.get("state_names", None),
                    ),
                )

                mark_variable_assigned(dataset["state_types"])
                mark_variable_assigned(dataset["state_names"])

        num_states = dataset.sizes["state"]
        default_states = list(range(1, num_states + 1))

        if "state" not in dataset.coords:
            dataset = dataset.assign_coords(
                {"state": ("state", default_states, default_attrs["state"])}
            )

        # Fix statecomb if missing:
        if "statecomb" in dataset.dims:
            if "from" not in dataset.coords or "to" not in dataset.coords:
                statecomb_coords = xr.Coordinates.from_pandas_multiindex(
                    pd.MultiIndex.from_tuples(
                        combinations(default_states, 2), names=["from", "to"]
                    ),
                    dim="statecomb",
                )
                dataset = dataset.assign_coords(statecomb_coords)
            dataset["statecomb"].attrs.update(default_attrs.get("statecomb", {}))
            mark_variable_assigned(dataset["statecomb"])
            dataset["from"].attrs.update(default_attrs.get("from", {}))
            mark_variable_assigned(dataset["from"])
            dataset["to"].attrs.update(default_attrs.get("to", {}))
            mark_variable_assigned(dataset["to"])

        if "full_statecomb" in dataset.dims:
            if (
                "full_statecomb_from" not in dataset.coords
                or "full_statecomb_to" not in dataset.coords
            ):
                full_statecombs_coords = xr.Coordinates.from_pandas_multiindex(
                    pd.MultiIndex.from_tuples(
                        permutations(default_states, 2),
                        names=["full_statecomb_from", "full_statecomb_to"],
                    ),
                    dim="full_statecomb",
                )
                dataset = dataset.assign_coords(full_statecombs_coords)
            dataset["full_statecomb"].attrs.update(
                default_attrs.get("full_statecomb", {})
            )
            mark_variable_assigned(dataset["full_statecomb"])
            dataset["full_statecomb_from"].attrs.update(
                default_attrs.get("full_statecomb_from", {})
            )
            mark_variable_assigned(dataset["full_statecomb_from"])
            dataset["full_statecomb_to"].attrs.update(
                default_attrs.get("full_statecomb_to", {})
            )
            mark_variable_assigned(dataset["full_statecomb_to"])

    # Restore `direction` coordinate if missing.
    if "direction" in dataset.dims and "direction" not in dataset.coords:
        dataset = dataset.assign_coords(
            direction=(
                "direction",
                ["x", "y", "z"],
                default_attrs.get("direction", None),
            )
        )

    # Set trajectory-level attributes
    if "misc_attrs" in shnitsel_meta:
        dataset.attrs.update(shnitsel_meta["misc_attrs"])

    # Perform a check of the dimension sizes specified in the metadata if present
    if "dims" in shnitsel_meta:
        for dimname, dimdict in shnitsel_meta["dims"].items():
            if dimname == "tmp" or dimname == "state_or_statecomb":
                # Skip artificial dimensions
                continue
            if dimname == "time" or dimname == "frame":
                # Skip leading dimensions that may have been split
                continue
            dim_length = dimdict["length"] if "length" in dimdict else -1
            if dim_length >= 0:
                if dimname not in dataset.dims or dim_length != dataset.sizes[dimname]:
                    msg = "Size of dimension %(dimname)s in dataset parsed from ASE database has length inconsistent with metadata of ASE file. Was %(ds_dim_size)d but metadata specifies %(dim_spec)d"
                    params = {
                        'dimname': dimname,
                        'ds_dim_size': dataset.sizes.get(dimname, -1),
                        'dim_spec': dim_length,
                    }
                    logging.info(msg, params)
                    # raise ValueError(msg % params)

    if "est_level" not in dataset.attrs:
        if 'ReferenceMethod' in db_meta:
            parts = str(db_meta['ReferenceMethod']).strip().split("/")
            if len(parts) >= 2:
                dataset.attrs["est_level"] = parts[0]
                if "theory_basis_set" not in dataset.attrs:
                    dataset.attrs["theory_basis_set"] = parts[1]
            else:
                dataset.attrs["est_level"] = str(db_meta['ReferenceMethod'])
                if "theory_basis_set" not in dataset.attrs:
                    dataset.attrs["theory_basis_set"] = "unknown"

    return dataset


# TODO: FIXME: Check the return type and tree interaction
def read_ase(
    db_path: pathlib.Path,
    db_format: Literal['spainn', 'schnet'] | None = None,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset | TreeNode[Any, Trajectory | Frames | ShnitselDataset | xr.Dataset]:
    """Reads an ASE DB containing data in the SPaiNN or SchNet format

    Parameters
    ----------
    db_path: pathlib.Path
        Path to the database
    db_format: Literal['spainn', 'schnet'] | None, optional
        Must be one of 'spainn' or 'schnet' or None; determines interpretation of array shapes If None is provided, no shape will be assumed
    loading_parameters: LoadingParameters
        Potentially configured parameters to overwrite loading behavior

    Returns
    -------
    xr.Dataset
        An `xr.Dataset` of frames. Potentially with a `time` coordinate.

    Raises
    ------
    ValueError
        If `db_format` is not one of 'spainn' or 'schnet'
    FileNotFoundError
        If `db_path` is not a file
    ValueError
        If `db_path` does not contain data corresponding to the format `db_format`
    """
    if db_format is not None and db_format not in ['spainn', 'schnet']:
        raise ValueError("'db_format' should be one of 'schnet' or 'spainn'")

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Could not find databse at {db_path}")

    ase_default_attrs = get_default_input_attributes("ase", loading_parameters)

    with connect(db_path) as db:
        data_vars = {}
        coord_vars = {}
        found_rows = 0
        # available_varnames = next(db.select()).data.keys()
        # print(available_varnames)

        tmp_data_in = {
            "atXYZ": [],
        }

        set_list: list[
            tuple[
                dict[str, list[np.ndarray | float | int | str]],
                dict[str, tuple[int, ...]],
            ]
        ] = []

        curr_set_vars: dict[str, list[np.ndarray | float | int | str]] | None = None
        curr_set_shapes: dict[str, tuple[int, ...]] | None = None

        contract_keys_coords = [
            'delta_t',
            't_max',
            'max_ts',
            "charge",
        ]
        contract_keys_attrs = [
            "input_format",
            "input_type",
            "input_format_version",
            "completed",
        ]
        keys_change_is_mismatch = (
            [
                'atNames',
                'atrajectory',
            ]
            + contract_keys_coords
            + contract_keys_attrs
        )

        for row in db.select():
            row_vars: dict[str, np.ndarray | float | int | str] = {}
            row_shapes: dict[str, tuple[int, ...]] = {}
            for key, value in row.data.items():
                var_tmp = np.array(value)
                row_vars[key] = var_tmp
                row_shapes[key] = var_tmp.shape

            kv_pairs = row.key_value_pairs
            row_atoms = row.toatoms()
            # TODO: FIXME: deal with different atoms/compounds in the same DB.
            if 'atNames' not in row_vars:
                # tmp_data_in['atNames'] = []
                row_vars['atNames'] = row_atoms.get_chemical_symbols()
                row_shapes['atNames'] = (len(row_vars['atNames']),)
            else:
                new_symbols = row_atoms.get_chemical_symbols()
                if row_vars['atNames'] != new_symbols:
                    raise ValueError(
                        f"Mismatch between symbols of different rows. Previously read: {row_vars['atNames']} now {new_symbols}."
                        "We currently do not support reading multiple different compounds from one ASE db."
                    )

            # if row_atoms.has("positions"):
            row_vars['atXYZ'] = row_atoms.get_positions()
            row_shapes['atXYZ'] = (len(row_vars['atXYZ']),)

            if row_atoms.has("momenta"):
                row_vars['velocities'] = row_atoms.get_velocities()
                row_shapes['velocities'] = (len(row_vars['velocities']),)

            # if "time" in row_atoms.info:
            #     row_vars["time"] = float(row_atoms.info["time"])
            #     row_shapes['time'] = (1,)

            if "time" in kv_pairs:
                row_vars["time"] = float(kv_pairs["time"])
                row_shapes['time'] = (1,)

            # if "trajid" in row_atoms.info:
            #     row_vars["atrajectory"] = int(row_atoms.info["trajid"])
            #     row_shapes['atrajectory'] = (1,)
            if "trajid" in kv_pairs:
                row_vars["trajid"] = float(kv_pairs["trajid"])
                row_shapes['trajid'] = (1,)
            # if "trajectory" in row_atoms.info:
            #     row_vars["atrajectory"] = int(row_atoms.info["trajectory"])
            #     row_shapes['atrajectory'] = (1,)
            if "trajectory" in kv_pairs:
                row_vars["atrajectory"] = float(kv_pairs["trajectory"])
                row_shapes['atrajectory'] = (1,)
            # if "atrajectory" in row_atoms.info:
            #     row_vars["atrajectory"] = int(row_atoms.info["atrajectory"])
            #     row_shapes['atrajectory'] = (1,)
            if "atrajectory" in kv_pairs:
                row_vars["atrajectory"] = float(kv_pairs["atrajectory"])
                row_shapes['atrajectory'] = (1,)

            # if "delta_t" in row_atoms.info:
            #     row_vars["delta_t"] = float(row_atoms.info["delta_t"])
            #     row_shapes['delta_t'] = (1,)
            if "delta_t" in kv_pairs:
                row_vars["delta_t"] = float(kv_pairs["delta_t"])
                row_shapes['delta_t'] = (1,)

            # if "t_max" in row_atoms.info:
            #     row_vars["t_max"] = float(row_atoms.info["t_max"])
            #     row_shapes['t_max'] = (1,)
            if "t_max" in kv_pairs:
                row_vars["t_max"] = float(kv_pairs["t_max"])
                row_shapes['t_max'] = (1,)

            # if "max_ts" in row_atoms.info:
            #     row_vars["max_ts"] = int(row_atoms.info["max_ts"])
            #     row_shapes['max_ts'] = (1,)
            if "max_ts" in kv_pairs:
                row_vars["max_ts"] = float(kv_pairs["max_ts"])
                row_shapes['max_ts'] = (1,)

            if row.charge != 0.0:
                row_vars["charge"] = float(row.charge)
                row_shapes['charge'] = (1,)
            if "charge" in kv_pairs:
                row_vars["charge"] = float(kv_pairs["charge"])
                row_shapes['charge'] = (1,)

            if "input_format" in kv_pairs:
                row_vars["input_format"] = str(kv_pairs["input_format"])
                row_shapes['input_format'] = (1,)
            if "input_type" in kv_pairs:
                row_vars["input_type"] = str(kv_pairs["input_type"])
                row_shapes['input_type'] = (1,)
            if "input_format_version" in kv_pairs:
                vals = str(kv_pairs["input_format_version"])
                guard_prefix = "v__"
                if vals.startswith(guard_prefix):
                    vals = vals[len(guard_prefix) :]
                row_vars["input_format_version"] = vals
                row_shapes['input_format_version'] = (1,)

            for sup_coord in SUPPORTED_COORD_KEYS:
                if sup_coord in kv_pairs:
                    if (
                        sup_coord
                        not in {
                            'trajid',
                            'trajid_',
                            'trajectory',
                            'atrajectory',
                            'time',
                        }
                        and sup_coord not in row_vars
                    ):
                        row_vars[sup_coord] = kv_pairs[sup_coord]
                        row_shapes[sup_coord] = (1,)

            found_rows += 1

            is_mismatched = False
            if curr_set_shapes is None or curr_set_vars is None:
                is_mismatched = True
            else:
                if curr_set_shapes != row_shapes:
                    is_mismatched = True
                elif "time" in row_vars and (
                    "time" not in curr_set_vars
                    or float(row_vars["time"]) <= float(curr_set_vars["time"][-1])
                ):
                    # Time reset, should indicate a trajectory change
                    is_mismatched = True
                else:
                    for key in keys_change_is_mismatch:
                        if key in row_vars and (
                            key not in curr_set_vars
                            or row_vars[key] != curr_set_vars[key][-1]
                        ):
                            logging.debug(
                                f" Mismatch in var {key=}: {row_vars[key]=} != {curr_set_vars[key][-1]=}"
                            )
                            # E.g. Atoms have changed. Indicates a completely different compound.
                            # Trajectory changed
                            is_mismatched = True
                            break
                        # else:
                        # print(
                        #     f"{key=}: {row_vars.get(key,'~')=} != {curr_set_vars.get(key,'~')=}"
                        # )

            if is_mismatched:
                if curr_set_shapes is not None and curr_set_vars is not None:
                    set_list.append((curr_set_vars, curr_set_shapes))

                curr_set_shapes = row_shapes
                curr_set_vars = {}

            if TYPE_CHECKING:
                assert curr_set_shapes is not None and curr_set_vars is not None

            for key, value in row_vars.items():
                if key not in curr_set_vars:
                    curr_set_vars[key] = []
                curr_set_vars[key].append(value)

    # If there are no valid rows, raise a ValueError
    if found_rows == 0:
        raise ValueError(
            f"No rows with the appropriate format for `{db_format=}` were found in {db_path}"
        )
    # Add the final dataset to the collection
    if curr_set_shapes is not None and curr_set_vars is not None:
        set_list.append((curr_set_vars, curr_set_shapes))

    metadata = db.metadata
    (
        global_var_shapes,
        global_coord_shapes,
        global_leading_dimension_name,
        global_db_format,
    ) = shapes_from_metadata(metadata, db_format)

    # First the compound string as key, then the dict of shapes as key to differentiate different setups
    res_data: dict[str, dict[tuple, list]] = {}

    for tmp_data_in, tmp_shapes_in in set_list:
        stacked_tmp_data = {k: np.stack(v) for k, v in tmp_data_in.items()}
        leading_dimension_rename_target = None

        if 'time' not in tmp_data_in:
            # No time means we can only index with `frame` coordinate
            leading_dimension_rename_target = "frame"
        else:
            if "atrajectory" in tmp_data_in:
                leading_dimension_rename_target = "frame"
            else:
                # We have time and no active trajectory data.
                # Must be time series
                leading_dimension_rename_target = "time"
        (
            local_var_shapes,
            local_coord_shapes,
            db_format,
            leading_dimension_rename_target,
        ) = complete_shapes_guesses_from_variables(
            stacked_tmp_data,
            global_var_shapes,
            global_coord_shapes,
            global_db_format,
            global_leading_dimension_name,
            leading_dimension_rename_target,
        )

        ds_attrs = {}

        int_data_keys = [
            "state_types",
            "astate",
            "sdiag",
            "trajid",
            "atrajectory",
            "trajectory",
            "trajid_",
            "state",
            "atom",
        ]

        for k, data_array in stacked_tmp_data.items():
            if k in contract_keys_coords:
                # Deal with coordinates that should be 0-dimensional
                coord_vars[k] = (
                    [],
                    data_array[0],
                    ase_default_attrs.get(k, None),
                )
            elif k in contract_keys_attrs:
                # Deal with attributes that should be 0-dimensional
                ds_attrs[k] = data_array[0]
            elif k in local_var_shapes:
                # if str(k) == "socs":
                #     raise ValueError(
                #         f"Read variable {k} with shape: {shapes[k]} and numpy shape: {data_array.shape}"
                #     )
                data_vars[k] = (
                    local_var_shapes[k],
                    data_array,
                    (ase_default_attrs[k] if k in ase_default_attrs else None),
                )
            elif k in local_coord_shapes:
                if k == "atNames":
                    # Deal with coordinates that should be 1-dimensional
                    coord_vars[k] = (
                        local_coord_shapes["atNames"],
                        data_array[0],
                        ase_default_attrs.get(k, None),
                    )

                    coord_vars["atNums"] = (
                        local_coord_shapes["atNums"],
                        np.array(
                            [get_atom_number_from_symbol(x) for x in data_array[0]]
                        ),
                        ase_default_attrs.get("atNums", None),
                    )
                elif k == 'atrajectory':
                    # Deal with trajectory index being a variable
                    if leading_dimension_rename_target == "frame":
                        coord_vars[k] = (
                            ['frame'],
                            data_array,
                            ase_default_attrs.get(k, None),
                        )
                    else:
                        coord_vars["trajectory"] = xr.DataArray(
                            dims=[],
                            data=data_array[0],
                            attrs=ase_default_attrs.get(k, None),
                        ).astype(int)
                else:
                    coord_vars[k] = (
                        local_coord_shapes[k],
                        data_array,
                        ase_default_attrs.get(k, None),
                    )
            else:
                logging.warning(
                    "Dropping data entry %(key)s due to missing shape information",
                    {'key': k},
                )
            if k in int_data_keys:
                if k in coord_vars:
                    ldims, ldata, lattrs = coord_vars[k]
                    # print(f"{k=} -> {coord_vars[k]=}")
                    coord_vars[k] = xr.DataArray(
                        dims=ldims,
                        data=ldata,
                        attrs=lattrs,
                    ).astype(int)

            # atXYZ = np.stack([row.positions for row in db.select()])
            # data_vars['atXYZ'] = ['frame', 'atom', 'direction'], atXYZ
            # atNames = ['atom'], next(db.select()).symbols
        nstates: int = -1
        if "dims" in metadata:
            if "state" in metadata["dims"]:
                nstates = metadata["dims"]["state"]["length"]

        if "states" in metadata:
            nstates = len(str(metadata["states"]).split())

        if nstates < 0:
            logging.debug("Extracting number of states from shape of energy array")
            if "energy" in data_vars:
                nstates = data_vars['energy'][0].shape[-1]

        if 'dipoles' in data_vars and nstates > 0:
            dipoles = data_vars['dipoles'][1]
            dip_perm = dipoles[:, :nstates, :]
            dip_trans = dipoles[:, nstates:, :]
            del data_vars['dipoles']

            data_vars['dip_perm'] = (
                [global_leading_dimension_name, 'state', 'direction'],
                dip_perm,
                ase_default_attrs["dip_perm"],
            )
            data_vars['dip_trans'] = (
                [global_leading_dimension_name, 'statecomb', 'direction'],
                dip_trans,
                ase_default_attrs["dip_perm"],
            )

        # print(data_vars["atXYZ"][1].shape)
        # print(data_vars["atXYZ"][1].shape)
        # print(coord_vars["frame"])
        frames = xr.Dataset(data_vars, coord_vars, attrs=ds_attrs)

        # Set flags to mark as assigned
        for k in coord_vars.keys():
            mark_variable_assigned(frames[k])
        for k in data_vars.keys():
            mark_variable_assigned(frames[k])

        if db_format == 'spainn':
            # Only squeeze if the tmp dimension is there
            if 'tmp' in frames.dims:
                frames = frames.squeeze('tmp')
            else:
                logging.warning(
                    "Input of type `spainn` did not yield a `tmp` dimension, indicating missing energy. Input file %(path)s may be malformed.",
                    {'path': db_path},
                )

        # Deal with us not identifying the leading dimension from metadata alone.
        if global_leading_dimension_name == dummy_leading_dim:
            # Only rename if a variable with the dimension was created. Otherwise an error would trigger in rename
            if global_leading_dimension_name in frames.dims:
                if leading_dimension_rename_target is None:
                    if (
                        "time" in frames.coords
                        and "trajid" not in frames.coords
                        and "trajid_" not in frames.coords
                        and "atrajectory" not in frames.coords
                        and "trajectory" not in frames.coords
                    ):
                        leading_dimension_rename_target = "time"
                    else:
                        leading_dimension_rename_target = "frame"

                frames = frames.swap_dims(
                    {global_leading_dimension_name: leading_dimension_rename_target}
                )

        # Restore missing coordinates, metadata, units, etc.
        frames = apply_dataset_meta_from_db_metadata(
            frames, metadata, ase_default_attrs
        )

        # Order dimensions in default shnitsel order
        shnitsel_default_order = [
            "trajectory",
            "frame",
            "time",
            ...,
            "state",
            "statecomb",
            "full_statecomb",
            "atom",
            "direction",
        ]

        # Try and guess charge or check if the charge is appropriate

        initial_charge = frames.attrs.get(
            "charge", frames.charge if "charge" in frames.coords else None
        )
        needs_charge_guess = initial_charge is None
        initial_positions = frames.atXYZ.isel(frame=0) * 2

        # print(f"Initial charge: {initial_charge}")
        # print(f"Initial positions: {initial_positions=}")

        if initial_charge is not None:
            logging.info(
                "Checking initial charge of molecule as read from ASE database: `%s`",
                str(initial_charge),
            )
            try:
                to_mol(initial_positions, charge=int(round(initial_charge)))
            except Exception as e:
                # print(e)
                logging.info(
                    "Initial charge guess `%s` did not yield correct molecule structure. Will be recalculated recalculation.",
                    str(initial_charge),
                )
                needs_charge_guess = True

        if needs_charge_guess:
            # print(f"Making a guess after: {initial_charge=}")
            if "charge" in frames.attrs:
                del frames.attrs["charge"]
            if "charge" not in frames.coords:
                frames = frames.assign_coords(
                    charge=([], 0.0, ase_default_attrs.get("charge", None))
                )

            attempt_order = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5]
            for charge_guess in attempt_order:
                try:
                    to_mol(initial_positions, charge=int(round(charge_guess)))
                    frames.coords["charge"].values = charge_guess
                    needs_charge_guess = False
                    break
                except Exception as e:
                    continue

            if needs_charge_guess:
                logging.warning(
                    "Could not derive charge for molecule read from ASE database initial guess was %s. Removing charge property for now. Please assign manually with `set_charge(charge)`",
                    str(initial_charge),
                )
                frames.drop_vars("charge", errors="ignore")
            else:
                logging.info(
                    "Successfully derived molecular charge `%s` from information in ASE database.",
                    str(frames.coords["charge"]),
                )

        frames = wrap_dataset(
            normalize_dataset(
                frames.transpose(*shnitsel_default_order, missing_dims="ignore")
            ),
            Trajectory | Frames | ShnitselDataset,
        )
        compound_string = ""
        dim_dict: dict = {
            k: v for k, v in frames.sizes.items() if k not in {"time", "frame"}
        }
        if "atNames" in frames.coords:
            compound_string = "-".join(frames.coords["atNames"].values)
            dim_dict["atNames"] = compound_string

        if isinstance(frames, Frames):
            if frames.has_coordinate("time") and not frames.is_multi_trajectory:
                frames = frames.as_trajectory

        if compound_string not in res_data:
            res_data[compound_string] = {}

        dim_dict_key = tuple(sorted(dim_dict.items()))

        if dim_dict_key not in res_data[compound_string]:
            res_data[compound_string][dim_dict_key] = []

        res_data[compound_string][dim_dict_key].append(frames)

    compounds = {}
    for compound_string, compound_groups in res_data.items():
        compound_subgroups = {}
        for i_g, group_list in enumerate(compound_groups.values()):
            leaves = {
                str(i): DataLeaf(name=str(i), data=f) for i, f in enumerate(group_list)
            }
            group = DataGroup(name=str(i_g), children=leaves)
            compound_subgroups[str(i_g)] = group

        compounds[compound_string] = CompoundGroup(
            compound_info=CompoundInfo(compound_name=compound_string),
            children=compound_subgroups,
        )

    return ShnitselDBRoot(compounds=compounds)
