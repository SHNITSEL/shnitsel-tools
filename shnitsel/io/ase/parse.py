from itertools import combinations, permutations
import json
import logging
import os
import pathlib
from typing import Any, Literal

from ase.db import connect
import numpy as np
import pandas as pd
import xarray as xr

from shnitsel.data.trajectory_format import Trajectory
from shnitsel.io.shared.helpers import LoadingParameters
from shnitsel.io.shared.trajectory_setup import (
    RequiredTrajectorySettings,
    assign_required_settings,
)
from shnitsel.io.shared.variable_flagging import mark_variable_assigned
from shnitsel.units.defaults import get_default_input_attributes
from shnitsel.data.atom_helpers import get_atom_number_from_symbol

dummy_leading_dim: str = "leading_dim_unknown"
multi_level_prefix: str = "_MultiIndex_levels_for_"


def shapes_from_metadata(
    db_meta: dict, db_format: Literal['spainn', 'schnet'] | None = None
) -> tuple[dict[str, list[str]], dict[str, list[str]], str]:
    """Function to assign shapes based on the chosen db_format and potential information in the metadata of a database.

    If conflicting information on the db_format is provided and present in the database, en error will be raised.

    Args:
        db_meta (dict): The metadata dict of an ASE database.
        db_format (Literal['spainn', 'schnet'] | None, optional): The requested format of the database. Defaults to None.

    Return:
        dict[str, list[str]]: Dict of data_var shapes
        dict[str, list[str]]: Dict of coordinate shapes
        str: The name of the leading dimension. Should be `frame` or `time`, but can be `leading_dim_unknown` if unknown

    Raises:
        ValueError: If a db_format of database was requested that conflicts with the format of the database.
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
        "trajid": [leading_dim_name],  # This only exists for frame dimensions
        "from": ["statecomb"],
        "to": ["statecomb"],
        "full_statecomb_from": ["full_statecomb"],
        "full_statecomb_to": ["full_statecomb"],
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
            logging.info(f"Automatically detected format: {db_format}")

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

    return shapes, coord_shapes, leading_dim_name


def _json_deserialize_ndarray(value: str) -> Any:
    if isinstance(value, str):
        value_d = json.loads(value)
    else:
        value_d = value

    try:
        # if isinstance(value_d, dict):
        #     print("Is a dict")
        # if "__ndarray" in value_d:
        config = value_d["__ndarray"]

        entries = config["entries"]
        dtype_descr = np.dtype([tuple(i) for i in config["dtype"]])

        return np.array(entries, dtype=dtype_descr)
    except TypeError as e:
        pass
    except KeyError as e:
        pass

    return value


def apply_dataset_meta_from_db_metadata(
    dataset: Trajectory,
    db_meta: dict,
    default_attrs: dict,
) -> Trajectory:
    """Apply attributes from db metadata and perform some validation checks on the result.

    Loads remaining missing coordinate variables from db metadata if available.
    Checks size of resulting dimensions if specified in db metadata.
    Further initializes the multi indices if specified in the metadata.

    Args:
        dataset (Trajectory): Trajectory dataset parsed from ASE db
        db_meta (dict): Metadata from the trajectory db file
        default_attrs (dict): Attributes to apply to variables by default


    Returns:
        Trajectory: Dataset with attributes set from from db metadata and dimension sizes asserted
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
        # Add dummy frame coordinate values treating all entries as initial conditions/static data in different trajectories
        if "frame" not in dataset:
            frame_vals = np.arange(0, dataset.sizes["frame"], 1)
            dataset = dataset.assign_coords(
                {
                    "frame": (["frame"], frame_vals, default_attrs["frame"]),
                    "trajid": (["frame"], frame_vals, default_attrs["trajid"]),
                    "time": (["frame"], frame_vals * 0.0, default_attrs["time"]),
                }
            )
    elif "time" in dataset.dims:
        # Fill in missing time coordinate with dummy values if no frame is set as dimension
        if "time" not in dataset:
            time_vals = np.arange(0, dataset.sizes["time"], 1) * (
                dataset.attrs["delta_t"] if "delta_t" in dataset.attrs else 1.0
            )
            dataset = dataset.assign_coords(
                {
                    "time": (["time"], time_vals, default_attrs["time"]),
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
        if "atXYZ" in dataset and "unit" not in dataset["atXYZ"].attrs:
            dataset["atXYZ"].attrs["unit"] = db_meta["_distance_unit"]

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

    delta_t = dataset.attrs["delta_t"] if "delta_t" in dataset.attrs else None
    if delta_t is None:
        # Try and extract from time info
        if "time" in dataset:
            # print(dataset["time"])
            # print(dataset["trajid"])
            diff_t = list(set(dataset["time"].values))
            diff_t = sorted(diff_t)
            if len(diff_t) > 1:
                delta_t = diff_t[1] - diff_t[0]
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
        t_max=dataset.attrs["t_max"] if "t_max" in dataset.attrs else -1,
        delta_t=delta_t,
        max_ts=dataset.attrs["max_ts"]
        if "max_ts" in dataset.attrs
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

    assign_required_settings(dataset, extract_settings)

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
                        default_attrs["state_types"],
                    ),
                    state_names=(
                        ["state"],
                        state_name_data,
                        default_attrs["state_names"],
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
            dataset["statecomb"].attrs.update(default_attrs["statecomb"])
            mark_variable_assigned(dataset["statecomb"])
            dataset["from"].attrs.update(default_attrs["from"])
            mark_variable_assigned(dataset["from"])
            dataset["to"].attrs.update(default_attrs["to"])
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
            dataset["full_statecomb"].attrs.update(default_attrs["full_statecomb"])
            mark_variable_assigned(dataset["full_statecomb"])
            dataset["full_statecomb_from"].attrs.update(
                default_attrs["full_statecomb_from"]
            )
            mark_variable_assigned(dataset["full_statecomb_from"])
            dataset["full_statecomb_to"].attrs.update(
                default_attrs["full_statecomb_to"]
            )
            mark_variable_assigned(dataset["full_statecomb_to"])

    # Set trajectory-level attributes
    if "misc_attrs" in shnitsel_meta:
        dataset.attrs.update(shnitsel_meta["misc_attrs"])

    # Perform a check of the dimension sizes specified in the metadata if present
    if "dims" in shnitsel_meta:
        for dimname, dimdict in shnitsel_meta["dims"].items():
            if dimname == "tmp" or dimname == "state_or_statecomb":
                # Skip artificial dimensions
                continue
            dim_length = dimdict["length"] if "length" in dimdict else -1
            if dim_length >= 0:
                if dim_length != dataset.sizes[dimname]:
                    msg = f"Size of dimension {dimname} in dataset parsed from ASE database has length inconsistent with metadata of ASE file. Was {dataset.sizes[dimname]} but metadata specifies {dim_length}"
                    logging.error(msg)
                    raise ValueError(msg)

    if "est_level" not in dataset.attrs:
        if 'ReferenceMethod' in db_meta:
            # TODO: FIXME: Possibly split up into basis and method?
            dataset.attrs["est_level"] = db_meta['ReferenceMethod']

    return dataset


def read_ase(
    db_path: pathlib.Path,
    db_format: Literal['spainn', 'schnet'] | None = None,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset:
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
        An `xr.Dataset` of frames

    Raises
    ------
    ValueError
        If `db_format` is not one of 'spainn' or 'schnet'
    FileNotFoundError
        If `db_path` is not a file
    ValueError
        If `db_path` does not contain data corresponding to the format `db_format`
    """

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Could not find databse at {db_path}")

    ase_default_attrs = get_default_input_attributes("ase", loading_parameters)

    with connect(db_path) as db:
        metadata = db.metadata
        shapes, coord_shapes, leading_dimension_name = shapes_from_metadata(
            metadata, db_format
        )
        leading_dimension_rename_target = None

        data_vars = {}
        coord_vars = {}
        found_rows = 0
        # available_varnames = next(db.select()).data.keys()
        # print(available_varnames)

        tmp_data_in = {
            "atXYZ": [],
        }

        for row in db.select():
            for key, value in row.data.items():
                if key not in tmp_data_in:
                    tmp_data_in[key] = []

                tmp_data_in[key].append(value)
            row_atoms = row.toatoms()
            # TODO: FIXME: deal with different atoms/compounds in the same DB.
            if 'atNames' not in tmp_data_in:
                # tmp_data_in['atNames'] = []
                tmp_data_in['atNames'] = row_atoms.get_chemical_symbols()
            else:
                new_symbols = row_atoms.get_chemical_symbols()
                if tmp_data_in['atNames'] != new_symbols:
                    raise ValueError(
                        f"Mismatch between symbols of different rows. Previously read: {tmp_data_in['atNames']} now {new_symbols}. We currently do not support reading multiple different compounds from one ASE db."
                    )

            # if row_atoms.has("positions"):
            tmp_data_in['atXYZ'].append(row_atoms.get_positions())

            if row_atoms.has("momenta"):
                if 'velocities' not in tmp_data_in:
                    tmp_data_in['velocities'] = []
                tmp_data_in['velocities'].append(row_atoms.get_velocities())

            if "time" in row_atoms.info:
                if "time" not in tmp_data_in:
                    tmp_data_in["time"] = []
                tmp_data_in["time"].append(float(row_atoms.info["time"]))

            if "trajid" in row_atoms.info:
                if "trajid" not in tmp_data_in:
                    tmp_data_in["trajid"] = []
                    leading_dimension_rename_target = "frame"

                tmp_data_in["trajid"].append(row_atoms.info["trajid"])

            found_rows += 1

    # If there are no valid rows, raise a ValueError
    if found_rows == 0:
        raise ValueError(
            f"No rows with the appropriate format for `{db_format=}` were found in {db_path}"
        )

    for k, v in tmp_data_in.items():
        data_array = np.stack(v)
        if k in shapes:
            # if str(k) == "socs":
            #     raise ValueError(
            #         f"Read variable {k} with shape: {shapes[k]} and numpy shape: {data_array.shape}"
            #     )
            data_vars[k] = (
                shapes[k],
                data_array,
                (ase_default_attrs[k] if k in ase_default_attrs else None),
            )
        elif k in coord_shapes:
            if k == "atNames":
                coord_vars[k] = (
                    coord_shapes[k],
                    data_array,
                    (ase_default_attrs[k] if k in ase_default_attrs else None),
                )

                coord_vars["atNums"] = (
                    coord_shapes["atNums"],
                    np.array([get_atom_number_from_symbol(x) for x in data_array]),
                    (
                        ase_default_attrs["atNums"]
                        if "atNums" in ase_default_attrs
                        else None
                    ),
                )
            else:
                coord_vars[k] = (
                    coord_shapes[k],
                    data_array,
                    (ase_default_attrs[k] if k in ase_default_attrs else None),
                )
        else:
            logging.warning(f"Dropping data entry {k} due to missing shape information")

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
            nstates = data_vars['energy'][1].shape[1]

    if 'dipoles' in data_vars and nstates > 0:
        dipoles = data_vars['dipoles'][1]
        dip_perm = dipoles[:, :nstates, :]
        dip_trans = dipoles[:, nstates:, :]
        del data_vars['dipoles']

        data_vars['dip_perm'] = (
            [leading_dimension_name, 'state', 'direction'],
            dip_perm,
            ase_default_attrs["dip_perm"],
        )
        data_vars['dip_trans'] = (
            [leading_dimension_name, 'statecomb', 'direction'],
            dip_trans,
            ase_default_attrs["dip_perm"],
        )

    # print(data_vars["atXYZ"][1].shape)
    # print(data_vars["atXYZ"][1].shape)
    # print(coord_vars["frame"])
    frames = xr.Dataset(data_vars).assign_coords(coord_vars)

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
                f"Input of type `spainn` did not yield a `tmp` dimension, indicating missing energy. Input file {db_path} may be malformed."
            )

    # Deal with us not identifying the leading dimension from metadata alone.
    if leading_dimension_name == dummy_leading_dim:
        # Only rename if a variable with the dimension was created. Otherwise an error would trigger in rename
        if leading_dimension_name in frames.dims:
            if leading_dimension_rename_target is None:
                if "time" in coord_vars and "trajid" not in coord_vars:
                    leading_dimension_rename_target = "time"
                else:
                    leading_dimension_rename_target = "frame"

            frames = frames.rename(
                {leading_dimension_name: leading_dimension_rename_target}
            )

    frames = apply_dataset_meta_from_db_metadata(frames, metadata, ase_default_attrs)

    # Order dimensions in default shnitsel order
    shnitsel_default_order = [
        "frame",
        "trajid",
        "time",
        "state",
        "statecomb",
        "full_statecomb",
        "atom",
        "direction",
    ]
    frames = frames.transpose(*shnitsel_default_order, missing_dims="ignore")
    return frames
