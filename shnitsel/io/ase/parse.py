import logging
import os
import pathlib
from typing import Literal

from ase.db import connect
import numpy as np
import xarray as xr

from shnitsel.data.trajectory_format import Trajectory
from shnitsel.io.helpers import LoadingParameters
from shnitsel.units.defaults import get_default_input_attributes
from shnitsel.io.helpers import get_atom_number_from_symbol

dummy_leading_dim: str = "leading_dim_unknown"
multi_level_prefix: str = "_MultiIndex_levels_for_"


def shapes_from_metadata(
    db_meta: dict, kind: Literal['spainn', 'schnet'] | None = None
) -> tuple[dict[str, list[str]], dict[str, list[str]], str]:
    """Function to assign shapes based on the chosen kind and potential information in the metadata of a database.

    If conflicting information on the format/kind is provided and present in the database, en error will be raised.

    Args:
        db_meta (dict): The metadata dict of an ASE database.
        kind (Literal['spainn', 'schnet'] | None, optional): The requested format of the database. Defaults to None.

    Return:
        dict[str, list[str]]: Dict of data_var shapes
        dict[str, list[str]]: Dict of coordinate shapes
        str: The name of the leading dimension. Should be `frame` or `time`, but can be `leading_dim_unknown` if unknown

    Raises:
        ValueError: If a kind of database was requested that conflicts with the format of the database.
    """

    if 'shnitsel_leading_dim' in db_meta:
        leading_dim_name = db_meta['shnitsel_leading_dim']
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
        'socs': [leading_dim_name, 'statecomb'],
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
        'socs': [leading_dim_name, 'statecomb'],
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
        "statecomb": ["statecomb"],
    }

    if "db_format" in db_meta:
        meta_format = db_meta["db_format"]
        if meta_format not in ["schnet", "spainn"]:
            raise ValueError(
                f"Database is of unsupported format: {meta_format}. Only `schnet` and `spainn` are supported."
            )

        if kind is None:
            kind = meta_format
            logging.info(f"Automatically detected format: {kind}")

        if meta_format != kind:
            raise ValueError(
                f"Database is of format: {meta_format} instead of requested format {kind}."
            )
    shapes: dict[str, list[str]]
    # Determine basis shapes based on the format
    if kind == 'schnet':
        shapes = schnet_shapes
    elif kind == 'spainn':
        shapes = spainn_shapes
    elif kind is None:
        shapes = {}
        logging.warning(
            "Correct format could not be extracted from the database metadata. No dimension names assigned"
        )
    else:
        raise ValueError(f"'kind' should be one of 'schnet' or 'spainn', not '{kind}'.")

    # Read further shape data from the database
    if "var_meta" in db_meta:
        variable_metadata = db_meta["var_meta"]
        for varname, vardict in variable_metadata.items():
            if "dims" in vardict:
                shapes[varname] = vardict["dims"]

    if "coords" in db_meta:
        coord_metadata = db_meta["coords"]
        for coordname, coorddict in coord_metadata.items():
            if "dims" in coorddict:
                coord_shapes[coordname] = coorddict["dims"]

    return shapes, coord_shapes, leading_dim_name


def apply_dataset_meta_from_db_metadata(
    dataset: Trajectory, db_meta: dict
) -> Trajectory:
    """Apply attributes from db metadata and perform some validation checks on the result.

    Loads remaining missing coordinate variables from db metadata if available.
    Checks size of resulting dimensions if specified in db metadata.
    Further initializes the multi indices if specified in the metadata.

    Args:
        dataset (Trajectory): Trajectory dataset parsed from ASE db
        db_meta (dict): Metadata from the trajectory db file

    Returns:
        Trajectory: Dataset with attributes set from from db metadata and dimension sizes asserted
    """
    # Restore missing coordinates
    if "coords" in db_meta:
        coords_data = db_meta["coords"]
        for coordname, coorddict in coords_data.items():
            if coordname not in dataset.coords:
                dataset = dataset.assign_coords(
                    {coordname: (coorddict["dims"], coorddict["values"])}
                )

    # Potentially reconstruct multiindex levels
    if (
        "_MultiIndex_levels_from_attrs" in db_meta
        and db_meta["_MultiIndex_levels_from_attrs"] == 1
    ):
        for k, v in db_meta["__multi_indices"].items():
            if str(k).startswith(multi_level_prefix):
                # index_name = str(k)[len(multi_level_prefix):]
                index_levels = v
                dataset = dataset.set_xindex(index_levels)

    # Apply variable metadata where available
    if "var_meta" in db_meta:
        vars_dict = db_meta["var_meta"]
        for varname, vardict in vars_dict.items():
            if "attrs" in vardict:
                var_attrs = vardict["attrs"]
                if varname == "dipoles" and (
                    "dip_perm" in dataset or "dip_trans" in dataset
                ):
                    # Dipoles should have been split back up
                    if "dip_perm" in dataset:
                        dataset["dip_perm"].attrs.update(var_attrs)
                    if "dip_trans" in dataset:
                        dataset["dip_trans"].attrs.update(var_attrs)
                else:
                    dataset[varname].attrs.update(var_attrs)

    # Set trajectory-level attributes
    if "misc_attrs" in db_meta:
        dataset.attrs.update(db_meta["misc_attrs"])

    # Perform a check of the dimension sizes specified in the metadata if present
    if "dims" in db_meta:
        for dimname, dimdict in db_meta["dims"].items():
            dim_length = dimdict["length"] if "length" in dimdict else -1
            if dim_length >= 0:
                if dim_length != dataset.sizes[dimname]:
                    msg = f"Size of dimension {dimname} in dataset parsed from ASE database has length inconsistent with metadata of ASE file. Was {dataset.sizes[dimname]} but metadata specifies {dim_length}"
                    logging.error(msg)
                    raise ValueError(msg)

    return dataset


def read_ase(
    db_path: pathlib.Path,
    kind: Literal['spainn', 'schnet'] | None = None,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset:
    """Reads an ASE DB containing data in the SPaiNN or SchNet format

    Parameters
    ----------
    db_path: pathlib.Path
        Path to the database
    kind: Literal['spainn', 'schnet'] | None, optional
        Must be one of 'spainn' or 'schnet' or None; determines interpretation of array shapes If None is provided, no shape will be assumed
    loading_parameters: LoadingParameters
        Potentially configured parameters to overwrite loading behavior

    Returns
    -------
        An `xr.Dataset` of frames

    Raises
    ------
    ValueError
        If `kind` is not one of 'spainn' or 'schnet'
    FileNotFoundError
        If `db_path` is not a file
    ValueError
        If `db_path` does not contain data corresponding to the format `kind`
    """

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Could not find databse at {db_path}")

    ase_default_attrs = get_default_input_attributes("ase", loading_parameters)

    with connect(db_path) as db:
        metadata = db.metadata
        shapes, coord_shapes, leading_dimension_name = shapes_from_metadata(
            metadata, kind
        )
        leading_dimension_rename_target = None

        data_vars = {}
        coord_vars = {}
        found_rows = 0
        available_varnames = next(db.select()).data.keys()
        print(available_varnames)

        tmp_data_in = {
            "atXYZ": [],
        }

        for row in db.select():
            for key, value in row.data.items():
                if key not in tmp_data_in:
                    tmp_data_in[key] = []

                tmp_data_in[key].append(value)
            # TODO: FIXME: deal with different atoms/compounds in the same DB.
            row_atoms = row.toatoms()
            if row_atoms.has("positions"):
                tmp_data_in['atXYZ'].append(row_atoms.get_positions())

            if row_atoms.has("symbols"):
                if 'atNames' not in tmp_data_in:
                    tmp_data_in['atNames'] = []
                tmp_data_in['atNames'].append(row_atoms.get_symbols())

            if row_atoms.has("momenta"):
                if 'velocities' not in tmp_data_in:
                    tmp_data_in['velocities'] = []
                tmp_data_in['velocities'].append(row_atoms.get_velocities())

            if "time" in row_atoms.info:
                if "time" not in tmp_data_in:
                    tmp_data_in["time"] = []
                tmp_data_in["time"].append(row_atoms.info["time"])

            if "trajid" in row_atoms.info:
                if "trajid" not in tmp_data_in:
                    tmp_data_in["trajid"] = []
                    leading_dimension_rename_target = "frame"

                tmp_data_in["trajid"].append(row_atoms.info["trajid"])

            found_rows += 1

    # If there are no valid rows, raise a ValueError
    if found_rows == 0:
        raise ValueError(
            f"No rows with the appropriate format for kind=`{kind}` were found in {db_path}"
        )

    for k, v in tmp_data_in.items():
        data_array = np.stack(v)
        if k in shapes:
            data_vars[k] = (
                shapes[k],
                data_array,
                (ase_default_attrs[k] if k in ase_default_attrs else None),
            )
        elif k in coord_shapes:
            if k == "atNames":
                data_array = data_array[0, :]
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

    frames = xr.Dataset(data_vars).assign_coords(coord_vars)
    if kind == 'spainn':
        assert 'tmp' in frames.dims
        frames = frames.squeeze('tmp')

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

    return apply_dataset_meta_from_db_metadata(frames, metadata)
