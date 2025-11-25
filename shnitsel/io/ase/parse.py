import logging
import os
import pathlib
from typing import Literal

from ase.db import connect
import numpy as np
import xarray as xr

from shnitsel.io.helpers import LoadingParameters
from shnitsel.units.defaults import get_default_input_attributes


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
    schnet_shapes = {
        'energy': ['frame', 'state'],
        'forces': ['frame', 'state', 'atom', 'direction'],
        'nacs': ['frame', 'statecomb', 'atom', 'direction'],
        'smooth_nacs': ['frame', 'statecomb', 'atom', 'direction'],
        'socs': ['frame', 'statecomb'],
        'dipoles': ['frame', 'state_or_statecomb', 'direction'],
    }
    spainn_shapes = {
        # Note the extra dim, removed below
        'energy': ['frame', 'tmp', 'state'],
        'forces': ['frame', 'atom', 'state', 'direction'],
        'nacs': ['frame', 'atom', 'statecomb', 'direction'],
        'smooth_nacs': ['frame', 'atom', 'statecomb', 'direction'],
        'socs': ['frame', 'statecomb'],
        'dipoles': ['frame', 'state_or_statecomb', 'direction'],
    }

    if kind == 'schnet':
        shapes = schnet_shapes
    elif kind == 'spainn':
        shapes = spainn_shapes
    elif kind is None:
        shapes = None
    else:
        raise ValueError(f"'kind' should be one of 'schnet' or 'spainn', not '{kind}'")

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Could not find databse at {db_path}")

    ase_default_attrs = get_default_input_attributes("ase", loading_parameters)

    with connect(db_path) as db:
        metadata = db.metadata

        if "db_format" in metadata:
            meta_format = metadata["db_format"]
            if meta_format not in ["schnet", "spainn"]:
                raise ValueError(
                    f"Database {db_path} is of unsupported format: {meta_format}. Only `schnet` and `spainn` are supported."
                )

            if kind is None:
                kind = meta_format
                logging.info(f"Automatically detected format: {kind}")
                if kind == 'schnet':
                    shapes = schnet_shapes
                elif kind == 'spainn':
                    shapes = spainn_shapes

            if meta_format != kind:
                raise ValueError(
                    f"Database {db_path} is of format: {meta_format} instead of requested format {kind}."
                )
           

        data_vars = {}
        found_rows = 0
        available_varnames = next(db.select()).data.keys()
        print(available_varnames)
        raise ValueError()

        for name, dims in shapes.items():
            try:
                data = np.stack([row.data[name] for row in db.select()])

                data_vars[name] = dims, data
                found_rows += 1
            except KeyError:
                pass

        # If there are no valid rows, raise a ValueError
        if found_rows == 0:
            raise ValueError(
                f"No rows with the appropriate format for kind={kind} were found in {db_path}"
            )

        atXYZ = np.stack([row.positions for row in db.select()])
        data_vars['atXYZ'] = ['frame', 'atom', 'direction'], atXYZ
        atNames = ['atom'], next(db.select()).symbols

    if 'dipoles' in data_vars:
        nstates = data_vars['energy'][1].shape[1]

        dipoles = data_vars['dipoles'][1]
        dip_perm = dipoles[:, :nstates, :]
        dip_trans = dipoles[:, nstates:, :]
        del data_vars['dipoles']

        data_vars['dip_perm'] = ['frame', 'state', 'direction'], dip_perm
        data_vars['dip_trans'] = ['frame', 'statecomb', 'direction'], dip_trans

    frames = xr.Dataset(data_vars).assign_coords(atNames=atNames)
    if kind == 'spainn':
        assert 'tmp' in frames.dims
        frames = frames.squeeze('tmp')

    return frames
