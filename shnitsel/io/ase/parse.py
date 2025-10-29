import os
from typing import Literal

from ase.db import connect
import numpy as np
import xarray as xr


def parse_ase(db_path: str, kind: Literal['spainn', 'schnet'] | None) -> xr.Dataset:
    """Reads an ASE DB containing data in the SPaiNN or SchNet format

    Parameters
    ----------
    db_path
        Path to the database
    kind
        Must be one of 'spainn' or 'schnet' or None; determines interpretation of array shapes If None is provided, no shape will be assumed

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
    if kind == 'schnet':
        shapes = {
            'energy': ['frame', 'state'],
            'forces': ['frame', 'state', 'atom', 'direction'],
            'nacs': ['frame', 'statecomb', 'atom', 'direction'],
            'dipoles': ['frame', 'state_or_statecomb', 'direction'],
        }
    elif kind == 'spainn':
        shapes = {
            # Note the extra dim, removed below
            'energy': ['frame', 'tmp', 'state'],
            'forces': ['frame', 'atom', 'state', 'direction'],
            'nacs': ['frame', 'atom', 'statecomb', 'direction'],
            'dipoles': ['frame', 'state_or_statecomb', 'direction'],
        }
    elif kind is None:
        shapes = {}
    else:
        raise ValueError(
            f"'kind' should be one of 'schnet' or 'spainn', not '{kind}'")

    if not os.path.isfile(db_path):
        raise FileNotFoundError(db_path)

    with connect(db_path) as db:
        data_vars = {}
        found_rows = 0
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
                f"No rows with the appropriate format for kind={kind} were found in {db_path}")

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
