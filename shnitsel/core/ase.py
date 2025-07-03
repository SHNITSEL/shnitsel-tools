import os
from typing import Collection

from ase import Atoms
from ase.db import connect
import numpy as np
import xarray as xr


def _prepare_for_write(frames: xr.Dataset) -> xr.Dataset:
    # Recombine permanent and transition dipoles, as schnetpack expects
    dipoles: np.ndarray | xr.DataArray | None = None
    frames = frames.copy(deep=False)
    if 'dipoles' in frames:
        dipoles = frames['dipoles']
    elif 'dip_perm' in frames and 'dip_trans' in frames:
        dip_perm = frames['dip_perm'].transpose('frame', 'state', 'direction').data
        dip_trans = (
            frames['dip_trans'].transpose('frame', 'statecomb', 'direction').data
        )
        dipoles = np.concat((dip_perm, dip_trans.data), axis=1)
        del frames['dip_perm'], frames['dip_trans']
    elif 'dip_perm' in frames:
        dipoles = frames['dip_perm']
        del frames['dip_perm']
    elif 'dip_trans' in frames:
        dipoles = frames['dip_trans']
        del frames['dip_trans']

    if dipoles is not None:
        frames['dipoles'] = ['frame', 'state_or_statecomb', 'direction'], dipoles

    return frames


def write_ase(
    frames: xr.Dataset,
    db_path: str,
    kind: str | None,
    keys: Collection | None = None,
    preprocess: bool = True,
):
    if preprocess:
        frames = _prepare_for_write(frames)

    statedims = ['state', 'statecomb', 'state_or_statecomb']
    if kind == 'schnet':
        order = ['frame', *statedims, 'atom', 'direction']
        frames = frames.transpose(*order, missing_dims='ignore')
    elif kind == 'spainn':
        frames['energy'] = frames['energy'].expand_dims('tmp', axis=1)
        order = ['frame', 'tmp', 'atom', *statedims, 'direction']
        frames = frames.transpose(*order, missing_dims='ignore')
    elif kind is None:
        # leave the axis orders as they are
        pass
    else:
        raise ValueError(
            f"'kind' should be one of 'schnet', 'spainn' or None, not '{kind}'"
        )

    if os.path.exists(db_path):
        os.remove(db_path)

    if not keys:
        keys = frames.data_vars.keys()
    keys = set(frames.data_vars).intersection(keys).difference({'atNames'})

    with connect(db_path, type='db') as db:
        for i, frame in frames.groupby('frame'):
            frame = frame.squeeze('frame')
            db.write(
                Atoms(symbols=frame['atNames'].data, positions=frame['atXYZ']),
                data={k: frame[k].data for k in keys},
            )


def read_ase(db_path: str, kind: str):
    """Reads an ASE DB containing data in the SPaiNN or SchNet format

    Parameters
    ----------
    db_path
        Path to the database
    kind
        Must be one of 'spainn' or 'schnet'; determines interpretation of array shapes

    Returns
    -------
        An `xr.Dataset` of frames

    Raises
    ------
    ValueError
        If `kind` is not one of 'spainn' or 'schnet'
    FileNotFoundError
        If `db_path` is not a file
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
            'energy': ['frame', 'tmp', 'state'],  # Note the extra dim, removed below
            'forces': ['frame', 'atom', 'state', 'direction'],
            'nacs': ['frame', 'atom', 'statecomb', 'direction'],
            'dipoles': ['frame', 'state_or_statecomb', 'direction'],
        }
    else:
        raise ValueError(f"'kind' should be one of 'schnet' or 'spainn', not '{kind}'")

    if not os.path.isfile(db_path):
        raise FileNotFoundError(db_path)

    with connect(db_path) as db:
        data_vars = {}
        for name, dims in shapes.items():
            try:
                data = np.stack([row.data[name] for row in db.select()])
                data_vars[name] = dims, data
            except KeyError:
                pass

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
