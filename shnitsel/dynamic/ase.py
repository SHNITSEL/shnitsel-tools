import os
from ase import Atoms
from ase.db import connect
import numpy as np
import xarray as xr
from typing import Collection


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


def write_ase_db(frames: xr.Dataset, db_path: str, keys: Collection | None = None):
    frames = _prepare_for_write(frames)

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


def read_ase_db(db_path: str):
    shapes = {
        'energy': ['frame', 'state'],
        'socs': ['frame', 'soc'],
        'forces': ['frame', 'state', 'atom', 'direction'],
        'nacs': ['frame', 'statecomb', 'atom', 'direction'],
        'dipoles': ['frame', 'state_or_statecomb', 'direction'],
    }

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
        data_vars['atNames'] = ['atom'], next(db.select()).symbols

    if 'dipoles' in data_vars:
        nstates = data_vars['energy'][1].shape[1]

        dipoles = data_vars['dipoles'][1]
        dip_perm = dipoles[:, :nstates, :]
        dip_trans = dipoles[:, nstates:, :]
        del data_vars['dipoles']

        data_vars['dip_perm'] = ['frame', 'state', 'direction'], dip_perm
        data_vars['dip_trans'] = ['frame', 'statecomb', 'direction'], dip_trans

    return xr.Dataset(data_vars)
