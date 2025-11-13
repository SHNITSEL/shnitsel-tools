import os
from typing import Collection, Literal

from ase import Atoms
from ase.db import connect
import numpy as np
import xarray as xr

from shnitsel._contracts import needs


def _prepare_for_write(frames: xr.Dataset) -> xr.Dataset:
    # Recombine permanent and transition dipoles, as schnetpack expects
    dipoles: np.ndarray | xr.DataArray | None = None
    frames = frames.copy(deep=False)
    if 'dipoles' in frames:
        dipoles = frames['dipoles']
    elif 'dip_perm' in frames and 'dip_trans' in frames:
        dip_perm = frames['dip_perm'].transpose(
            'frame', 'state', 'direction').data
        dip_trans = (
            frames['dip_trans'].transpose(
                'frame', 'statecomb', 'direction').data
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
        frames['dipoles'] = [
            'frame', 'state_or_statecomb', 'direction'], dipoles

    return frames


@needs(dims={'frame'})
def write_ase_db(
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
