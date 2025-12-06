from numbers import Number

import numpy as np
import xarray as xr

from shnitsel.core.geom import get_bond_lengths
from shnitsel.geo.geomatch import flag_bats_multiple
from shnitsel.bridges import default_mol
from shnitsel.core.filtration2 import omit, truncate, transect


def lengths_for_searches(atXYZ, searches):
    mol = default_mol(atXYZ)
    matches = flag_bats_multiple(mol, searches)
    bonds = xr.concat(
        [
            (x := get_bond_lengths(atXYZ, v)).assign_coords(
                {'bond_search': ('descriptor', np.full(x.sizes['descriptor'], k))}
            )
            for k, v in matches.items()
        ],
        dim='descriptor',
    )
    return bonds


def bond_length_filtranda(frames, search_dict):
    bonds = lengths_for_searches(frames['atXYZ'], list(search_dict))
    return (
        bonds.groupby('bond_search')
        .max()
        .rename({'bond_search': 'criterion'})
        .assign_coords({'thresholds': ('criterion', list(search_dict.values()))})
    )


def filter_by_length(
    frames,
    search_dict,
    cut=False,
):
    frames = frames.assign(
        filtranda=bond_length_filtranda(frames, search_dict=search_dict)
    )
    if cut == 'truncate':
        return truncate(frames)
    elif cut == 'omit':
        return omit(frames)
    else:
        assert isinstance(cut, Number)
        return transect(frames, cut)