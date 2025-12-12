from numbers import Number
from typing import Literal

import numpy as np
import xarray as xr

from shnitsel.geo.geocalc import get_bond_lengths
from shnitsel.geo.geomatch import flag_bats_multiple
from shnitsel.bridges import default_mol
from shnitsel.clean.common import dispatch_cut
from shnitsel.units.conversion import convert_length

_default_bond_length_thresholds_angstrom = {'[#6,#7][H]': 2.0, '[*]~[*]': 3.0}


def _dict_to_thresholds(keys: list[str], d: dict, units: str) -> xr.DataArray:
    data = [d.get(c, np.nan) for c in keys]
    res = xr.DataArray(list(data), coords={'criterion': keys}, attrs={'units': units})
    return res.astype(float)


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


def bond_length_filtranda(
    frames, search_dict: dict[str, Number] | None = None, units='angstrom'
):
    if search_dict is None:
        search_dict = {}
    criteria = list(search_dict | _default_bond_length_thresholds_angstrom)
    default_thresholds = _dict_to_thresholds(
        criteria, _default_bond_length_thresholds_angstrom, units='angstrom'
    )
    default_thresholds = convert_length(default_thresholds, to=units)
    user_thresholds = _dict_to_thresholds(criteria, search_dict, units=units)
    thresholds = user_thresholds.where(~np.isnan(user_thresholds), default_thresholds)

    bonds = lengths_for_searches(
        convert_length(frames['atXYZ'], to=units),
        list(thresholds.coords['criterion'].data),
    )
    return (
        bonds.groupby('bond_search')
        .max()
        .rename({'bond_search': 'criterion'})
        .assign_coords({'thresholds': thresholds})
    )


def filter_by_length(
    frames,
    cut: Literal['truncate', 'omit', False] | Number = 'truncate',
    search_dict: dict[str, Number] | None = None,
):

    frames = frames.assign(
        filtranda=bond_length_filtranda(frames, search_dict=search_dict)
    )
    return dispatch_cut(frames, cut)