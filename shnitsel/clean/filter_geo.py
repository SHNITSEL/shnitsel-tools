from numbers import Number
from typing import Literal, Sequence
import logging

import numpy as np
import xarray as xr
from rdkit.Chem import Mol

from shnitsel.geo.geocalc import get_bond_lengths
from shnitsel.geo.geomatch import flag_bats_multiple
from shnitsel.bridges import default_mol
from shnitsel.clean.common import dispatch_cut
from shnitsel.units.conversion import convert_length
from shnitsel.clean.dispatch_plots import dispatch_plots

_default_bond_length_thresholds_angstrom = {'[#6,#7][H]': 2.0, '[*]~[*]': 3.0}


def _dict_to_thresholds(keys: list[str], d: dict, units: str) -> xr.DataArray:
    data = [d.get(c, np.nan) for c in keys]
    res = xr.DataArray(list(data), coords={'criterion': keys}, attrs={'units': units})
    return res.astype(float)


def _lengths_for_searches(atXYZ, searches, mol=None):
    if mol is None:
        mol = default_mol(atXYZ)

    logging.disable(logging.INFO)
    matches = flag_bats_multiple(mol, searches)
    logging.disable(logging.NOTSET)

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
    frames,
    search_dict: dict[str, Number] | None = None,
    units='angstrom',
    mol: Mol | None = None,
):
    """Derive bond length filtration targets from an xr.Dataset

    Parameters
    ----------
    frames
        A xr.Dataset with an ``atXYZ`` variable
    search_dict, optional
        A mapping from SMARTS-strings to length-thresholds.

            - The SMARTS-strings describe bonds which are searched
              for in an RDKit Mol object obtained via :py:func:`shnitsel.bridges.default_mol`
            - The thresholds describe maximal tolerable bond-lengths; if there are multiple matches
              for a given search, the longest bond-length will be considered for each frame
    units, optional
        Units in which custom thresholds are given, and to which defaults and data will be converted, by default
        'angstrom'

    Returns
    -------
        An xr.DataArray of filtration targets stacked along the ``criterion`` dimension;
        one criterion per ``search_dict`` entry.
    """
    if search_dict is None:
        search_dict = {}
        criteria = list(_default_bond_length_thresholds_angstrom)
        default_thresholds = _dict_to_thresholds(
                criteria, _default_bond_length_thresholds_angstrom, units='angstrom'
                )
        default_thresholds = convert_length(default_thresholds, to=units)
        thresholds = default_thresholds
    else:
        criteria = list(search_dict.keys())
        user_thresholds = _dict_to_thresholds(criteria, search_dict, units=units)
        thresholds = user_thresholds
        #user_thresholds.where(~np.isnan(user_thresholds), default_thresholds)

    convert_coords = convert_length(frames['atXYZ'], to=units)

    bonds = _lengths_for_searches(
        convert_coords,
        list(thresholds.coords['criterion'].data),
        mol=mol,
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
    units: str = 'angstrom',
    plot_thresholds: bool | Sequence[float] = False,
    plot_populations: bool | Literal['independent', 'intersections'] = False,
    mol: Mol | None = None,
):
    """Filter trajectories according to bond length

    Parameters
    ----------
    frames
        A xr.Dataset with an ``atXYZ`` variable (NB. this function takes an xr.Dataset as
        opposed to an xr.DataArray for consistency with :py:func:`shnitsel.clean.sanity_check`)
    cut
        Specifies the manner in which to remove data;

            - if 'omit', drop trajectories unless all frames meet criteria (:py:func:`shnitsel.clean.omit`)
            - if 'truncate', cut each trajectory off just before the first frame that doesn't meet criteria
              (:py:func:`shnitsel.clean.truncate`)
            - if a number, interpret this number as a time, and cut all trajectories off at this time,
              discarding those which violate criteria before reaching the given limit,
              (:py:func:`shnitsel.clean.transect`)
            - if ``False``, merely annotate the data;
        see :py:func:`shnitsel.clean.dispatch_cut`.
    search_dict
        A mapping from SMARTS-strings to length-thresholds.

            - The SMARTS-strings describe bonds which are searched
              for in an RDKit Mol object obtained via :py:func:`shnitsel.bridges.default_mol`
            - The thresholds describe maximal tolerable bond-lengths; if there are multiple matches
              for a given search, the longest bond-length will be considered for each frame
    plot_thresholds
        See :py:func:`shnitsel.vis.plot.filtration.check_thresholds`.

        - If ``True``, will plot using ``check_thresholds`` with
        default quantiles
        - If a ``Sequence``, will plot using ``check_thresholds``
        with specified quantiles
        - If ``False``, will not plot threshold plot
    plot_populations
        See :py:func:`shnitsel.vis.plot.filtration.validity_populations`.

        - If ``True`` or ``'intersections'``, will plot populations of
        trajectories satisfying intersecting conditions
        - If ``'independent'``, will plot populations of
        trajectories satisfying conditions taken independently
        - If ``False``, will not plot populations plot
    mol
        An rdkit mol object, if not provided it will be generated from the XYZ coordinates in the data
    units
        Units in which custom thresholds are given, and to which defaults and data will be converted, by default
        'angstrom'

    Returns
    -------
        The filtered Dataset

    Notes
    -----
    The resulting object has a ``filtranda`` data_var, representing the values by which the data were filtered.
    If the input has a ``filtranda`` data_var, it is overwritten.
    """
    filtranda = bond_length_filtranda(
        frames, search_dict=search_dict, units=units, mol=mol
    )
    frames = frames.drop_dims(['criterion'], errors='ignore').assign(
        filtranda=filtranda
    )
    dispatch_plots(filtranda, plot_thresholds, plot_populations)

    return dispatch_cut(frames, cut)


