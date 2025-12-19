"""\
    This module will contain two types of closely-related functionality:
    - generic geometry functions
    - functions that use RDKit to identify internal coordinates, and then the above to calculate the values
    Currently, the first category of function is found under postprocess
"""

from typing import Literal, Iterable, Sequence


import numpy as np
import rdkit.Chem as rc
from rdkit.Chem import Mol
import xarray as xr

from shnitsel.core._api_info import API, internal
from shnitsel.filtering.structure_selection import (
    FeatureLevelType,
    StructureSelection,
)
from ..analyze.generic import norm
from ..bridges import construct_default_mol
from .._contracts import needs

from ..core.typedefs import AtXYZ

from shnitsel.geo.geomatch import flag_bats, flag_bla_chromophor

from .geocalc_.positions import get_positions
from .geocalc_.distances import get_distances
from .geocalc_.angles import get_angles
from .geocalc_.dihedrals import get_dihedrals
from .geocalc_.pyramids import get_pyramidalization


@API()
@needs(dims={'atom', 'direction'})
def get_bats(
    atXYZ: AtXYZ,
    structure_selection: StructureSelection | None = None,
    default_features: Sequence[FeatureLevelType] = ['bonds', 'angles', 'dihedrals'],
    signed: bool = False,
    deg: bool | Literal['trig'] = True,
) -> xr.DataArray:
    """Get bond lengths, angles and torsions/dihedrals.

    Parameters
    ----------
    atXYZ
        The positional data of atoms to use.
    structure_selection, optional: StructureSelection
        A feature selection to use. Can specify which features (positions, distances,
        angles, torsions or pyramidalizations) to include in the result.
        If not set, will be initialized to a default selection of all molecule-internal features
        as specified by the structure in the first frame of `atXYZ` and the features
        listed in `default_features`.
    default_features, optional: Sequence[FeatureLevelType]
        If no `structure_selection` object is provided, will select all features of these levels
        within the structure encoded in `atXYZ`.
        Options are
        - `atoms` for positional data,
        - `bonds` for distances between pairs of atoms (defaults to only bonds)
        - `angles` for angles between pairs of bonds between atoms.
        - `dihedrals` for torsion angles of bonds
        - `pyramids` for pyramidalization angles in the molecule.
        Defaults to using bonds, angles and dihedrals/torsions.
    signed, optional
        Whether to distinguish between clockwise and anticlockwise rotation,
        when returning angles as opposed to cosine & sine values;
        by default, do not distinguish.
        NB. This applies only to the dihedrals, not to the three-center angles.
        The latter are always unsigned.
    deg, optional: bool or Literal['trig']
        If True (the default), returns angles in degrees.
        If False, returns angles in radians.
        If set to 'trig' returns sines and cosines;

    Returns
    -------
        An :py:class:`xarray.DataArray` containing bond lengths, angles and tensions.

    Examples
    --------
        >>> import shnitsel as st
        >>> from shnitsel.core import geom
        >>> frames = st.open_frames('/tmp/A03_filtered.nc')
        >>> geom.get_bats(frames['atXYZ'])
    """
    # TODO: FIXME: Example is not up to date
    if matches_or_mol is None:
        mol = default_mol(atXYZ)
        matches_or_mol = flag_bats(mol)[0]

    d = {
        'bonds': get_bond_lengths(atXYZ, matches_or_mol=matches_or_mol),
        'angles': get_bond_angles(atXYZ, matches_or_mol=matches_or_mol, ang=ang),
        'dihedrals': get_bond_torsions(
            atXYZ, matches_or_mol=matches_or_mol, signed=signed, ang=ang
        ),
    }

    if pyr:
        d['pyr'] = get_pyramids(
            atXYZ, matches_or_mol=matches_or_mol, ang=ang, signed=signed
        )

    for k in d:
        d[k] = d[k].drop_vars(
            ['atom0', 'atom1', 'atom2', 'atom3', 'atNames', 'atNums'], errors='ignore'
        )

    return xr.concat(list(d.values()), dim='descriptor')
