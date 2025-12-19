"""\
    This module will contain two types of closely-related functionality:
    - generic geometry functions
    - functions that use RDKit to identify internal coordinates, and then the above to calculate the values
    Currently, the first category of function is found under postprocess
"""

import logging
from typing import Literal, Sequence


import xarray as xr

from shnitsel.core._api_info import API
from shnitsel.filtering.structure_selection import (
    FeatureLevelType,
    StructureSelection,
)
from shnitsel.geo.geocalc_.helpers import (
    _empty_descriptor_results,
    _get_default_selection,
)
from .._contracts import needs

from ..core.typedefs import AtXYZ

from .geocalc_.positions import get_positions
from .geocalc_.distances import get_distances
from .geocalc_.angles import get_angles
from .geocalc_.dihedrals import get_dihedrals
from .geocalc_.pyramids import get_pyramidalization
from .geocalc_.bla_chromophor import get_max_chromophor_BLA


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
    structure_selection = _get_default_selection(
        structure_selection, atXYZ=atXYZ, default_levels=default_features
    )

    feature_data: list[xr.DataArray] = []
    if len(structure_selection.atoms_selected) > 0:
        feature_data.append(
            get_positions(atXYZ, structure_selection=structure_selection)
        )
    if len(structure_selection.bonds_selected) > 0:
        feature_data.append(
            get_distances(atXYZ, structure_selection=structure_selection)
        )
    if len(structure_selection.angles_selected) > 0:
        feature_data.append(
            get_angles(
                atXYZ, structure_selection=structure_selection, deg=deg, signed=signed
            )
        )
    if len(structure_selection.dihedrals_selected) > 0:
        feature_data.append(
            get_dihedrals(
                atXYZ,
                structure_selection=structure_selection,
                deg=deg if isinstance(deg, bool) else True,
                signed=signed,
            )
        )
    if len(structure_selection.pyramids_selected) > 0:
        feature_data.append(
            get_pyramidalization(
                atXYZ, structure_selection=structure_selection, signed=signed
            )
        )

    if len(feature_data) > 0:
        return xr.concat(feature_data, dim='descriptor')
    else:
        logging.warning(
            "No feature data could be calculated. Did you provide an empty selection?"
        )
        return _empty_descriptor_results(atXYZ)
