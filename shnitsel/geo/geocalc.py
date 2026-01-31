"""\
    This module contains functionality to calculate certain geometric features from positional
    data found within datasets. 
    Through the use of structural selection support from the `shnitsel.filtration` module,
    specifically `StructureSelection`, the set of featueres to calculate can be restricted.
"""

import logging
from typing import Any, Literal, Sequence, overload


import xarray as xr

from shnitsel.core._api_info import API
from shnitsel.data.dataset_containers import Frames, Trajectory, wrap_dataset
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree import ShnitselDB

from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import (
    FeatureLevelType,
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc_.helpers import (
    _empty_descriptor_results,
)
from shnitsel.filtering.helpers import _get_default_structure_selection

from .._contracts import needs

from ..core.typedefs import AtXYZ

from .geocalc_.positions import get_positions
from .geocalc_.distances import get_distances
from .geocalc_.angles import get_angles
from .geocalc_.dihedrals import get_dihedrals
from .geocalc_.pyramids import get_pyramidalization
from .geocalc_.bla_chromophor import get_max_chromophor_BLA
from .alignment import get_centered_geometry, kabsch

__all__ = [
    "get_bats",
    "get_positions",
    "get_distances",
    "get_angles",
    "get_dihedrals",
    "get_pyramidalization",
    "get_max_chromophor_BLA",
    "get_centered_geometry",
    "kabsch",
]


@overload
def get_bats(
    atXYZ: ShnitselDataset | xr.Dataset | AtXYZ,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    default_features: Sequence[FeatureLevelType] = ['bonds', 'angles', 'dihedrals'],
    signed: bool = False,
    deg: bool | Literal['trig'] = True,
) -> xr.DataArray: ...


@overload
def get_bats(
    atXYZ: ShnitselDB[ShnitselDataset | xr.Dataset | AtXYZ],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    default_features: Sequence[FeatureLevelType] = ['bonds', 'angles', 'dihedrals'],
    signed: bool = False,
    deg: bool | Literal['trig'] = True,
) -> ShnitselDB[xr.DataArray]: ...


@API()
@needs(dims={'atom', 'direction'})
def get_bats(
    atXYZ: ShnitselDataset
    | xr.Dataset
    | TreeNode[Any, ShnitselDataset | xr.Dataset]
    | AtXYZ,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    default_features: Sequence[FeatureLevelType] = ['bonds', 'angles', 'dihedrals'],
    signed: bool = False,
    deg: bool | Literal['trig'] = True,
) -> xr.DataArray | TreeNode[Any, xr.DataArray]:
    """Get bond lengths, angles and torsions/dihedrals.

    Parameters
    ----------
    atXYZ : Trajectory | Frames | TreeNode[Any, Trajectory | Frames] | AtXYZ
        The positional data of atoms to use.
    structure_selection: StructureSelection | StructureSelectionDescriptor, optional
        A feature selection to use. Can specify which features (positions, distances,
        angles, torsions or pyramidalizations) to include in the result.
        If not set, will be initialized to a default selection of all molecule-internal features
        as specified by the structure in the first frame of `atXYZ` and the features
        listed in `default_features`.
    default_features : Sequence[FeatureLevelType], optional
        If no `structure_selection` object is provided, will select all features of these levels
        within the structure encoded in `atXYZ`.
        Options are
        - `atoms` for positional data,
        - `bonds` for distances between pairs of atoms (defaults to only bonds)
        - `angles` for angles between pairs of bonds between atoms.
        - `dihedrals` for torsion angles of bonds
        - `pyramids` for pyramidalization angles in the molecule.
        Defaults to using bonds, angles and dihedrals/torsions.
    signed: bool, optional
        Whether to distinguish between clockwise and anticlockwise rotation,
        when returning angles as opposed to cosine & sine values;
        by default, do not distinguish.
        NB. This applies only to the dihedrals, not to the three-center angles.
        The latter are always unsigned.
    deg : bool or Literal['trig'], optional
        If True (the default), returns angles in degrees.
        If False, returns angles in radians.
        If set to 'trig' returns sines and cosines;

    Returns
    -------
    xr.DataArray | TreeNode[Any, xr.DataArray]
        An :py:class:`xarray.DataArray` containing bond lengths, angles and tensions.

    Examples
    --------
        >>> import shnitsel as st
        >>> from shnitsel.geo import geocalc
        >>> data = st.read('/test_data/shnitsel/traj_I02.nc')
        >>> geocalc.get_bats(data)
    """

    if isinstance(atXYZ, TreeNode):
        return atXYZ.map_data(
            lambda x: get_bats(
                x,
                structure_selection=structure_selection,
                default_features=default_features,
                signed=signed,
                deg=deg,
            ),
            dtype=xr.DataArray,
        )
    else:
        # TODO: Fix other geocalc structure selection extractions once this works for partial submodules
        position_data: xr.DataArray
        position_source: xr.DataArray | ShnitselDataset
        charge_info: int | None
        if isinstance(atXYZ, xr.DataArray):
            position_data = atXYZ
            position_source = atXYZ
            charge_info = None
        else:
            wrapped_ds = wrap_dataset(atXYZ, (Trajectory | Frames))
            position_source = wrapped_ds
            position_data = wrapped_ds.atXYZ
            charge_info = int(wrapped_ds.charge)

        # TODO: FIXME: Example is not up to date
        structure_selection = _get_default_structure_selection(
            structure_selection,
            atXYZ_source=position_source,
            default_levels=default_features,
            charge_info=charge_info,
        )

        feature_data: list[xr.DataArray] = []
        if len(structure_selection.atoms_selected) > 0:
            feature_data.append(
                get_positions(position_data, structure_selection=structure_selection)
            )
        if len(structure_selection.bonds_selected) > 0:
            feature_data.append(
                get_distances(position_data, structure_selection=structure_selection)
            )
        if len(structure_selection.angles_selected) > 0:
            feature_data.append(
                get_angles(
                    position_data,
                    structure_selection=structure_selection,
                    deg=deg,
                    signed=signed,
                )
            )
        if len(structure_selection.dihedrals_selected) > 0:
            feature_data.append(
                get_dihedrals(
                    position_data,
                    structure_selection=structure_selection,
                    deg=deg,
                    signed=signed,
                )
            )
        if len(structure_selection.pyramids_selected) > 0:
            feature_data.append(
                get_pyramidalization(
                    position_data,
                    structure_selection=structure_selection,
                    deg=deg,
                    signed=signed,
                )
            )

        if len(feature_data) > 0:
            tmp_res = xr.concat(
                feature_data, dim='descriptor', combine_attrs='drop_conflicts'
            )
        else:
            logging.warning(
                "No feature data could be calculated. Did you provide an empty selection?"
            )
            tmp_res = _empty_descriptor_results(position_data)

        tmp_res.name = "BATs(+P)"
        tmp_res.attrs["long_name"] = (
            "Positions, bonds, angles, tortions, and pyramidalizations"
        )

        return tmp_res
