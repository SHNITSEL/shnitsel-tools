from typing import Literal, Sequence, TypeAlias

import numpy as np

from shnitsel.core.typedefs import AtXYZ
from shnitsel.filtering.structure_selection import (
    FeatureDescriptor,
    FeatureLevelOptions,
    FeatureTypeLabel,
    StructureSelection,
)
import rdkit.Chem as rc
import xarray as xr
from shnitsel.bridges import construct_default_mol


def _get_default_selection(
    structure_selection: StructureSelection | None = None,
    mol: rc.Mol | None = None,
    atXYZ: xr.DataArray | None = None,
    default_levels: Sequence[FeatureLevelOptions] = ['atoms', 'bonds'],
) -> StructureSelection:
    """Get a default selection object from any accessible data if possible.

    Args:
        structure_selection (StructureSelection | None, optional): A potential already provided structure/feature selection. Defaults to None. If provided, will be returned back.
        mol (rc.Mol | None, optional): An optional instance of an RDKit molecule. Used to construct a StructureSelection instance if `structure_selection` is None. Defaults to None.
        atXYZ (xr.DataArray | None, optional): An xr.DataArray holding positional data of a molecule. Only the first frame will be used to create a default mol if `structure_selection` and `mol` were not provided. Defaults to None.
        default_levels (Sequence[FeatureLevelOptions], optional): the desired default levels included in the selection that may be recreated if none was provided. Defaults to 'atoms' and 'bonds'.

    Raises:
        ValueError: _description_

    Returns:
        StructureSelection: The initialized default structure selection.
    """
    if structure_selection is not None:
        return structure_selection

    if mol is not None and isinstance(mol, rc.Mol):
        return StructureSelection.init_from_mol(mol, default_selection=default_levels)
    elif mol is None and atXYZ is not None:
        mol = construct_default_mol(atXYZ)
    elif mol is None:
        raise ValueError(
            "You did not provide sufficient data to construct a default feature selection. Please provide your own StructureSelection object."
        )

    return StructureSelection.init_from_mol(mol, default_selection=default_levels)


def _empty_descriptor_results(atXYZ: AtXYZ) -> xr.DataArray:
    """Get an empty result with no 'atom' dimension but a zero-length 'descriptor' dimension in its stead.

    Args:
        atXYZ (AtXYZ): The positional input to model the descriptor result after (shape- and coordinate-wise)

    Returns:
        xr.DataArray: The empty descriptor DataArray of the desired shape.
    """
    res = atXYZ.isel(atom=[], direction=0, drop=True).rename({'atom': 'descriptor'})
    return _assign_descriptor_coords(res, [], [], [], [])


def _assign_descriptor_coords(
    obj: xr.DataArray,
    feature_descriptors: Sequence[FeatureDescriptor],
    feature_type: Sequence[FeatureTypeLabel],
    feature_tex_label: Sequence[str],
    feature_name: Sequence[str],
) -> xr.DataArray:
    """Helper function to assign all descriptor coordinates to keys of the same name

    Args:
        obj (xr.DataArray): The object to assign the new coordinates to. Should have a 'descriptor' dimension.
        feature_descriptors (Sequence[FeatureDescriptor]): The indices of the atoms involved as a sequence of arbitrary tuples (see `FeatureDescriptor`).
        feature_type (Sequence[FeatureTypeLabel]): The list of feature types.
        feature_tex_label (Sequence[str]): Tex labels for the respective feature for more elaborate output.
        feature_name (Sequence[str]): A simple ascii type name to be able to lookup features easier.

    Returns:
        xr.DataArray: The resulting DataArray with the new coordinates assigned.
    """

    feature_descriptors_np = np.empty((len(feature_descriptors),), dtype='O')
    feature_descriptors_np[:] = feature_descriptors
    coords = xr.Coordinates(
        {
            'descriptor': ('descriptor', feature_name),
            'descriptor_tex': ('descriptor', feature_tex_label),
            'descriptor_type': ('descriptor', feature_type),
            'feature_indices': ('descriptor', feature_descriptors_np),
        }
    )

    return obj.assign_coords(coords)  # .set_xindex('descriptor')


def _remove_atom_coords(da: xr.DataArray) -> xr.DataArray:
    """Helper function to remove all standard atom-related coordinates.

    Args:
        da (xr.DataArray): The data array to clean of all the atom-related coordinates.

    Returns:
        xr.DataArray: The cleaned array.
    """

    if 'atom' in da.dims:
        da = da.squeeze('atom', drop=True)

    atom_keys = ['atNames', 'atNums', 'atom']

    atom_keys = set(da.coords.keys()).intersection(atom_keys)

    return da.reset_coords(atom_keys, drop=True)
