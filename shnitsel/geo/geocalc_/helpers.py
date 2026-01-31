from typing import Literal, Sequence, TypeAlias

import numpy as np

from shnitsel.core.typedefs import AtXYZ
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.filtering.structure_selection import (
    FeatureDescriptor,
    FeatureLevelOptions,
    FeatureTypeLabel,
    StructureSelection,
    StructureSelectionDescriptor,
)
import xarray as xr


def _empty_descriptor_results(atXYZ: AtXYZ) -> xr.DataArray:
    """Get an empty result with no 'atom' dimension but a zero-length 'descriptor' dimension in its stead.

    Parameters
    ----------
    atXYZ : AtXYZ
        The positional input to model the descriptor result after (shape- and coordinate-wise)

    Returns
    -------
    xr.DataArray
        The empty descriptor DataArray of the desired shape.
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

    Parameters
    ----------
    obj : xr.DataArray
        The object to assign the new coordinates to. Should have a 'descriptor' dimension.
    feature_descriptors : Sequence[FeatureDescriptor]
        The indices of the atoms involved as a sequence of arbitrary tuples (see `FeatureDescriptor`).
    feature_type : Sequence[FeatureTypeLabel]
        The list of feature types.
    feature_tex_label : Sequence[str]
        Tex labels for the respective feature for more elaborate output.
    feature_name : Sequence[str]
        A simple ascii type name to be able to lookup features easier.

    Returns
    -------
    xr.DataArray
        The resulting DataArray with the new coordinates assigned.
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

    return (
        obj.assign_coords(coords)
        .set_xindex('feature_indices')
        .set_xindex('descriptor_type')
    )  # .set_xindex('descriptor')


# We export this now
assign_descriptor_coords = _assign_descriptor_coords


def _remove_atom_coords(da: xr.DataArray) -> xr.DataArray:
    """Helper function to remove all standard atom-related coordinates.

    Parameters
    ----------
    da : xr.DataArray
        The data array to clean of all the atom-related coordinates.

    Returns
    -------
    xr.DataArray
        The cleaned array.
    """

    if 'atom' in da.dims:
        da = da.squeeze('atom', drop=True)

    atom_keys = ['atNames', 'atNums', 'atom']

    atom_keys = set(da.coords.keys()).intersection(atom_keys)

    return da.reset_coords(atom_keys, drop=True)
