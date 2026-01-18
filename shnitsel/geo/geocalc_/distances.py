from typing import Any, overload
from shnitsel._contracts import needs
from shnitsel.core._api_info import API
from shnitsel.core.typedefs import AtXYZ

import xarray as xr

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import FeatureTypeLabel, StructureSelection
from .helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
    _get_default_selection,
)

from .algebra import dnorm


@API()
@needs(dims={'atom'})
def distance(atXYZ: AtXYZ, i: int, j: int) -> xr.DataArray:
    """Method to calculate the various distances between atoms i and j throughout time

    Args:
        atXYZ (AtXYZ): Array with atom positions
        i (int): Index of the first atom
        j (int): Index of the second atom

    Returns:
        xr.DataArray: The resulting array holding the pairwise distance between i and j.
    """
    a = atXYZ.isel(atom=i, drop=True)
    b = atXYZ.isel(atom=j, drop=True)
    with xr.set_options(keep_attrs=True):
        result: xr.DataArray = dnorm(a - b)
    # result.name = 'distance'
    # result.attrs['long_name'] = r"$\|\mathbf{r}_{%d,%d}\|$" % (i, j)
    return result


@overload
def get_distances(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray],
    structure_selection: StructureSelection | None = None,
) -> TreeNode[Any, xr.DataArray]: ...
@overload
def get_distances(
    atXYZ_source: Trajectory | Frames | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection | None = None,
) -> xr.DataArray: ...


@API()
@needs(dims={'atom', 'direction'})
def get_distances(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray]
    | Trajectory
    | Frames
    | xr.Dataset
    | xr.DataArray,
    structure_selection: StructureSelection | None = None,
) -> TreeNode[Any, xr.DataArray] | xr.DataArray:
    """Identify bonds (using RDKit) and find the length of each bond in each
    frame.

    Parameters
    ----------
    atXYZ_source
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions ``atom`` and
        ``direction`` or another source of positional data like a trajectory, a frameset,
        a dataset representing either of those or a tree structure holding such data.
    structure_selection, optional
        Object encapsulating feature selection on the structure whose positional information is provided in `atXYZ`.
        If this argument is omitted altogether, a default selection for all bonds within the structure is created.
    Returns
    -------
        An :py:class:`xarray.DataArray` of bond lengths/distances with dimension `descriptor` to index the distances along.

    Raises
    ------
    UserWarning
        If both `matches` and `mol` are specified.
    """
    if isinstance(atXYZ_source, TreeNode):
        return atXYZ_source.map_data(
            lambda x: get_distances(
                x,
                structure_selection=structure_selection,
            ),
            keep_empty_branches=True,
            dtype=xr.DataArray,
        )

    structure_selection = _get_default_selection(
        structure_selection, atXYZ_source=atXYZ_source, default_levels=['bonds']
    )

    atXYZ: xr.DataArray
    if isinstance(atXYZ_source, xr.DataArray):
        atXYZ = atXYZ_source
    else:
        atXYZ = atXYZ_source.atXYZ

    bond_indices = list(structure_selection.bonds_selected)

    if len(bond_indices) == 0:
        return _empty_descriptor_results(atXYZ)

    distance_arrs = [
        distance(atXYZ, a, b).expand_dims('descriptor') for a, b in bond_indices
    ]

    distance_res = xr.concat(distance_arrs, dim='descriptor')

    descriptor_tex = [r'|\vec{r}_{%d,%d}|' % (a, b) for a, b in bond_indices]
    descriptor_name = [r'dist(%d,%d)' % (a, b) for a, b in bond_indices]
    descriptor_type: list[FeatureTypeLabel] = ['dist'] * len(descriptor_tex)

    return _assign_descriptor_coords(
        distance_res,
        feature_descriptors=bond_indices,
        feature_type=descriptor_type,
        feature_tex_label=descriptor_tex,
        feature_name=descriptor_name,
    )
