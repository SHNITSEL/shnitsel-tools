from typing import Any, overload
from shnitsel._contracts import needs
from shnitsel.core._api_info import API
from shnitsel.core.typedefs import AtXYZ

import xarray as xr

from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import (
    FeatureTypeLabel,
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
)
from shnitsel.filtering.helpers import _get_default_structure_selection

from shnitsel.geo.geocalc_.algebra import dnorm


@API()
@needs(dims={'atom'})
def distance(atXYZ: AtXYZ, i: int, j: int) -> xr.DataArray:
    """Method to calculate the various distances between atoms i and j throughout time

    Parameters
    ----------
    atXYZ : AtXYZ
        Array with atom positions
    i : int
        Index of the first atom
    j : int
        Index of the second atom

    Returns
    -------
    xr.DataArray
        The resulting array holding the pairwise distance between i and j.
    """
    if isinstance(atXYZ, TreeNode):
        return atXYZ.map_data(
            distance,
            i=i,
            j=j,
        )

    a = atXYZ.sel(atom=i, drop=True)
    b = atXYZ.sel(atom=j, drop=True)
    with xr.set_options(keep_attrs=True):
        result: xr.DataArray = dnorm(a - b)
    result.name = 'distance'
    result.attrs['long_name'] = r"\|\mathbf{r}_{%d,%d}\|" % (i, j)
    return result


@overload
def get_distances(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
) -> TreeNode[Any, xr.DataArray]: ...
@overload
def get_distances(
    atXYZ_source: Trajectory | Frames | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
) -> xr.DataArray: ...


@API()
@needs(dims={'atom', 'direction'})
def get_distances(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray]
    | Trajectory
    | Frames
    | xr.Dataset
    | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
) -> TreeNode[Any, xr.DataArray] | xr.DataArray:
    """Identify bonds (using RDKit) and find the length of each bond in each
    frame.

    Parameters
    ----------
    atXYZ_source : TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray] | Trajectory | Frames | xr.Dataset | xr.DataArray
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions ``atom`` and
        ``direction`` or another source of positional data like a trajectory, a frameset,
        a dataset representing either of those or a tree structure holding such data.
    structure_selection: StructureSelection | StructureSelectionDescriptor, optional
        Object encapsulating feature selection on the structure whose positional information is provided in `atXYZ`.
        If this argument is omitted altogether, a default selection for all bonds within the structure is created.

    Returns
    -------
    TreeNode[Any, xr.DataArray] | xr.DataArray
        An :py:class:`xarray.DataArray` of bond lengths/distances with dimension `descriptor` to index the distances along.
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

    position_data: xr.DataArray
    charge_info: int | None
    if isinstance(atXYZ_source, xr.DataArray):
        position_data = atXYZ_source
        charge_info = None
    else:
        wrapped_ds = wrap_dataset(atXYZ_source, (Trajectory | Frames))
        position_data = wrapped_ds.atXYZ
        charge_info = int(wrapped_ds.charge)

    structure_selection = _get_default_structure_selection(
        structure_selection,
        atXYZ_source=position_data,
        default_levels=['bonds'],
        charge_info=charge_info,
    )

    bond_indices = list(structure_selection.bonds_selected)

    if len(bond_indices) == 0:
        return _empty_descriptor_results(position_data)

    distance_arrs = [
        distance(position_data, a, b).expand_dims('descriptor') for a, b in bond_indices
    ]

    distance_res = xr.concat(distance_arrs, dim='descriptor')

    descriptor_tex = [r'|\vec{r}_{%d,%d}|' % (a, b) for a, b in bond_indices]
    descriptor_name = [r'dist(%d,%d)' % (a, b) for a, b in bond_indices]
    descriptor_type: list[FeatureTypeLabel] = ['dist'] * len(descriptor_tex)

    distance_res.name = "distances"

    return _assign_descriptor_coords(
        distance_res,
        feature_descriptors=bond_indices,
        feature_type=descriptor_type,
        feature_tex_label=descriptor_tex,
        feature_name=descriptor_name,
    )
