from typing import Any, overload
from shnitsel._contracts import needs
from shnitsel.core._api_info import API
import xarray as xr

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import FeatureTypeLabel, StructureSelection
from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
    _get_default_selection,
)


@overload
def get_positions(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray],
    structure_selection: StructureSelection | None = None,
) -> TreeNode[Any, xr.DataArray]: ...
@overload
def get_positions(
    atXYZ_source: Trajectory | Frames | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection | None = None,
) -> xr.DataArray: ...


@API()
@needs(dims={'atom', 'direction'})
def get_positions(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray]
    | Trajectory
    | Frames
    | xr.Dataset
    | xr.DataArray,
    structure_selection: StructureSelection | None = None,
) -> TreeNode[Any, xr.DataArray] | xr.DataArray:
    """Return a descriptor-indexed set of positions for bats calculation.

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
        An :py:class:`xarray.DataArray` of positions with dimension `descriptor` to index the positions along.

    """

    if isinstance(atXYZ_source, TreeNode):
        return atXYZ_source.map_data(
            lambda x: get_positions(
                x,
                structure_selection=structure_selection,
            ),
            keep_empty_branches=True,
            dtype=xr.DataArray,
        )

    structure_selection = _get_default_selection(
        structure_selection, atXYZ_source=atXYZ_source, default_levels=['atoms']
    )

    atXYZ: xr.DataArray
    if isinstance(atXYZ_source, xr.DataArray):
        atXYZ = atXYZ_source
    else:
        atXYZ = atXYZ_source.atXYZ

    structure_selection = _get_default_selection(
        structure_selection, atXYZ_source=atXYZ, default_levels=['atoms']
    )

    position_indices = list(structure_selection.atoms_selected)

    if len(position_indices) == 0:
        return _empty_descriptor_results(atXYZ)

    positions_arrs = [atXYZ.sel(atom=a, drop=True) for a in position_indices]
    coordinates_x: list[xr.DataArray]
    coordinates_y: list[xr.DataArray]
    coordinates_z: list[xr.DataArray]

    coordinates_x, coordinates_y, coordinates_z = [
        list(a)
        for a in zip(
            *[
                (
                    arr.sel(direction='x', drop=True),
                    arr.sel(direction='y', drop=True),
                    arr.sel(direction='z', drop=True),
                )
                for arr in positions_arrs
            ]
        )
    ]

    coordinates_res = xr.concat(
        coordinates_x + coordinates_y + coordinates_z, dim='descriptor'
    )

    descriptor_tex = (
        [r'x_{%d}' % (a) for a in position_indices]
        + [r'y_{%d}' % (a) for a in position_indices]
        + [r'z_{%d}' % (a) for a in position_indices]
    )
    descriptor_name = (
        [r'pos_x(%d)' % (a) for a in position_indices]
        + [r'pos_y(%d)' % (a) for a in position_indices]
        + [r'pos_z(%d)' % (a) for a in position_indices]
    )
    descriptor_type: list[FeatureTypeLabel] = ['pos'] * len(descriptor_tex)

    return _assign_descriptor_coords(
        coordinates_res,
        feature_descriptors=position_indices * 3,
        feature_type=descriptor_type,
        feature_tex_label=descriptor_tex,
        feature_name=descriptor_name,
    )
