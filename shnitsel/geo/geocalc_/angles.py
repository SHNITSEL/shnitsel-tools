from typing import Any, Literal, Sequence, overload
from shnitsel._contracts import needs
from shnitsel.core._api_info import API
import xarray as xr

from shnitsel.core.typedefs import AtXYZ
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import (
    FeatureTypeLabel,
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc_.algebra import angle_, angle_cos_sin_
import numpy as np

from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
)
from shnitsel.filtering.helpers import _get_default_structure_selection


@API()
@needs(dims={'atom'})
def angle(
    atXYZ: AtXYZ, a_index: int, b_index: int, c_index: int, *, deg: bool = False
) -> xr.DataArray:  # noqa: F821
    """Method to calculate the angle between atoms with indices `a_index`, `b_index`, and `c_index` in the positions DataArray throughout time.
    The `b_index` specifies the center atom at which the angle is located.
    The other two indices specify the legs of the angle.

    Can return results in radian (default) and degrees (if `deg=True`)

    Parameters
    ----------
    atXYZ : AtXYZ
        DataArray with positions
    a_index : int
        Index of first atom.
    b_index : int
        Index of second center atom comprising the angle.
    c_index : int
        Index of third atom.
    deg : bool, optional
        Flag whether the results should be in degrees instead of radian. Defaults to False.

    Returns
    -------
    xr.DataArray
        The resulting angles between the denoted atoms.
    """
    if isinstance(atXYZ, TreeNode):
        return atXYZ.map_data(
            angle,
            a_index=a_index,
            b_index=b_index,
            c_index=c_index,
            deg=deg,
        )
    a = atXYZ.sel(atom=a_index, drop=True)
    b = atXYZ.sel(atom=b_index, drop=True)
    c = atXYZ.sel(atom=c_index, drop=True)
    ab = a - b
    cb = c - b
    result: xr.DataArray = angle_(ab, cb)
    if deg:
        result = result * 180 / np.pi
        result.attrs['units'] = 'degrees'
    else:
        result.attrs['units'] = 'rad'
    result.name = 'angle'
    result.attrs['long_name'] = r"\theta_{%d,%d,%d}" % (a_index, b_index, c_index)
    return result


@API()
@needs(dims={'atom'})
def angle_cos_sin(
    atXYZ: AtXYZ, a_index: int, b_index: int, c_index: int, *, deg: bool = False
) -> tuple[xr.DataArray, xr.DataArray]:
    """Method to calculate the cosine and sine of the angle between atoms with indices `a_index`, `b_index`, and `c_index` in the positions DataArray throughout time.
    The `b_index` specifies the center atom at which the angle is located.
    The other two indices specify the legs of the angle.

    Parameters
    ----------
    atXYZ : AtXYZ
        DataArray with positions
    a_index : int
        Index of first atom.
    b_index : int
        Index of second center atom comprising the angle.
    c_index : int
        Index of third atom.

    Returns
    -------
    xr.DataArray
        The resulting angles between the denoted atoms.
    """
    if isinstance(atXYZ, TreeNode):
        return atXYZ.map_data(
            angle_cos_sin,
            a_index=a_index,
            b_index=b_index,
            c_index=c_index,
            deg=deg,
        )
    a = atXYZ.sel(atom=a_index, drop=True)
    b = atXYZ.sel(atom=b_index, drop=True)
    c = atXYZ.sel(atom=c_index, drop=True)
    ab = a - b
    cb = c - b

    res_cos, res_sin = angle_cos_sin_(ab, cb)
    res_cos.name = 'cos'
    res_cos.attrs['units'] = 'trig'
    res_sin.name = 'sin'
    res_sin.attrs['units'] = 'trig'
    return res_cos, res_sin


@overload
def get_angles(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    deg: bool | Literal['trig'] = True,
    signed=True,
) -> TreeNode[Any, xr.DataArray]: ...


@overload
def get_angles(
    atXYZ_source: Trajectory | Frames | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    deg: bool | Literal['trig'] = True,
    signed=True,
) -> xr.DataArray: ...


@API()
@needs(dims={'atom', 'direction'})
def get_angles(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray]
    | Trajectory
    | Frames
    | xr.Dataset
    | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    deg: bool | Literal['trig'] = True,
    signed: bool = True,
) -> TreeNode[Any, xr.DataArray] | xr.DataArray:
    """Identify triples of bonded atoms (using RDKit) and calculate bond angles for each frame.

    Parameters
    ----------
    atXYZ_source : xr.DataArray | TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray] | Trajectory | Frames | xr.Dataset
         An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions ``atom`` and
         ``direction`` or another source of positional data like a trajectory, a frameset,
         a dataset representing either of those or a tree structure holding such data.
    structure_selection: StructureSelection | StructureSelectionDescriptor, optional
         Object encapsulating feature selection on the structure whose positional information is provided in `atXYZ`.
         If this argument is omitted altogether, a default selection for all bonds within the structure is created.
    deg: bool | Literal['trig'], optional
         Whether to return angles in degrees (as opposed to radians), by default True.
         Can also be set to the string literal `trig` if sin and cos of the calculated angle should be returned instead.
    signed: bool, optional
         Whether the result should be returned with a sign or just as an absolute value in the range. Only relevant for `trig` option in `deg`.


    Returns
    -------
    TreeNode[Any, xr.DataArray] | xr.DataArray
        An :py:class:`xarray.DataArray` of bond angles with dimension `descriptor` to index the angles along.

    """

    if isinstance(atXYZ_source, TreeNode):
        return atXYZ_source.map_data(
            lambda x: get_angles(
                x, structure_selection=structure_selection, deg=deg, signed=signed
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
        default_levels=['angles'],
        charge_info=charge_info,
    )

    angle_indices = list(structure_selection.angles_selected)

    if len(angle_indices) == 0:
        return _empty_descriptor_results(position_data)

    if isinstance(deg, bool):
        angle_arrs = [
            angle(position_data, a, b, c, deg=deg).expand_dims('descriptor')
            for a, b, c in angle_indices
        ]

        angle_res = xr.concat(angle_arrs, dim='descriptor')
        angle_res.name = "angles"

        descriptor_tex = [r"\theta_{%d,%d,%d}" % (a, b, c) for a, b, c in angle_indices]
        descriptor_name = [r'angle(%d,%d,%d)' % (a, b, c) for a, b, c in angle_indices]
        descriptor_type: list[FeatureTypeLabel] = ['angle'] * len(descriptor_tex)

        return _assign_descriptor_coords(
            angle_res,
            feature_descriptors=angle_indices,
            feature_type=descriptor_type,
            feature_tex_label=descriptor_tex,
            feature_name=descriptor_name,
        )
    else:
        # Trigonometric results requested
        cos_res: Sequence[xr.DataArray]
        sin_res: Sequence[xr.DataArray]
        cos_res, sin_res = zip(
            *[angle_cos_sin(position_data, a, b, c) for a, b, c in angle_indices]
        )

        cos_res = [
            x.expand_dims('descriptor')  # .squeeze('atom', drop=True,)
            for x in cos_res
        ]
        sin_res = [
            x.expand_dims('descriptor')  # .squeeze('atom', drop=True)
            for x in sin_res
        ]
        all_res: Sequence[xr.DataArray] = cos_res + sin_res

        if not signed:
            all_res = [np.abs(x) for x in all_res]

        descriptor_tex = [
            r"\cos\theta_{%d,%d,%d}" % (a, b, c) for a, b, c in angle_indices
        ] + [r"\sin\theta_{%d,%d,%d}" % (a, b, c) for a, b, c in angle_indices]

        descriptor_name = [
            r'cos(%d,%d,%d)' % (a, b, c) for a, b, c in angle_indices
        ] + [r'sin(%d,%d,%d)' % (a, b, c) for a, b, c in angle_indices]

        descriptor_type: list[FeatureTypeLabel] = ['cos_angle'] * len(cos_res) + [
            'sin_angle'
        ] * len(sin_res)  # pyright: ignore[reportAssignmentType] # Is allowed string, but not generally advertised

        angle_res: xr.DataArray = xr.concat(all_res, dim='descriptor')  # type: ignore
        angle_res.name = "angles"
        return _assign_descriptor_coords(
            angle_res,
            feature_descriptors=angle_indices + angle_indices,
            feature_type=descriptor_type,
            feature_tex_label=descriptor_tex,
            feature_name=descriptor_name,
        )
