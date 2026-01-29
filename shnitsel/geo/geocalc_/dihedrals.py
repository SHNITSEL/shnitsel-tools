from typing import Any, Literal, Sequence, overload
from shnitsel._contracts import needs
from shnitsel.core._api_info import API, internal
import xarray as xr
import numpy as np

from shnitsel.core.typedefs import AtXYZ, DataArrayOrVar
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import (
    FeatureTypeLabel,
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc_.algebra import angle_, angle_cos_sin_, dcross, ddot, normal
from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
)
from shnitsel.filtering.helpers import _get_default_structure_selection


def _dihedral_deg(
    atXYZ: AtXYZ,
    a_index: int,
    b_index: int,
    c_index: int,
    d_index: int,
    full: bool = False,
) -> xr.DataArray:
    """Function to calculate the limited (0 to pi radian) dihedral angle between the points in arrays a,b,c and d.

    if `full=True`, calculates the signed/full dihedral angle (up to +-\\pi radian) between the points in arrays a,b,c and d.

    Parameters
    ----------
    a_index : int
        The first atom index
    b_index : int
        The second atom index
    c_index : int
        The third atom index
    d_index : int
        The fourth atom index
    full : bool, optional
        Flag to enforce calculation of the full dihedral in the range (-pi,pi)

    Returns
    -------
    xr.DataArray | xr.Variable
        The array of dihedral angels between the four input indices.
    """
    a = atXYZ.sel(atom=a_index, drop=True)
    b = atXYZ.sel(atom=b_index, drop=True)
    c = atXYZ.sel(atom=c_index, drop=True)
    d = atXYZ.sel(atom=d_index, drop=True)
    abc_normal = normal(a, b, c)
    bcd_normal = normal(b, c, d)
    if full:
        sign = np.sign(ddot(dcross(abc_normal, bcd_normal), (c - b)))
        res = angle_(abc_normal, bcd_normal) * sign
    else:
        res = angle_(abc_normal, bcd_normal)

    res.attrs['units'] = 'rad'
    return res


def _dihedral_trig_(
    atXYZ: AtXYZ,
    a_index: int,
    b_index: int,
    c_index: int,
    d_index: int,
    full: bool = False,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Function to calculate the sine and cosine of the dihedral between the points in arrays a,b,c and d.

    Parameters
    ----------
    a_index : int
        The first atom index
    b_index : int
        The second atom index
    c_index : int
        The third atom index
    d_index : int
        The fourth atom index

    Returns
    -------
    tuple[xr.DataArray, xr.DataArray]
        First the array of cosines and then the array of sines of the dihedral angle
    """
    a = atXYZ.sel(atom=a_index, drop=True)
    b = atXYZ.sel(atom=b_index, drop=True)
    c = atXYZ.sel(atom=c_index, drop=True)
    d = atXYZ.sel(atom=d_index, drop=True)
    abc_normal = normal(a, b, c)
    bcd_normal = normal(b, c, d)
    res = angle_cos_sin_(abc_normal, bcd_normal)
    res[0].attrs['units'] = 'trig'
    res[1].attrs['units'] = 'trig'
    return res


@overload
@needs(dims={'atom'})
def dihedral(
    atXYZ: AtXYZ,
    a_index: int,
    b_index: int,
    c_index: int,
    d_index: int,
    *,
    deg: Literal['trig'],
    full: bool = False,
) -> tuple[xr.DataArray, xr.DataArray]: ...


@overload
@needs(dims={'atom'})
def dihedral(
    atXYZ: AtXYZ,
    a_index: int,
    b_index: int,
    c_index: int,
    d_index: int,
    *,
    deg: bool = True,
    full: bool = False,
) -> xr.DataArray: ...


@API()
@needs(dims={'atom'})
def dihedral(
    atXYZ: AtXYZ,
    a_index: int,
    b_index: int,
    c_index: int,
    d_index: int,
    *,
    deg: bool | Literal['trig'] = True,
    full: bool = False,
) -> xr.DataArray | tuple[xr.DataArray, xr.DataArray]:
    """Calculate all dihedral angles between the atoms specified.
    The atoms specified need to be bonded in this sequence (a-b), (b-c), (c-d).

    Parameters
    ----------
    atXYZ : AtXYZ
        A ``DataArray`` of coordinates, with ``atom`` and ``direction`` dimensions
    a_index, b_index, c_index, d_index : int
        The four atom indices, where successive atoms should be bonded in this order.
    deg :  bool | Literal['trig'],optional
        Whether to return angles in degrees (True) or radians (False) or as cosine and sine ('trig'), by default False
    full : bool, optional
        Whether to return signed full dihedrals or unsigned (positive) dihedrals if False, by default False

    Returns
    -------
    xr.DataArray
        A ``DataArray`` containing dihedral angles (or the sin and cos thereof)
    """
    if isinstance(atXYZ, TreeNode):
        return atXYZ.map_data(
            dihedral,
            a_index=a_index,
            b_index=b_index,
            c_index=c_index,
            d_index=d_index,
            deg=deg,
            full=full,
        )
    if deg == 'trig':
        result_cos, result_sin = _dihedral_trig_(
            atXYZ, a_index, b_index, c_index, d_index, full=full
        )
        result_cos.name = 'cos(dihedral)'
        result_cos.attrs['long_name'] = r"\cos(\varphi_{%d,%d,%d,%d})" % (
            a_index,
            b_index,
            c_index,
            d_index,
        )
        result_sin.name = 'cos(dihedral)'
        result_sin.attrs['long_name'] = r"\cos(\varphi_{%d,%d,%d,%d})" % (
            a_index,
            b_index,
            c_index,
            d_index,
        )
        return result_cos, result_sin
    if isinstance(deg, bool):
        result: xr.DataArray = _dihedral_deg(
            atXYZ, a_index, b_index, c_index, d_index, full=full
        )
        if deg:
            result = result * 180 / np.pi
            result.attrs['units'] = 'degrees'
        else:
            result.attrs['units'] = 'rad'
        result.name = 'dihedral'
        result.attrs['long_name'] = r"\varphi_{%d,%d,%d,%d}" % (
            a_index,
            b_index,
            c_index,
            d_index,
        )
        return result


@overload
def get_dihedrals(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    deg: bool | Literal['trig'] = True,
    signed=True,
) -> TreeNode[Any, xr.DataArray]: ...


@overload
def get_dihedrals(
    atXYZ_source: Trajectory | Frames | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    deg: bool | Literal['trig'] = True,
    signed=True,
) -> xr.DataArray: ...


@API()
@needs(dims={'atom', 'direction'})
def get_dihedrals(
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
    """Identify quadruples of bonded atoms (using RDKit) and calculate the corresponding proper bond torsion for each
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
    deg: bool | Literal['trig'], optional
        Whether to return angles in degrees (as opposed to radians), by default True. Alternatively, return cos and sin (option `trig`) for each dihedral
    signed, optional
        Whether the result should be returned with a sign or just as an absolute value in the range. Triggers calculation of 'full' i.e. signed dihedrals.


    Returns
    -------
    TreeNode[Any, xr.DataArray] | xr.DataArray
        An :py:class:`xarray.DataArray` of bond torsions/dihedrals with dimension `descriptor` to index the dihedrals along.
    """

    if isinstance(atXYZ_source, TreeNode):
        return atXYZ_source.map_data(
            lambda x: get_dihedrals(
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
        default_levels=['dihedrals'],
        charge_info=charge_info,
    )

    dihedral_indices = list(structure_selection.dihedrals_selected)

    if len(dihedral_indices) == 0:
        return _empty_descriptor_results(position_data)

    dihedral_arrs = [
        dihedral(position_data, a, b, c, d, deg=deg, full=signed)
        for a, b, c, d in dihedral_indices
    ]

    if deg == 'trig':
        dih_angles_cos, dih_angles_sin = zip(*dihedral_arrs)
        descriptor_tex_cos = [
            r'\cos(\varphi_{%d,%d,%d,%d})' % (a, b, c, d)
            for (a, b, c, d) in dihedral_indices
        ]
        descriptor_tex_sin = [
            r'\sin(\varphi_{%d,%d,%d,%d})' % (a, b, c, d)
            for (a, b, c, d) in dihedral_indices
        ]
        descriptor_name_cos = [
            r'cos(%d,%d,%d,%d)' % (a, b, c, d) for (a, b, c, d) in dihedral_indices
        ]
        descriptor_name_sin = [
            r'sin(%d,%d,%d,%d)' % (a, b, c, d) for (a, b, c, d) in dihedral_indices
        ]
        descriptor_type_cos: list[FeatureTypeLabel] = ['cos_dih'] * len(
            descriptor_tex_cos
        )
        descriptor_type_sin: list[FeatureTypeLabel] = ['sin_dih'] * len(
            descriptor_tex_sin
        )

        dihedral_angles_extended: list[xr.DataArray] = [
            arr.expand_dims('descriptor') for arr in dih_angles_cos
        ] + [arr.expand_dims('descriptor') for arr in dih_angles_sin]

        dih_concatenated = xr.concat(dihedral_angles_extended, 'descriptor')

        dihedral_res = _assign_descriptor_coords(
            dih_concatenated,
            feature_name=descriptor_name_cos + descriptor_name_sin,
            feature_tex_label=descriptor_tex_cos + descriptor_tex_sin,
            feature_type=descriptor_type_cos + descriptor_type_sin,
            feature_descriptors=dihedral_indices + dihedral_indices,
        )

        dihedral_res: xr.DataArray = dihedral_res
        dihedral_res.name = "dihedrals"
        dihedral_res.attrs['units'] = 'trig'
        return dihedral_res
    else:
        dihedral_arrs_extended: list[xr.DataArray] = [
            arr.expand_dims('descriptor') for arr in dihedral_arrs
        ]
        dihedral_res = xr.concat(dihedral_arrs_extended, dim='descriptor')

        descriptor_tex = [
            r"\varphi_{%d,%d,%d,%d}" % (a, b, c, d) for a, b, c, d in dihedral_indices
        ]
        descriptor_name = [
            r'dih(%d,%d,%d,%d)' % (a, b, c, d) for a, b, c, d in dihedral_indices
        ]
        descriptor_type: list[FeatureTypeLabel] = ['dih'] * len(descriptor_tex)
        dihedral_res.name = "dihedrals"

        return _assign_descriptor_coords(
            dihedral_res,
            feature_descriptors=dihedral_indices,
            feature_type=descriptor_type,
            feature_tex_label=descriptor_tex,
            feature_name=descriptor_name,
        )
    # else:
    #     raise ValueError(
    #         "We only support boolean values for `deg` parameter in dihedral/torsion calculation."
    #     )

    # matches = _check_matches(matches_or_mol, atXYZ)['dihedrals']
    # if len(matches) == 0:
    #     return _empty_results(atXYZ)

    # _, atom_idxs, bond_idxs, bond_types, fragment_objs = zip(*matches)

    # assert all(len(x) == 4 for x in atom_idxs)
    # assert all(len(x) == 3 for x in bond_idxs)
    # assert all(len(x) == 3 for x in bond_types)

    # atom_positions = _positions(atXYZ, atom_idxs)
    # std_args = (atom_idxs, bond_idxs, bond_types, fragment_objs)

    # if ang:
    #     if signed:
    #         data = full_dihedral_(*atom_positions)
    #     else:
    #         data = dihedral_(*atom_positions)
    #     if ang == 'deg':
    #         data *= 180 / np.pi
    #     return _assign_descriptor_coords(
    #         data,
    #         *std_args,
    #         r"$\varphi_{%d,%d,%d,%d}$",
    #     )
    # else:
    #     if signed is not None:
    #         raise ValueError("Can't use `signed` parameter when ang==False")

    #     r0, r1, r2, r3 = atom_positions
    #     n012 = normal(r0, r1, r2)
    #     n123 = normal(r1, r2, r3)
    #     cos, sin = angle_cos_sin_(n012, n123)
    #     cos = _assign_descriptor_coords(cos, *std_args, r"$\cos\varphi_{%d,%d,%d,%d}$")
    #     sin = _assign_descriptor_coords(sin, *std_args, r"$\sin\varphi_{%d,%d,%d,%d}$")
    #     return xr.concat([cos, sin], dim='descriptor')
