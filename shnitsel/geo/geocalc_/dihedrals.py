from typing import Literal, Sequence
from shnitsel._contracts import needs
from shnitsel.core._api_info import API, internal
import xarray as xr
import numpy as np

from shnitsel.core.typedefs import AtXYZ, DataArrayOrVar
from shnitsel.filtering.structure_selection import FeatureTypeLabel, StructureSelection
from shnitsel.geo.geocalc_.algebra import angle_, dcross, ddot, normal
from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
    _get_default_selection,
)


def dihedral_(
    atXYZ: AtXYZ,
    a_index: int,
    b_index: int,
    c_index: int,
    d_index: int,
    full: bool = False,
) -> xr.DataArray:
    """Function to calculate the limited (0 to pi radian) dihedral angle between the points in arrays a,b,c and d.

    if `full=True`, calculates the signed/full dihedral angle (up to +-\\pi radian) between the points in arrays a,b,c and d.

    Args:
        a_index (int): The first atom index
        b_index (int): The second atom index
        c_index (int): The third atom index
        d_index (int): The fourth atom index
        full (bool, optional): Flag to enforce calculation of the full dihedral in the range (-pi,pi)

    Returns:
        xr.DataArray | xr.Variable: The array of dihedral angels between the four input indices.
    """
    a = atXYZ.sel(atom=a_index, drop=True)
    b = atXYZ.sel(atom=b_index, drop=True)
    c = atXYZ.sel(atom=c_index, drop=True)
    d = atXYZ.sel(atom=d_index, drop=True)
    abc_normal = normal(a, b, c)
    bcd_normal = normal(b, c, d)
    if full:
        sign = np.sign(ddot(dcross(abc_normal, bcd_normal), (c - b)))
        return angle_(abc_normal, bcd_normal) * sign
    else:
        return angle_(abc_normal, bcd_normal)


@API()
@needs(dims={'atom'})
def dihedral(
    atXYZ: AtXYZ,
    a_index: int,
    b_index: int,
    c_index: int,
    d_index: int,
    *,
    deg: bool = False,
    full: bool = False,
) -> xr.DataArray:
    """Calculate all dihedral angles between the atoms specified.
    The atoms specified need to be bonded in this sequence (a-b), (b-c), (c-d).

    Parameters
    ----------
    atXYZ
        A ``DataArray`` of coordinates, with ``atom`` and ``direction`` dimensions
    a_index, b_index, c_index, d_index
        The four atom indices, where successive atoms should be bonded in this order.
    deg
        Whether to return angles in degrees (True) or radians (False), by default False
    full
        Whether to return signed full dihedrals or unsigned (positive) dihedrals if False, by default False

    Returns
    -------
        A ``DataArray`` containing dihedral angles
    """
    result: xr.DataArray = dihedral_(
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


@API()
@needs(dims={'atom', 'direction'})
def get_dihedrals(
    atXYZ: xr.DataArray,
    structure_selection: StructureSelection | None = None,
    deg: bool = True,
    signed: bool = True,
) -> xr.DataArray:
    """Identify quadruples of bonded atoms (using RDKit) and calculate the corresponding proper bond torsion for each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        at least `atom` and `direction`
    structure_selection, optional
        Object encapsulating feature selection on the structure whose positional information is provided in `atXYZ`.
        If this argument is omitted altogether, a default selection for all bonds within the structure is created.

    deg, optional
        Whether to return angles in degrees (as opposed to radians), by default True.
    signed, optional
        Whether the result should be returned with a sign or just as an absolute value in the range. Triggers calculation of 'full' i.e. signed dihedrals.


    Returns
    -------
        An :py:class:`xarray.DataArray` of bond torsions with dimension `descriptor` to index the dihedrals along.
    """

    structure_selection = _get_default_selection(
        structure_selection, atXYZ_source=atXYZ, default_levels=['dihedrals']
    )

    dihedral_indices = list(structure_selection.dihedrals_selected)

    if len(dihedral_indices) == 0:
        return _empty_descriptor_results(atXYZ)

    if isinstance(deg, bool):
        dihedral_arrs = [
            dihedral(atXYZ, a, b, c, d, deg=deg, full=signed)
            .expand_dims('descriptor')
            for a, b, c, d in dihedral_indices
        ]

        dihedral_res = xr.concat(dihedral_arrs, dim='descriptor')

        descriptor_tex = [
            r"\varphi_{%d,%d,%d,%d}" % (a, b, c, d) for a, b, c, d in dihedral_indices
        ]
        descriptor_name = [
            r'dih(%d,%d,%d,%d)' % (a, b, c, d) for a, b, c, d in dihedral_indices
        ]
        descriptor_type: list[FeatureTypeLabel] = ['dih'] * len(descriptor_tex)

        return _assign_descriptor_coords(
            dihedral_res,
            feature_descriptors=dihedral_indices,
            feature_type=descriptor_type,
            feature_tex_label=descriptor_tex,
            feature_name=descriptor_name,
        )
    else:
        raise ValueError(
            "We only support boolean values for `deg` parameter in dihedral/torsion calculation."
        )

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
