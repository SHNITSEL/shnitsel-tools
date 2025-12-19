from typing import Literal, Sequence
from shnitsel._contracts import needs
from shnitsel.core._api_info import API
import xarray as xr

from shnitsel.core.typedefs import AtXYZ
from shnitsel.filtering.structure_selection import StructureSelection
from shnitsel.geo.geocalc_.algebra import angle_, angle_cos_sin_
import numpy as np

from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _get_default_selection,
    _empty_descriptor_results,
)


@API()
@needs(dims={'atom'})
def angle(
    atXYZ: AtXYZ, a_index: int, b_index: int, c_index: int, *, deg: bool = False
) -> xr.DataArray:  # noqa: F821
    """Method to calculate the angle between atoms with indices `a_index`, `b_index`, and `c_index` in the positions DataArray throughout time.
    The `b_index` specifies the center atom at which the angle is located.
    The other two indices specify the legs of the angle.

    Can return results in radian (default) and degrees (if `deg=True`)

    Args:
        atXYZ (AtXYZ): DataArray with positions
        a_index (int): Index of first atom.
        b_index (int): Index of second center atom comprising the angle.
        c_index (int): Index of third atom.
        deg (bool, optional): Flag whether the results should be in degrees instead of radian. Defaults to False.

    Returns:
        xr.DataArray: The resulting angles between the denoted atoms.
    """
    a = atXYZ.isel(atom=a_index)
    b = atXYZ.isel(atom=b_index)
    c = atXYZ.isel(atom=c_index)
    ab = a - b
    cb = c - b
    result: xr.DataArray = angle_(ab, cb)
    if deg:
        result = result * 180 / np.pi
        result.attrs['units'] = 'degrees'
    else:
        result.attrs['units'] = 'rad'
    result.name = 'angle'
    result.attrs['long_name'] = r"$\theta_{%d,%d,%d}$" % (a_index, b_index, c_index)
    return result


@API()
@needs(dims={'atom'})
def angle_cos_sin(
    atXYZ: AtXYZ, a_index: int, b_index: int, c_index: int, *, deg: bool = False
) -> tuple[xr.DataArray, xr.DataArray]:
    """Method to calculate the cosine and sine of the angle between atoms with indices `a_index`, `b_index`, and `c_index` in the positions DataArray throughout time.
    The `b_index` specifies the center atom at which the angle is located.
    The other two indices specify the legs of the angle.

    Args:
        atXYZ (AtXYZ): DataArray with positions
        a_index (int): Index of first atom.
        b_index (int): Index of second center atom comprising the angle.
        c_index (int): Index of third atom.

    Returns:
        xr.DataArray: The resulting angles between the denoted atoms.
    """
    a = atXYZ.isel(atom=a_index)
    b = atXYZ.isel(atom=b_index)
    c = atXYZ.isel(atom=c_index)
    ab = a - b
    cb = c - b
    return angle_cos_sin_(ab, cb)


@API()
@needs(dims={'atom', 'direction'})
def get_angles(
    atXYZ: xr.DataArray,
    structure_selection: StructureSelection | None = None,
    deg: bool | Literal['trig'] = True,
    signed: bool = True,
) -> xr.DataArray:
    """Identify triples of bonded atoms (using RDKit) and calculate every bond angle for each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with at least dimensions
        `atom` and `direction`
    structure_selection, optional
        Object encapsulating feature selection on the structure whose positional information is provided in `atXYZ`.
        If this argument is omitted altogether, a default selection for all bonds within the structure is created.
    deg, optional
        Whether to return angles in degrees (as opposed to radians), by default True.
        Can also be set to the string literal `trig` if sin and cos of the calculated angle should be returned instead.
    signed, optional
        Whether the result should be returned with a sign or just as an absolute value in the range. Only relevant for `trig` option in `deg`.


    Returns
    -------
        An :py:class:`xarray.DataArray` of bond angles with dimensions `frame` and `angle`.

    Raises
    ------
    UserWarning
        If both `matches` and `mol` are specified.
    """
    structure_selection = _get_default_selection(
        structure_selection, atXYZ=atXYZ, default_levels=['angles']
    )

    angle_indices = list(structure_selection.angles_selected)

    if len(angle_indices) == 0:
        return _empty_descriptor_results(atXYZ)

    if isinstance(deg, bool):
        angle_arrs = [
            angle(atXYZ, a, b, c, deg=deg)
            .squeeze('atom', drop=True)
            .expand_dims('descriptor')
            for a, b, c in angle_indices
        ]

        angle_res = xr.concat(angle_arrs, dim='descriptor')

        descriptor_tex = [
            r"$\theta_{%d,%d,%d}$" % (a, b, c) for a, b, c in angle_indices
        ]
        descriptor_name = [r'angle(%d,%d,%d)' % (a, b, c) for a, b, c in angle_indices]
        descriptor_type = ['angle'] * len(descriptor_tex)

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
            *[[angle_cos_sin(atXYZ, a, b, c) for a, b, c in angle_indices]]
        )

        cos_res = [
            x.squeeze('atom', drop=True).expand_dims('descriptor') for x in cos_res
        ]
        sin_res = [
            x.squeeze('atom', drop=True).expand_dims('descriptor') for x in sin_res
        ]
        all_res: Sequence[xr.DataArray] = cos_res + sin_res

        if not signed:
            all_res = [np.abs(x) for x in all_res]

        descriptor_tex = [
            r"$\cos\theta_{%d,%d,%d}$" % (a, b, c) for a, b, c in angle_indices
        ] + [r"$\sin\theta_{%d,%d,%d}$" % (a, b, c) for a, b, c in angle_indices]

        descriptor_name = [
            r'cos(%d,%d,%d)' % (a, b, c) for a, b, c in angle_indices
        ] + [r'sin(%d,%d,%d)' % (a, b, c) for a, b, c in angle_indices]

        descriptor_type = ['cos'] * len(descriptor_tex) + ['sin'] * len(descriptor_tex)

        angle_res: xr.DataArray = xr.concat(all_res, dim='descriptor')  # type: ignore
        return _assign_descriptor_coords(
            angle_res,
            feature_descriptors=angle_indices,
            feature_type=descriptor_type,
            feature_tex_label=descriptor_tex,
            feature_name=descriptor_name,
        )
