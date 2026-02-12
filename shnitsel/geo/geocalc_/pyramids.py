from typing import Any, Literal, overload
import xarray as xr
import numpy as np

from shnitsel._contracts import needs
from shnitsel.core._api_info import API
from shnitsel.core.typedefs import AtXYZ
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import (
    FeatureTypeLabel,
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc_.algebra import angle_cos_sin_, normal, angle_, normalize

from shnitsel.geo.geocalc_.helpers import (
    AngleOptions,
    _assign_descriptor_coords,
    _at_XYZ_subset_to_descriptor,
    _empty_descriptor_results,
)
from shnitsel.filtering.helpers import _get_default_structure_selection


@overload
@needs(dims={'atom', 'direction'})
def pyramidalization_angle(
    atXYZ: AtXYZ,
    x_index: int | list[int],
    a_index: int | list[int],
    b_index: int | list[int],
    c_index: int | list[int],
    angles: Literal['trig'],
) -> tuple[xr.DataArray, xr.DataArray]: ...


@overload
@needs(dims={'atom', 'direction'})
def pyramidalization_angle(
    atXYZ: AtXYZ,
    x_index: int | list[int],
    a_index: int | list[int],
    b_index: int | list[int],
    c_index: int | list[int],
    angles: AngleOptions = 'deg',
) -> xr.DataArray: ...


@needs(dims={'atom', 'direction'})
def pyramidalization_angle(
    atXYZ: AtXYZ,
    x_index: int | list[int],
    a_index: int | list[int],
    b_index: int | list[int],
    c_index: int | list[int],
    angles: AngleOptions = 'deg',
) -> xr.DataArray | tuple[xr.DataArray, xr.DataArray]:
    """Method to calculate the pyramidalization angle of a quadruple of atoms.

    The result will be $\\pi/2$ minus the angle between the normal of ABC and the vector BX.
    (I.e.: the pyramidalization at atom b)
    We choose the universal (independent of the choice of `b` within the tuple) way of calculating a unique
    pyramidalization angle accordingt to `https://doi.org/10.1063/5.0008368`, where the existence and derivation
    of such an angle is shown.

    Parameters
    ----------
    atXYZ : AtXYZ
        Array with atom positions
    x_index : int | list[int]
        Index of the center atom in the pyramidalization
    a_index : int | list[int]
        Index of the first atom bonded to the `x`-atom
    b_index : int | list[int]
        Index of the second atom bonded to the `x`-atom
    c_index : int | list[int]
        Index of the third atom bonded to the `x`-atom
    angles : Literal['deg', 'rad', 'trig'], default='deg'
        Option parameter to control the unit/representation of the angle.
        Use 'deg' for results in 'degrees', 'rad' for results in 'radian'
        and 'trig' for a representation as a sin and a cos

    Returns
    -------
    xr.DataArray
        The array of pyramidalization angles
    tuple[xr.DataArray, xr.DataArray]
        If `deg='trig'`, return the pair of xr.DataArray sets of cosine and then sines of the pyramidalization angle.

    Notes
    -----
    According to https://doi.org/10.1063/5.0008368 this should yield a unique angle independent of the permutation of a,b,c if the distances are normalized first.
    This should give the p-orbital-aligned perpendicular normal.
    """
    if isinstance(atXYZ, TreeNode):
        return atXYZ.map_data(
            pyramidalization_angle,
            x_index=x_index,
            a_index=a_index,
            b_index=b_index,
            c_index=c_index,
            angles=angles,
        )
    # NOTE: According to https://doi.org/10.1063/5.0008368 this should yield a unique angle independent of the permutation of a,b,c if the distances are normalized first.
    # This should give the p-orbital-aligned perpendicular normal.

    x: xr.DataArray = _at_XYZ_subset_to_descriptor(atXYZ.sel(atom=x_index))
    a: xr.DataArray = _at_XYZ_subset_to_descriptor(atXYZ.sel(atom=a_index))
    b: xr.DataArray = _at_XYZ_subset_to_descriptor(atXYZ.sel(atom=b_index))
    c: xr.DataArray = _at_XYZ_subset_to_descriptor(atXYZ.sel(atom=c_index))

    da_norm = normalize(a - x)
    db_norm = normalize(b - x)
    dc_norm = normalize(c - x)
    orbital_aligned_normal = normal(da_norm, db_norm, dc_norm)

    if angles == 'trig':
        cos_raw, sin_raw = angle_cos_sin_(orbital_aligned_normal, x - b)

        if isinstance(a, int):
            cos_raw.attrs['long_name'] = r"\cos(\chi_{%d}^{%d,%d,%d})" % (x, a, b, c)
            sin_raw.attrs['long_name'] = r"\sin(\chi_{%d}^{%d,%d,%d})" % (x, a, b, c)
        else:
            cos_raw.attrs['long_name'] = r"\cos(\chi_{x}^{a,b,c})"
            sin_raw.attrs['long_name'] = r"\sin(\chi_{%d}^{%d,%d,%d})" % (x, a, b, c)
        # The 90-x swaps cos and sin
        return (sin_raw, cos_raw)  # type: ignore # Cannot be a variable if provided with a DataArray

    angle_rad = 0.5 * np.pi - angle_(orbital_aligned_normal, x - b)

    if angles == 'deg':
        angle_rad *= 180 / np.pi
        angle_rad.attrs['units'] = 'degrees'
    else:
        angle_rad.attrs['units'] = 'rad'

    if isinstance(a, int):
        angle_rad.attrs['long_name'] = r"\chi_{%d}^{%d,%d,%d}" % (x, a, b, c)
    else:
        angle_rad.attrs['long_name'] = r"\chi_{x}^{a,b,c}"

    return angle_rad  # type: ignore # Cannot be a variable if provided with a DataArray


@overload
def get_pyramidalization(
    atXYZ_source: TreeNode[Any, ShnitselDataset | xr.Dataset | xr.DataArray],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    angles: AngleOptions = 'deg',
    signed=True,
) -> TreeNode[Any, xr.DataArray]: ...


@overload
def get_pyramidalization(
    atXYZ_source: ShnitselDataset | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    angles: AngleOptions = 'deg',
    signed=True,
) -> xr.DataArray: ...


@API()
def get_pyramidalization(
    atXYZ_source: TreeNode[Any, ShnitselDataset | xr.Dataset | xr.DataArray]
    | ShnitselDataset
    | xr.Dataset
    | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    angles: AngleOptions = 'deg',
    signed: bool = True,
) -> TreeNode[Any, xr.DataArray] | xr.DataArray:
    """Identify atoms with three bonds (using RDKit) and calculate the corresponding pyramidalization angles
    for each frame.

    Each 'pyramid' consists of four atoms. Three of these (the "plane" atoms) are consecutive, and the fourth (the "bending" atom)
    is bonded to the middle atom of the plane atoms. The pyramidalization is the the angle between the plane of the plane atoms
    and the bond from the middle plane atom to the bending atom.

    Two sorts of pyramids are currently handled: terminal and chain-internal.

    - Terminal pyramids are those where the central atom is bonded to two hydrogens and a single non-hydrogen;
      for these, the central atom and the **hydrogens** constitute the plane and the non-hydrogen becomes the bending atom.
    - Chain-internal pyramids are those where the central atom is bonded to non-hydrogens and a single hydrogen;
      for these, the central atom and the **non-hydrogens** constitute the plane and the hydrogen becomes the bending atom.

    Parameters
    ----------
    atXYZ_source : TreeNode[Any, ShnitselDataset| xr.Dataset | xr.DataArray] | ShnitselDataset| xr.Dataset | xr.DataArray
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions ``atom`` and ``direction`` or another source of positional data like a trajectory, a frameset, a dataset representing either of those or a tree structure holding such data.
    structure_selection : StructureSelection | StructureSelectionDescriptor, optional
        An optional argument to specify the substructures for which pyramidalization angles should be calculated.
        If not provided, will be generated using `_get_default_structure_selection()` using the atXYZ data for the pyramids level.
    deg : bool | Literal['trig'] = True, optional
        Whether to return angles in degrees (as opposed to radians), by default False.
    angles : Literal['deg', 'rad', 'trig'], default='deg'
        Option parameter to control the unit/representation of the pyramidalization angle.
        Use 'deg' for results in 'degrees', 'rad' for results in 'radian'.
        Alternatively, with the option `trig`, this will yield the sin and cos of each pyramidalization angle instead.
    signed : bool, optional
        Whether the result should be returned with a sign or just as an absolute value. Defaults to True, yielding the signed pyramidalization.

    Returns
    -------
    TreeNode[Any, xr.DataArray] | xr.DataArray
        An :py:class:`xarray.DataArray` of pyramidalizations with dimensions all dimensions but `atom` still intact and a new `descriptor` dimension introduced to index all the chosen quadruples for pyramidalization instead of the `atom` dimension.

    """

    if isinstance(atXYZ_source, TreeNode):
        return atXYZ_source.map_data(
            lambda x: get_pyramidalization(
                x, structure_selection=structure_selection, angles=angles, signed=signed
            ),
            keep_empty_branches=True,
            dtype=xr.DataArray,
        )

    position_data: xr.DataArray
    position_source: ShnitselDataset | xr.Dataset | xr.DataArray
    charge_info: int | None
    if isinstance(atXYZ_source, xr.DataArray):
        position_data = atXYZ_source
        position_source = atXYZ_source
        charge_info = None
    else:
        wrapped_ds = wrap_dataset(atXYZ_source, ShnitselDataset)
        position_data = wrapped_ds.atXYZ
        position_source = wrapped_ds
        charge_info = int(wrapped_ds.charge)

    structure_selection = _get_default_structure_selection(
        structure_selection,
        atXYZ_source=position_source,
        default_levels=['pyramids'],
        charge_info=charge_info,
    )

    pyramid_descriptors = list(structure_selection.pyramids_selected)

    if len(pyramid_descriptors) == 0:
        return _empty_descriptor_results(position_data)

    x_indices, a_indices, b_indices, c_indices = [
        list(x) for x in zip(*[[x, a, b, c] for x, (a, b, c) in pyramid_descriptors])
    ]

    pyr_angles = pyramidalization_angle(
        position_data, x_indices, a_indices, b_indices, c_indices, angles=angles
    )

    if angles == 'trig':
        pyr_angles_cos, pyr_angles_sin = pyr_angles
        descriptor_tex_cos = [
            r'\cos(\chi_{%d}^{%d,%d,%d})' % (x, a, b, c)
            for x, (a, b, c) in pyramid_descriptors
        ]
        descriptor_tex_sin = [
            r'\sin(\chi_{%d}^{%d,%d,%d})' % (x, a, b, c)
            for x, (a, b, c) in pyramid_descriptors
        ]
        descriptor_name_cos = [
            r'cos(%d,(%d,%d,%d))' % (x, a, b, c) for x, (a, b, c) in pyramid_descriptors
        ]
        descriptor_name_sin = [
            r'sin(%d,(%d,%d,%d))' % (x, a, b, c) for x, (a, b, c) in pyramid_descriptors
        ]
        descriptor_type_cos: list[FeatureTypeLabel] = ['cos_pyr'] * len(
            descriptor_tex_cos
        )
        descriptor_type_sin: list[FeatureTypeLabel] = ['sin_pyr'] * len(
            descriptor_tex_sin
        )

        pyr_angles_extended: list[xr.DataArray] = [pyr_angles_cos, pyr_angles_sin]

        pyr_concatenated = xr.concat(pyr_angles_extended, 'descriptor')

        pyr_res = _assign_descriptor_coords(
            pyr_concatenated,
            feature_name=descriptor_name_cos + descriptor_name_sin,
            feature_tex_label=descriptor_tex_cos + descriptor_tex_sin,
            feature_type=descriptor_type_cos + descriptor_type_sin,
            feature_descriptors=pyramid_descriptors + pyramid_descriptors,
        )

        pyr_res: xr.DataArray = pyr_res
        pyr_res.name = "pyramids"
        pyr_res.attrs['units'] = 'trig'
        pyr_res.attrs['unitdim'] = 'angles'
        return pyr_res

    else:
        descriptor_tex = [
            r'\chi_{%d}^{%d,%d,%d}' % (x, a, b, c)
            for x, (a, b, c) in pyramid_descriptors
        ]
        descriptor_name = [
            r'pyr(%d,(%d,%d,%d))' % (x, a, b, c) for x, (a, b, c) in pyramid_descriptors
        ]
        descriptor_type: list[FeatureTypeLabel] = ['pyr'] * len(descriptor_tex)

        pyr_res = _assign_descriptor_coords(
            pyr_angles,
            feature_name=descriptor_name,
            feature_tex_label=descriptor_tex,
            feature_type=descriptor_type,
            feature_descriptors=pyramid_descriptors,
        )
        pyr_res.name = "pyramids"

        pyr_res: xr.DataArray = pyr_res
        if angles == 'deg':
            pyr_res.attrs['units'] = 'degrees'
        else:
            pyr_res.attrs['units'] = 'rad'
        pyr_res.attrs['unitdim'] = 'angles'

        if not signed:
            pyr_res.attrs['sign'] = 'unsigned'
            pyr_res = np.abs(pyr_res)

        return pyr_res
