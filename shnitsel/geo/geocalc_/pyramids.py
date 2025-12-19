import xarray as xr
import numpy as np

from shnitsel._contracts import needs
from shnitsel.core._api_info import API
from shnitsel.core.typedefs import AtXYZ
from shnitsel.filtering.structure_selection import StructureSelection
from shnitsel.geo.geocalc_.algebra import normal, angle_, normalize

from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
    _get_default_selection,
)


@needs(dims={'atom', 'direction'})
def pyramidalization_angle_(
    atXYZ: AtXYZ, x_index: int, a_index: int, b_index: int, c_index: int
) -> xr.DataArray:
    """Method to calculate the pyramidalization angle of a quadruple of atoms.

    The result will be $\\pi/2$ minus the angle between the normal of ABC and the vector BX.
    (I.e.: the pyramidalization at atom b)
    We choose the universal (independent of the choice of `b` within the tuple) way of calculating a unique
    pyramidalization angle accordingt to `https://doi.org/10.1063/5.0008368`, where the existence and derivation
    of such an angle is shown.

    Args:
        atXYZ (AtXYZ): Array with atom positions
        x_index (int): Index of the center atom in the pyramidalization
        a_index (int): Index of the first atom bonded to the `x`-atom
        b_index (int): Index of the second atom bonded to the `x`-atom
        c_index (int): Index of the third atom bonded to the `x`-atom

    Returns:
        xr.DataArray: The array of pyramidalization angles
    """
    # TODO: FIXME: Not sure if this is correct. This definition depends on the choice abc. Should b be the center molecule?

    # NOTE: According to https://doi.org/10.1063/5.0008368 this should yield a unique angle independent of the permutation of a,b,c if the distances are normalized first.
    # This should give the p-orbital-aligned perpendicular normal.

    x: xr.DataArray = atXYZ.sel(atom=x_index)
    a: xr.DataArray = atXYZ.sel(atom=a_index)
    b: xr.DataArray = atXYZ.sel(atom=b_index)
    c: xr.DataArray = atXYZ.sel(atom=c_index)

    da_norm = normalize(a - x)
    db_norm = normalize(b - x)
    dc_norm = normalize(c - x)
    orbital_aligned_normal = normal(da_norm, db_norm, dc_norm)

    return 0.5 * np.pi - angle_(orbital_aligned_normal, x - b)  # type: ignore # Cannot be a variable if provided with a DataArray


@API()
def get_pyramidalization(
    atXYZ: xr.DataArray,
    structure_selection: StructureSelection | None = None,
    deg: bool = False,
    signed=True,
) -> xr.DataArray:
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
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        ``frame``, ``atom`` and ``direction``
    structure_selection, optional
        An optional argument to specify the substructures for which pyramidalization angles should be calculated.
        If not provided, will be generated using `_get_default_selection()` using the atXYZ data for the pyramids level.
    deg, optional
        Whether to return angles in degrees (as opposed to radians), by default False
    signed, optional
        Whether the result should be returned with a sign or just as an absolute value. Defaults to True, yielding the signed pyramidalization.

    Returns
    -------
        An :py:class:`xarray.DataArray` of pyramidalizations with dimensions all dimensions but `atom` still intact and a new `descriptor` dimension introduced to index all the chosen quadruples for pyramidalization instead of the `atom` dimension.

    """
    structure_selection = _get_default_selection(
        structure_selection, atXYZ=atXYZ, default_levels=['pyramids']
    )

    pyramid_descriptors = list(structure_selection.pyramids_selected)

    if len(pyramid_descriptors) == 0:
        return _empty_descriptor_results(atXYZ)

    pyr_angles = [
        pyramidalization_angle_(atXYZ, x, a, b, c)
        for x, (a, b, c) in pyramid_descriptors
    ]

    descriptor_tex = [
        r'$\chi_{%d,%d}^{%d,%d}$' % (b, x, a, c) for x, (a, b, c) in pyramid_descriptors
    ]
    descriptor_name = [
        r'pyr(%d,(%d,%d,%d))' % (x, a, b, c) for x, (a, b, c) in pyramid_descriptors
    ]
    descriptor_type = ['pyr'] * len(descriptor_tex)

    pyr_angles_extended = [arr.expand_dims('descriptor') for arr in pyr_angles]

    pyr_concatenated = xr.concat(pyr_angles_extended, 'descriptor')

    pyr_res = _assign_descriptor_coords(
        pyr_concatenated,
        feature_name=descriptor_name,
        feature_tex_label=descriptor_tex,
        feature_type=descriptor_type,
        feature_descriptors=pyramid_descriptors,
    )

    pyr_res: xr.DataArray = pyr_res.set_xindex('descriptor')
    if deg:
        pyr_res *= 180 / np.pi
        pyr_res.attrs['units'] = 'degrees'
    else:
        pyr_res.attrs['units'] = 'rad'

    if not signed:
        pyr_res.attrs['sign'] = 'unsigned'
        pyr_res = np.abs(pyr_res)

    return pyr_res
