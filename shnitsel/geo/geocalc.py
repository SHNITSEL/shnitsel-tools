"""\
    This module will contain two types of closely-related functionality:
    - generic geometry functions
    - functions that use RDKit to identify internal coordinates, and then the above to calculate the values
    Currently, the first category of function is found under postprocess
"""

from typing import Literal, Iterable, Sequence


import numpy as np
import rdkit.Chem as rc
from rdkit.Chem import Mol
import xarray as xr

from shnitsel.core._api_info import API, internal
from shnitsel.filtering.structure_selection import (
    FeatureLevelOptions,
    FeatureLevelType,
    StructureSelection,
)
from ..analyze.generic import norm
from ..bridges import default_mol
from .._contracts import needs

from ..core.typedefs import AtXYZ

from shnitsel.geo.geomatch import flag_bats, flag_bla_chromophor


@needs(dims={'atom', 'direction'})
def get_bats(
    atXYZ: AtXYZ,
    structure_selection: StructureSelection | None = None,
    default_features: Sequence[FeatureLevelType] = ['bonds', 'angles', 'dihedrals'],
    signed: bool = False,
    deg: bool | Literal['trig'] = True,
) -> xr.DataArray:
    """Get bond lengths, angles and torsions/dihedrals.

    Parameters
    ----------
    atXYZ
        The positional data of atoms to use.
    structure_selection, optional: StructureSelection
        A feature selection to use. Can specify which features (positions, distances,
        angles, torsions or pyramidalizations) to include in the result.
        If not set, will be initialized to a default selection of all molecule-internal features
        as specified by the structure in the first frame of `atXYZ` and the features
        listed in `default_features`.
    default_features, optional: Sequence[FeatureLevelType]
        If no `structure_selection` object is provided, will select all features of these levels
        within the structure encoded in `atXYZ`.
        Options are
        - `atoms` for positional data,
        - `bonds` for distances between pairs of atoms (defaults to only bonds)
        - `angles` for angles between pairs of bonds between atoms.
        - `dihedrals` for torsion angles of bonds
        - `pyramids` for pyramidalization angles in the molecule.
        Defaults to using bonds, angles and dihedrals/torsions.
    signed, optional
        Whether to distinguish between clockwise and anticlockwise rotation,
        when returning angles as opposed to cosine & sine values;
        by default, do not distinguish.
        NB. This applies only to the dihedrals, not to the three-center angles.
        The latter are always unsigned.
    deg, optional: bool or Literal['trig']
        If True (the default), returns angles in degrees.
        If False, returns angles in radians.
        If set to 'trig' returns sines and cosines;

    Returns
    -------
        An :py:class:`xarray.DataArray` containing bond lengths, angles and tensions.

    Examples
    --------
        >>> import shnitsel as st
        >>> from shnitsel.core import geom
        >>> frames = st.open_frames('/tmp/A03_filtered.nc')
        >>> geom.get_bats(frames['atXYZ'])
    """
    # TODO: FIXME: Example is not up to date
    if matches_or_mol is None:
        mol = default_mol(atXYZ)
        matches_or_mol = flag_bats(mol)[0]

    d = {
        'bonds': get_bond_lengths(atXYZ, matches_or_mol=matches_or_mol),
        'angles': get_bond_angles(atXYZ, matches_or_mol=matches_or_mol, ang=ang),
        'dihedrals': get_bond_torsions(
            atXYZ, matches_or_mol=matches_or_mol, signed=signed, ang=ang
        ),
    }

    if pyr:
        d['pyr'] = get_pyramids(
            atXYZ, matches_or_mol=matches_or_mol, ang=ang, signed=signed
        )

    for k in d:
        d[k] = d[k].drop_vars(
            ['atom0', 'atom1', 'atom2', 'atom3', 'atNames', 'atNums'], errors='ignore'
        )

    return xr.concat(list(d.values()), dim='descriptor')


# @needs(dims={'atom', 'direction'})
# def identify_pyramids(mol: Mol) -> dict[int, list[int]]:
#     """Identify atoms with three bonds (using RDKit), specifically chain-internal atoms with
#     a single carbon, and chain-terminal atoms with two carbons.

#     Each 'pyramid' consists of four atoms. Three of these (the "plane" atoms) are consecutive, and the fourth (the "bending" atom)
#     is bonded to the middle atom of the plane atoms. The pyramidalization is the the angle between the plane of the plane atoms
#     and the bond from the middle plane atom to the bending atom.

#     Two sorts of pyramids are currently handled: terminal and chain-internal.

#     - Terminal pyramids are those where the central atom is bonded to two hydrogens and a single non-hydrogen;
#       for these, the central atom and the **hydrogens** constitute the plane and the non-hydrogen becomes the bending atom.
#     - Chain-internal pyramids are those where the central atom is bonded to non-hydrogens and a single hydrogen;
#       for these, the central atom and the **non-hydrogens** constitute the plane and the hydrogen becomes the bending atom.


#     Parameters
#     ----------
#     mol
#         An RDKit Mol object

#     Returns
#     -------
#     dict
#         {
#             x1: (a1, b1, c1),
#             x2: (a2, b2, c2),
#             ...
#         }
#     where

#         - a, b, c are indices of contiguous atoms used to determine the plane
#         - x is the index of an atom bonded to atom b, considered to be bending
#         above and beneath the plane
#     """
#     res = {}
#     for a in mol.GetAtoms():
#         bonds = a.GetBonds()
#         if len(bonds) != 3:
#             continue

#         current_idx = a.GetIdx()
#         hydrogens = []
#         non_hydrogens = []
#         for b in bonds:
#             a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
#             if a1.GetIdx() == current_idx:
#                 other = a2
#             else:
#                 assert a2.GetIdx() == current_idx
#                 other = a1
#             if other.GetAtomicNum() == 1:
#                 hydrogens.append(other.GetIdx())
#             else:
#                 non_hydrogens.append(other.GetIdx())
#         if len(hydrogens) == 2:  # Terminal double bond
#             assert len(non_hydrogens) == 1
#             plane_idxs = [hydrogens[0], current_idx, hydrogens[1]]
#             bending_idx = non_hydrogens[0]
#         else:  # Chain-internal double bond
#             assert len(hydrogens) == 1
#             assert len(non_hydrogens) == 2
#             plane_idxs = [non_hydrogens[0], current_idx, non_hydrogens[1]]
#             bending_idx = hydrogens[0]
#         res[bending_idx] = plane_idxs
#     return res


@API()
@needs(dims={'atom'})
def get_centered_geometry(atXYZ: AtXYZ, by_mass: Literal[False] = False) -> AtXYZ:
    """Helper function to set the center of the geometry (i.e. mean along the `atom` axis) to zero.

    Args:
        atXYZ (AtXYZ): Array of positional data
        by_mass (Literal[False], optional): Flag whether the centering/average should be center of mass or just plain average of positions. Defaults to False.

    Raises:
        NotImplementedError: Centering the COM instead of the mean is currently not implemented.

    Returns:
        AtXYZ: Resulting positions after centering.
    """
    if by_mass:
        raise NotImplementedError
    return atXYZ - atXYZ.mean('atom')


def rotational_procrustes_(
    A: np.ndarray, B: np.ndarray, weight: Iterable[float] | None = None
) -> np.ndarray:
    """Rotationally align the geometrie(s) in A to the single geometry in B.

    Parameters
    ----------
    A
        The geometries to process with shape
        ``(n_geometries, n_points, n_coordinates)``
    B
        The reference geometry with shape
        ``(n_points, n_coordinates)``
    weight, optional
        How much importance should be given to the alignment of
        each point, by default equal importance

    Returns
    -------
        An array with the same shape as A
    """
    from scipy.linalg import svd

    if weight is not None:
        A = np.diag(weight) @ A

    # np.matrix_transpose always swaps last two axes, whereas
    # NDArray.T reverses the order of all axes.
    t = np.matrix_transpose
    # The following uses a double transpose in imitation of
    # scipy's orthogonal_procrustes, where this is said to
    # save memory. t(t(B) @ A) == t(A) @ B.
    u, _, vt = svd(t(t(B) @ A))
    # Flip the sign of the last row of each stacked vt matrix
    # depending on the sign of the corresponding determinant.
    # This is an alternative implementation of the algorithm
    # used in qcdev's procrustes.rotation.
    vt[..., -1, :] *= np.sign(np.linalg.det(u @ vt))[:, None]
    R = u @ vt
    return A @ R


def rotational_procrustes(
    A: xr.DataArray,
    B: xr.DataArray,
    dim0: str = 'atom',
    dim1: str = 'direction',
    weight: Iterable[float] | None = None,
) -> xr.DataArray:
    """Rotationally align the geometrie(s) in A to the single geometry in B.

    Parameters
    ----------
    A
        The geometries to process
    B
        The reference geometry
    dim0, optional
        The name of the dimension over points to be rotated;
        must be present in ``A`` and ``B`; by default 'atom'
    dim1, optional
        The name of the dimension over the coordinates of the aforementioned
        points; must be present in ``A`` and ``B`; by default 'direction'
    weight, optional
        How much importance should be given to the alignment of
        each point (atom), by default equal importance

    Returns
    -------
        An xr.DataArray with the same shape as ``A``
    """
    return xr.apply_ufunc(
        rotational_procrustes_,
        A,
        B,
        input_core_dims=[[dim0, dim1], [dim0, dim1]],
        output_core_dims=[[dim0, dim1]],
        kwargs={'weight': weight},
    )


@needs(dims={'atom', 'direction'})
def kabsch(
    atXYZ: xr.DataArray,
    reference_or_indexers: xr.DataArray | dict | None = None,
    **indexers_kwargs,
) -> xr.DataArray:
    """Rotationally align the molecular geometries in ``atXYZ`` to a single molecular geometry

    Parameters
    ----------
    atXYZ
        The geometries to process (with dims 'atom', 'direction')

    reference_or_indexers, optional
        Either a reference geometry (with dims 'atom', 'direction')
        or an indexer dictionary which will be passed to ``atXYZ.sel()``

    **indexer_kwargs
        The keyword-argument form of the indexer to be passed to ``atXYZ.sel()``

    Returns
    -------
        The aligned geometries

    Raises
    ------
    ValueError
        If nothing is done to indicate a reference geometry, i.e.
        neither reference_or_indexers nor indexer_kwargs are passed
    """

    if isinstance(reference_or_indexers, xr.DataArray):
        reference = reference_or_indexers
    elif isinstance(reference_or_indexers, dict):
        reference = atXYZ.sel(reference_or_indexers)
    elif len(indexers_kwargs) != 0:
        reference = atXYZ.sel(indexers_kwargs)
    elif 'frame' in atXYZ.dims:
        reference = atXYZ.isel(frame=0)
    else:
        raise ValueError("Please specify a reference geometry")

    # TODO: is it ever necessary to center the molecule?
    # If so, should this always be done using the physical center-of-mass,
    # or is it ever appropriate to use the unweighted mean of points?

    return rotational_procrustes(atXYZ, reference)
