"""\
    This module will contain two types of closely-related functionality:
    - generic geometry functions
    - functions that use RDKit to identify internal coordinates, and then the above to calculate the values
    Currently, the first category of function is found under postprocess
"""

from collections.abc import Collection
from itertools import combinations, product
from logging import warning
from typing import Literal, TypeAlias


import numpy as np
import rdkit.Chem as rc
from rdkit.Chem import Mol
import xarray as xr

from shnitsel.io.helpers import get_symbol_from_atom_number
from .generic import subtract_combinations, norm
from ..bridges import default_mol
from .xrhelpers import expand_midx
from .._contracts import needs

from shnitsel.core.geo.geomatch import flag_bats

AtXYZ: TypeAlias = xr.DataArray


def dnorm(a):
    return norm(a, dim='direction')


def dcross(a, b):
    return xr.cross(a, b, dim='direction')


def ddot(a, b):
    return xr.dot(a, b, dim='direction')


def angle_(a, b):
    return np.arccos(ddot(a, b) / (dnorm(a) * dnorm(b)))

def angle_cos_sin_(a, b):
    """
    returns the cosine and sine of the angle between two vectors
    """
    prod = dnorm(a) * dnorm(b)
    return (
        ddot(a, b) / prod,
        dnorm(dcross(a, b)) / prod,
    )


def normal(a, b, c):
    return dcross(a - b, c - b)


def dihedral_(a, b, c, d):
    abc = normal(a, b, c)
    bcd = normal(b, c, d)
    return angle_(abc, bcd)


def full_dihedral_(a, b, c, d):
    abc = normal(a, b, c)
    bcd = normal(b, c, d)
    sign = np.sign(ddot(dcross(abc, bcd), (c - b)))
    return sign * angle_(abc, bcd)


@needs(dims={'atom'})
def dihedral(
    atXYZ: AtXYZ,
    i: int,
    j: int,
    k: int,
    l: int,
    *,
    deg: bool = False,
    full: bool = False,
) -> xr.DataArray:
    """Calculate all dihedral angles between the atoms specified.
    The atoms specified needed be bonded.

    Parameters
    ----------
    atXYZ
        A ``DataArray`` of coordinates, with ``atom`` and ``direction`` dimensions
    i, j, k, l
        The four atom indices
    deg
        Whether to return angles in degrees (True) or radians (False), by default False
    full
        Whether to return signed angles between -180° and 180° (True) or unsigned angles between 0 and 180° (False), by default False

    Returns
    -------
        A ``DataArray`` containing dihedral angles
    """
    a = atXYZ.isel(atom=i)
    b = atXYZ.isel(atom=j)
    c = atXYZ.isel(atom=k)
    d = atXYZ.isel(atom=l)
    result: xr.DataArray = full_dihedral_(a, b, c, d) if full else dihedral_(a, b, c, d)
    if deg:
        result = result * 180 / np.pi
    result.name = 'dihedral'
    result.attrs['long_name'] = r"$\varphi_{%d,%d,%d,%d}$" % (i, j, k, l)
    return result


@needs(dims={'atom'})
def angle(atXYZ: AtXYZ, i: int, j: int, k: int, *, deg: bool = False) -> xr.DataArray:
    a = atXYZ.isel(atom=i)
    b = atXYZ.isel(atom=j)
    c = atXYZ.isel(atom=k)
    ab = a - b
    cb = c - b
    result: xr.DataArray = angle_(ab, cb)
    if deg:
        result = result * 180 / np.pi
    result.name = 'angle'
    result.attrs['long_name'] = r"$\theta_{%d,%d,%d}$" % (i, j, k)
    return result


@needs(dims={'atom'})
def distance(atXYZ: AtXYZ, i: int, j: int) -> xr.DataArray:
    a = atXYZ.isel(atom=i)
    b = atXYZ.isel(atom=j)
    result: xr.DataArray = dnorm(a - b)
    result.name = 'distance'
    result.attrs['long_name'] = r"$\|\mathbf{r}_{%d,%d}\|$" % (i, j)
    return result

def ndarray_of_tuples(da, dim):
    da = da.transpose(..., dim)
    res = np.empty(len(da), dtype=object)
    list_of_tuples = [tuple(int(i) for i in x.data) for x in da]
    res[:] = list_of_tuples
    return res


def _check_matches(matches_or_mol, atXYZ):
    if matches_or_mol is None:
        mol = default_mol(atXYZ)
        matches = flag_bats(mol)[0]
    elif isinstance(matches_or_mol, rc.Mol):
        matches = flag_bats(mol)[0]
    elif isinstance(matches_or_mol, dict):
        matches = matches_or_mol
    else:
        raise TypeError(f"`matches_or_mol` of wrong type: {type(matches_or_mol)}")

    matches = {
        k: [t for t in v if t[0]] for k, v in matches.items()
    }  # remove unflagged
    return matches


def _positions(atXYZ, atom_idxs):
    return [
        atXYZ.sel(atom=list(idxs))
        .drop(['atom', 'atNames'], errors='raise')
        .rename(atom='descriptor')
        for idxs in zip(*atom_idxs)
    ]


def _assign_descriptor_coords(
    obj,
    atom_idxs: list[list[int]],
    bond_idxs: list[list[int]],
    bond_types,
    fragment_objs,
    tex_pattern,
):
    smiles = []
    smarts = []
    # without_bonds = []
    # with_bonds = []
    for m in fragment_objs:
        mol = rc.Mol(m)
        rc.RemoveStereochemistry(mol)
        smiles.append(rc.MolToSmiles(mol))
        smarts.append(rc.MolToSmarts(mol))
    #     without_bonds.append(
    #         ''.join(
    #             [get_symbol_from_atom_number(a.GetAtomicNum()) for a in mol.GetAtoms()]
    #         )
    #     )
    coords = xr.Coordinates(
        {
            'atom_indices': ('descriptor', np.fromiter(atom_idxs, dtype=object)),
            'bond_indices': ('descriptor', np.fromiter(bond_idxs, dtype=object)),
            'bond_orders': ('descriptor', np.fromiter(bond_types, dtype=object)),
            'fragment_smiles': ('descriptor', smiles),
            'fragment_smarts': ('descriptor', smarts),
            # 'descriptor_type': ('descriptor', bt),  # TODO
            # 'descriptor_type': ('descriptor', without_bonds), # TODO
            'descriptor_tex': (
                'descriptor',
                [tex_pattern % atom for atom in atom_idxs],
            ),
        }
        | {f'atom{i}': ('descriptor', list(x)) for i, x in enumerate(zip(*atom_idxs))}
    )

    return obj.assign_coords(coords).set_xindex('descriptor_tex')


@needs(dims={'atom', 'direction'})
def get_bond_lengths(
    atXYZ: xr.DataArray, matches_or_mol: dict | Mol | None = None
) -> xr.DataArray:
    """Identify bonds (using RDKit) and find the length of each bond in each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        `frame`, `atom` and `direction`
    matches_or_mol, optional
        A list containing information for each internal coordinate to be calculated.
        It may be convenient to use :py:func:`shnitsel.geo.geomatch.flag_dihedrals`
        to create a dictionary in the correct format, and then customize it.
        Alternatively, you may supply an RDKit ``Mol`` object, which is passed to
        :py:func:`shnitsel.geo.geomatch.flag_bats`.
        If this argument is omitted altogether, an RDKit ``Mol`` object is generated
        using :py:func:`shnitsel.bridges.default_mol` and used as above.
    Returns
    -------
        An :py:class:`xarray.DataArray` of bond lengths with dimensions `frame` and `bond`.

    Raises
    ------
    UserWarning
        If both `matches` and `mol` are specified.
    """
    matches = _check_matches(matches_or_mol, atXYZ)['bonds']

    _, atom_idxs, bond_idxs, bond_types, fragment_objs = zip(*matches)

    assert all(len(x) == 2 for x in atom_idxs)
    assert all(len(x) == 1 for x in bond_idxs)
    assert all(len(x) == 1 for x in bond_types)

    r0, r1 = _positions(atXYZ, atom_idxs)

    data = dnorm(r0 - r1)

    return _assign_descriptor_coords(
        data, atom_idxs, bond_idxs, bond_types, fragment_objs, '$r_{%d,%d}$'
    )


@needs(dims={'atom', 'direction'})
def get_bond_angles(
    atXYZ: xr.DataArray,
    matches_or_mol: dict | Mol | None = None,
    mol: Mol | None = None,
    ang: Literal[False, 'deg', 'rad'] = False,
):
    """Identify triples of bonded atoms (using RDKit) and calculate every bond angle for each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        `frame`, `atom` and `direction`
    matches_or_mol, optional
        A list containing information for each internal coordinate to be calculated.
        It may be convenient to use :py:func:`shnitsel.geo.geomatch.flag_angles`
        to create a dictionary in the correct format, and then customize it.
        Alternatively, you may supply an RDKit ``Mol`` object, which is passed to
        :py:func:`shnitsel.geo.geomatch.flag_bats`.
        If this argument is omitted altogether, an RDKit ``Mol`` object is generated
        using :py:func:`shnitsel.bridges.default_mol` and used as above.
    mol, optional
        An RDKit `Mol` object, which is generated from `atXYZ` if this argument is omitted.
    ang, optional
        If False (the default), returns sines and cosines;
        if set to 'deg', returns angles in degrees
        if set to 'rad', returns angles in radians

    Returns
    -------
        An :py:class:`xarray.DataArray` of bond angles with dimensions `frame` and `angle`.

    Raises
    ------
    UserWarning
        If both `matches` and `mol` are specified.
    """
    matches = _check_matches(matches_or_mol, atXYZ)['angles']

    _, atom_idxs, bond_idxs, bond_types, fragment_objs = zip(*matches)

    assert all(len(x) == 3 for x in atom_idxs)
    assert all(len(x) == 2 for x in bond_idxs)
    assert all(len(x) == 2 for x in bond_types)

    r0, r1, r2 = _positions(atXYZ, atom_idxs)

    std_args = (atom_idxs, bond_idxs, bond_types, fragment_objs)
    if ang:
        data = angle_(r0 - r1, r2 - r1)
        if ang == 'deg':
            data *= 180 / np.pi
        return _assign_descriptor_coords(data, *std_args, r"$\theta_{%d,%d,%d}$")
    else:
        cos, sin = angle_cos_sin_(r0 - r1, r2 - r1)
        cos = _assign_descriptor_coords(cos, *std_args, r"$\cos\theta_{%d,%d,%d}$")
        sin = _assign_descriptor_coords(sin, *std_args, r"$\sin\theta_{%d,%d,%d}$")
        return xr.concat([cos, sin], dim='descriptor')


@needs(dims={'atom', 'direction'})
def get_bond_torsions(
    atXYZ: xr.DataArray,
    matches_or_mol: dict | None = None,
    signed: bool | None = None,
    ang: Literal[False, 'deg', 'rad'] = False,
):
    """Identify quadruples of bonded atoms (using RDKit) and calculate the corresponding proper bond torsion for each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        `frame`, `atom` and `direction`
    matches_or_mol, optional
        A list containing information for each internal coordinate to be calculated.
        It may be convenient to use :py:func:`shnitsel.geo.geomatch.flag_dihedrals`
        to create a dictionary in the correct format, and then customize it.
        Alternatively, you may supply an RDKit ``Mol`` object, which is passed to
        :py:func:`shnitsel.geo.geomatch.flag_bats`.
        If this argument is omitted altogether, an RDKit ``Mol`` object is generated
        using :py:func:`shnitsel.bridges.default_mol` and used as above.
    signed, optional
        Whether to distinguish between clockwise and anticlockwise rotation,
        when returning angles; by default, do not distinguish.
    ang, optional
        If False (the default), returns sines and cosines;
        if set to 'deg', returns angles in degrees
        if set to 'rad', returns angles in radians

    Returns
    -------
        An :py:class:`xarray.DataArray` of bond torsions with dimensions `frame` and `torsion`.
    """
    matches = _check_matches(matches_or_mol, atXYZ)['dihedrals']

    _, atom_idxs, bond_idxs, bond_types, fragment_objs = zip(*matches)

    assert all(len(x) == 4 for x in atom_idxs)
    assert all(len(x) == 3 for x in bond_idxs)
    assert all(len(x) == 3 for x in bond_types)

    atom_positions = _positions(atXYZ, atom_idxs)
    std_args = (atom_idxs, bond_idxs, bond_types, fragment_objs)

    if ang:
        if signed:
            data = full_dihedral_(*atom_positions)
        else:
            data = dihedral_(*atom_positions)
        if ang == 'deg':
            data *= 180 / np.pi
        return _assign_descriptor_coords(
            data,
            *std_args,
            r"$\varphi_{%d,%d,%d,%d}$",
        )
    else:
        if signed is not None:
            raise ValueError("Can't use `signed` parameter when ang==False")

        r0, r1, r2, r3 = atom_positions
        n012 = normal(r0, r1, r2)
        n123 = normal(r1, r2, r3)
        cos, sin = angle_cos_sin_(n012, n123)
        cos = _assign_descriptor_coords(cos, *std_args, r"$\cos\varphi_{%d,%d,%d,%d}$")
        sin = _assign_descriptor_coords(sin, *std_args, r"$\sin\varphi_{%d,%d,%d,%d}$")
        return xr.concat([cos, sin], dim='descriptor')


@needs(dims={'atom', 'direction'})
def get_bats(
    atXYZ: xr.DataArray,
    matches_or_mol: dict | Mol | None = None,
    signed: bool | None = None,
    ang: Literal[False, 'deg', 'rad'] = False,
    pyr=False,
):
    """Get bond lengths, angles and torsions.

    Parameters
    ----------
    atXYZ
        The coordinates to use.
    matches_or_mol, optional
        An rdkit Mol object used to determine connectivity; by default this is
        determined automatically based on the first frame of ``atXYZ``.
        Alternatively, a dictionary containing match information, as produced
        by one of the ``flag_*`` functions.
    signed, optional
        Whether to distinguish between clockwise and anticlockwise rotation,
        when returning angles as opposed to cosine & sine values;
        by default, do not distinguish.
        NB. This applies only to the dihedrals, not to the three-center angles.
        The latter are always unsigned.
    ang, optional
        If False (the default), returns sines and cosines;
        if set to 'deg', returns angles in degrees
        if set to 'rad', returns angles in radians
    pyr
        Whether to include pyramidalizations from :py:func:`shnitsel.core.geom.get_pyramids`

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
            ['atom0', 'atom1', 'atom2', 'atom3', 'atNames'], errors='ignore'
        )

    return xr.concat(list(d.values()), dim='descriptor')


@needs(dims={'atom', 'direction'})
def identify_pyramids(mol: Mol) -> dict[int, list[int]]:
    res = {}
    for a in mol.GetAtoms():
        bonds = a.GetBonds()
        if len(bonds) != 3:
            continue

        current_idx = a.GetIdx()
        hydrogens = []
        non_hydrogens = []
        for b in bonds:
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            if a1.GetIdx() == current_idx:
                other = a2
            else:
                assert a2.GetIdx() == current_idx
                other = a1
            if other.GetAtomicNum() == 1:
                hydrogens.append(other.GetIdx())
            else:
                non_hydrogens.append(other.GetIdx())
        if len(hydrogens) == 2:  # Terminal double bond
            assert len(non_hydrogens) == 1
            plane_idxs = [hydrogens[0], current_idx, hydrogens[1]]
            bending_idx = non_hydrogens[0]
        else:  # Chain-internal double bond
            assert len(hydrogens) == 1
            assert len(non_hydrogens) == 2
            plane_idxs = [non_hydrogens[0], current_idx, non_hydrogens[1]]
            bending_idx = hydrogens[0]
        res[bending_idx] = plane_idxs
    return res


@needs(dims={'atom', 'direction'})
def pyramid_(a, b, c, x):
    abc = normal(a, b, c)
    return 0.5 * np.pi - angle_(abc, x - b)


def get_pyramids(
    atXYZ: xr.DataArray,
    pyramid_idxs: dict[int, list[int]] | None = None,
    mol: Mol | None = None,
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
    pyramid_idxs
        A dictionary containing types of bonds as keys, and lists of atom index pair
        as values. It may be convenient to use :py:func:`shnitsel.core.geom.identify_pyramids`
        to create a dictionary in the correct format, and then customize it. If omitted,
        relevant sets of atoms are identified automatically based on the ``mol`` argument.
    mol
        An RDKit ``Mol`` object, which is generated from ``atXYZ`` if this argument is omitted.
    deg
        Whether to return angles in degrees (as opposed to radians), by default False

    Returns
    -------
        An :py:class:`xarray.DataArray` of pyramidalizations with dimensions ``frame`` and ``atom``.

    Raises
    ------
    UserWarning
        If both `pyramid_idxs` and `mol` are specified.
    """
    if pyramid_idxs is None:
        if mol is None:
            mol = default_mol(atXYZ)
        pyramid_idxs = identify_pyramids(mol)
    elif mol is not None:
        raise UserWarning("pyramid_idxs passed, so mol will not be used")
    x, a, b, c = zip(*[[x, a, b, c] for x, (a, b, c) in pyramid_idxs.items()])
    data = [atXYZ.sel(atom=list(idxs)) for idxs in (a, b, c, x)]
    for i in range(len(data)):
        if 'atom' in data[i].indexes:
            data[i] = data[i].reset_index('atom')

    res = pyramid_(*data)
    res = res.rename(atom='descriptor')
    descriptor_tex = [
        r'$\chi_{%d,%d}^{%d,%d}$' % (b, x, a, c)
        for x, (a, b, c) in pyramid_idxs.items()
    ]

    atom_indices = np.empty(len(pyramid_idxs), dtype=object)
    atom_indices[:] = [(b, x, a, c) for x, (a, b, c) in pyramid_idxs.items()]

    res = res.assign_coords(
        descriptor_tex=('descriptor', descriptor_tex),
        descriptor_type=('descriptor', np.full(res.sizes['descriptor'], 'pyr')),
        atom_indices=('descriptor', atom_indices),
    ).set_xindex('descriptor_tex')
    if deg:
        res *= 180 / np.pi
    if not signed:
        res = np.abs(res)
    return res


def center_geoms(atXYZ, by_mass: Literal[False] = False):
    if by_mass:
        raise NotImplementedError
    return atXYZ - atXYZ.mean('atom')


def rotational_procrustes_(A, B, weight=None):
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


def rotational_procrustes(A, B, dim0='atom', dim1='direction', weight=None):
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
    atXYZ, reference_or_indexers: xr.DataArray | dict | None = None, **indexers_kwargs
):
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

    # atXYZ = center_geoms(atXYZ)
    # reference = center_geoms(reference)

    return rotational_procrustes(atXYZ, reference)
