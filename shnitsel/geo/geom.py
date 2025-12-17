"""\
    This module will contain two types of closely-related functionality:
    - generic geometry functions
    - functions that use RDKit to identify internal coordinates, and then the above to calculate the values
    Currently, the first category of function is found under postprocess
"""

from itertools import combinations, product
from typing import Literal

import numpy as np
from rdkit.Chem import Mol
import xarray as xr

from shnitsel.core._api_info import API, internal
from shnitsel.data.multi_indices import expand_midx
from shnitsel.data.atom_helpers import (
    get_atom_number_from_symbol,
    get_symbol_from_atom_number,
)
from ..analyze.generic import subtract_combinations, norm
from ..bridges import default_mol
from .._contracts import needs

from ..core.typedefs import AtXYZ


def dnorm(a: xr.DataArray | xr.Variable) -> xr.DataArray | xr.Variable:
    """Calculate the norm along the `direction` dimension. All other dimensions are maintaned

    Args:
        a (xr.DataArray | xr.Variable): The data array to perform the norming on

    Returns:
        xr.DataArray | xr.Variable: Resulting dataarray after calculation of the norm in the dimension `direction.
    """
    return norm(a, dim='direction')


def dcross(
    a: xr.DataArray | xr.Variable, b: xr.DataArray | xr.Variable
) -> xr.DataArray | xr.Variable:
    """Generalized cross vector product in the dimension of `direction`.

    Args:
        a (xr.DataArray | xr.Variable): The first array to use for the binary operation
        b (xr.DataArray | xr.Variable): The second array to use for the binary operation

    Returns:
        xr.DataArray | xr.Variable: The resulting array of the cross-product
    """
    return xr.cross(a, b, dim='direction')


def ddot(
    a: xr.DataArray | xr.Variable, b: xr.DataArray | xr.Variable
) -> xr.DataArray | xr.Variable:
    """Dot product in the dimension of `direction`.

    Args:
        a (xr.DataArray | xr.Variable): The first array to use for the binary operation
        b (xr.DataArray | xr.Variable): The second array to use for the binary operation

    Returns:
        xr.DataArray | xr.Variable: The resulting array of the dot-product still retaining all other dimensions except `direction`.
    """
    return xr.dot(a, b, dim='direction')


def angle_(
    a: xr.DataArray | xr.Variable, b: xr.DataArray | xr.Variable
) -> xr.DataArray | xr.Variable:
    """Helper function to calculate the angle between the entries in a and b based on their coordinates in the `direction` dimension.

    Args:
        a (xr.DataArray | xr.Variable): The first array to use for the binary operation
        b (xr.DataArray | xr.Variable): The second array to use for the binary operation

    Returns:
        xr.DataArray | xr.Variable: The resulting array of the angle calculation still retaining all other dimensions except `direction`.
    """
    return np.arccos(ddot(a, b) / (dnorm(a) * dnorm(b)))


def normal(
    a: xr.DataArray | xr.Variable,
    b: xr.DataArray | xr.Variable,
    c: xr.DataArray | xr.Variable,
) -> xr.DataArray | xr.Variable:
    """Calculate normal vectors on the planes through corresponding points in a, b and c.

    The normal vector will be calculated based on the position in `direction` dimension.

    Args:
        a (xr.DataArray | xr.Variable): The first array to use for the ternary operation
        b (xr.DataArray | xr.Variable): The second array to use for the ternary operation
        c (xr.DataArray | xr.Variable): The third array to use for the ternary operation

    Returns:
        xr.DataArray | xr.Variable: An array with all dimensions equal to those of a, b, and c but holding normal vectors along the `direction` dimension.
    """
    return dcross(a - b, c - b)


def dihedral_(
    a: xr.DataArray | xr.Variable,
    b: xr.DataArray | xr.Variable,
    c: xr.DataArray | xr.Variable,
    d: xr.DataArray | xr.Variable,
) -> xr.DataArray | xr.Variable:
    """Function to calculate the limited (up to +-\\pi/2 radian) dihedral angle between the points in arrays a,b,c and d.

    Args:
        a (xr.DataArray | xr.Variable): The first array to use for the operation
        b (xr.DataArray | xr.Variable): The second array to use for the operation
        c (xr.DataArray | xr.Variable): The third array to use for the operation
        d (xr.DataArray | xr.Variable): The fourth array to use for the operation

    Returns:
        xr.DataArray | xr.Variable: The array of dihedral angels between the four input arrays.
    """
    abc = normal(a, b, c)
    bcd = normal(b, c, d)
    return angle_(abc, bcd)


@internal()
def full_dihedral_(
    a: xr.DataArray | xr.Variable,
    b: xr.DataArray | xr.Variable,
    c: xr.DataArray | xr.Variable,
    d: xr.DataArray | xr.Variable,
) -> xr.DataArray | xr.Variable:
    """Function to calculate the signed/full dihedral angle (up to +-\\pi radian) between the points in arrays a,b,c and d.

    Args:
        a (xr.DataArray | xr.Variable): The first array to use for the operation
        b (xr.DataArray | xr.Variable): The second array to use for the operation
        c (xr.DataArray | xr.Variable): The third array to use for the operation
        d (xr.DataArray | xr.Variable): The fourth array to use for the operation

    Returns:
        xr.DataArray | xr.Variable: The array of full signed dihedral angels between the four input arrays.
    """
    abc = normal(a, b, c)
    bcd = normal(b, c, d)
    sign = np.sign(ddot(dcross(abc, bcd), (c - b)))
    return sign * angle_(abc, bcd)


@API()
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


@API()
@needs(dims={'atom'})
def angle(atXYZ: AtXYZ, i: int, j: int, k: int, *, deg: bool = False) -> xr.DataArray:
    """Method to calculate the angle between atoms i,j, and k in the positions DataArray throughout time.

    Can return results in radian (default) and degrees (if `deg=True`)

    Args:
        atXYZ (AtXYZ): DataArray with positions
        i (int): Index of first atom
        j (int): Index of second atom
        k (int): Index of third atom
        deg (bool, optional): Flag whether the results should be in degrees instead of radian. Defaults to False.

    Returns:
        xr.DataArray: The resulting angles between the denoted atoms.
    """
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


@API()
@needs(dims={'atom'})
def distance(atXYZ: AtXYZ, i: int, j: int) -> xr.DataArray:
    """Method to calculate the various distances between atoms i and j throughout time

    Args:
        atXYZ (AtXYZ): Array with atom positions
        i (int): Index of the first atom
        j (int): Index of the second atom

    Returns:
        xr.DataArray: The resulting array holding the pairwise distance between i and j.
    """
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


@internal()
def bond_type_to_symbols(element1: int, element2: int) -> str:
    """Helper function to get the Symbol string of the bond type denoted by elements 1 and 2.

    Args:
        element1 (int): First element in the bond
        element2 (int): Second element in the bond.

    Returns:
        str: Resulting String representation of elements in the bond type.
    """
    s1 = get_symbol_from_atom_number(element1)
    s2 = get_symbol_from_atom_number(element2)
    return s1 + s2


@internal()
def identify_bonds(
    mol: Mol, symbols: bool = True
) -> dict[tuple[int, int], list[tuple[int, int]]] | dict[str, list[tuple[int, int]]]:
    """Function to get a dict of all bond types and associated bonds.

    The bond type is either identified as a tuple of atom numbers/elements or as a string of the element names in the bond.

    Args:
        mol (Mol): The RDkit Mol object
        symbols (bool, optional): Flag whether the keys for the bond dictionary should be string represenations of the bond atoms instead. Defaults to True.

    Returns:
        dict[tuple[int, int], list[tuple[int, int]]] | dict[str, list[tuple[int, int]]]: The dictionary with the bond type (tuple of element numbers or string representation) as keys and the pairs of indices of the respective bonds as lists in the values.
    """
    bond_types: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for b in mol.GetBonds():
        a1 = b.GetBeginAtom()
        a2 = b.GetEndAtom()
        indices = (a1.GetIdx(), a2.GetIdx())
        elements = tuple(sorted([a1.GetAtomicNum(), a2.GetAtomicNum()]))
        if elements not in bond_types:
            bond_types[elements] = []
        bond_types[elements].append(indices)
    if symbols:
        return {bond_type_to_symbols(*k): v for k, v in bond_types.items()}
    return bond_types


@API()
@needs(dims={'atom', 'direction'})
def get_bond_lengths(
    atXYZ: xr.DataArray, bond_types=None, mol: Mol | None = None
) -> xr.DataArray:
    """Identify bonds (using RDKit) and find the length of each bond in each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        `frame`, `atom` and `direction`
    bond_types, optional
        A dictionary containing types of bonds as keys, and lists of atom index pair
        as values. It may be convenient to use :py:func:`shnitsel.core.geom.identify_bonds`
        to create a dictionary in the correct format, and then customize it. If omitted,
        bonds are identified automatically based on the `mol` argument.
    mol, optional
        An RDKit `Mol` object, which is generated from `atXYZ` if this argument is omitted.

    Returns
    -------
        An :py:class:`xarray.DataArray` of bond lengths with dimensions `frame` and `bond`.

    Raises
    ------
    UserWarning
        If both `bond_types` and `mol` are specified.
    """
    # dists = atXYZ.pipe(subtract_combinations, 'atom', labels=True).pipe(norm)
    if bond_types is None:
        if mol is None:
            mol = default_mol(atXYZ)
        bond_types = identify_bonds(mol, symbols=True)
    elif mol is not None:
        raise UserWarning("bond_types passed, so mol will not be used")

    a0 = []
    a1 = []
    bt = []
    for bond_type, bonds in bond_types.items():
        for bond in bonds:
            a0.append(bond[0])
            a1.append(bond[1])
            bt.append(bond_type)

    atom_indices = np.empty(len(bt), dtype=object)
    atom_indices[:] = list(zip(a0, a1))

    sel = atXYZ.sel(atom=a0)
    return (
        dnorm(sel.copy(data=atXYZ.sel(atom=a0).data - atXYZ.sel(atom=a1).data))
        .drop(['atom', 'atNames'], errors='ignore')
        .rename(atom='descriptor')
        .assign_coords(
            {
                'atom0': ('descriptor', a0),
                'atom1': ('descriptor', a1),
                'atom_indices': ('descriptor', atom_indices),
                'descriptor_type': ('descriptor', bt),
                'descriptor_tex': (
                    'descriptor',
                    [r'$r_{%d,%d}$' % atom for atom in zip(a0, a1)],
                ),
            }
        )
        .set_xindex('descriptor_tex')
    )


def identify_angles(mol: Mol) -> xr.Dataset:
    """Helper function to generate a dataset containing the angles within the molecules,
    bond types, angle types, atom types, etc. for a rigorous angle analysis.

    Args:
        mol (Mol): The molecule as RDkit Molecule to derive the angles for.

    Returns:
        xr.Dataset: The dataset holding the full angular information.
            Will contain variables:
             - `at_idx`: The ids of the atoms in the molecule
             - `at_nums`: The elements of the atoms in the molecule
             - `bond_type`: The types of bonds involved in the respective angle
             - `angle_type`: str representation of the bond types and elements involved in the angle.
             - `angle_symbol`: LaTeX-printable angle label including the atom indices.
            Will have the dimensions:
             - `angle`: to iterate through all angles
             - `atom`: To iterate through the atoms as parts of the respective angle (not all atoms in the molecule). Length 3.
             - `bond`: To iterate through the bonds involved in the respective angle (not all bonds in the molecule). Length 2.
    """
    triples = []
    at_nums = []
    bond_types = []
    angle_types = []
    angle_symbols = []
    for a in mol.GetAtoms():
        bonds = a.GetBonds()
        if len(bonds) < 2:
            continue
        a0 = a.GetIdx()
        neighbors = set()
        for b in bonds:
            a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            neighbors |= {a1, a2}
        neighbors -= {a0}
        for a1, a2 in combinations(neighbors, 2):
            triple = a1, a0, a2
            n1, n0, n2 = (mol.GetAtomWithIdx(i).GetAtomicNum() for i in triple)
            b10 = int(mol.GetBondBetweenAtoms(a1, a0).GetBondType())
            b02 = int(mol.GetBondBetweenAtoms(a0, a2).GetBondType())
            # Swap order so end-atom with lowest atomic number is first;
            # if end-atoms have same atomic number, order by bond-type:
            if (n1 > n2) or (n1 == n2 and b10 > b02):
                a1, a2 = a2, a1
                n1, n2 = n2, n1
                b10, b02 = b02, b10
            triples.append((a1, a0, a2))
            at_nums.append((n1, n0, n2))
            bond_types.append((b10, b02))
            s = get_symbol_from_atom_number
            angle_types.append(f"{s(n1)}{b10}{s(n0)}{b02}{s(n2)}")
            angle_symbols.append(r"$\theta_{%d,%d,%d}$" % (a1, a0, a2))

    return xr.Dataset(
        {
            'at_idx': (('descriptor', 'atom'), np.array(triples)),
            'at_num': (('descriptor', 'atom'), np.array(at_nums)),
            'bond_type': (('descriptor', 'bond'), np.array(bond_types)),
            'descriptor_type': ('descriptor', angle_types),
            'descriptor_tex': ('descriptor', angle_symbols),
        }
    )


@API()
@needs(dims={'atom', 'direction'})
def get_bond_angles(
    atXYZ: xr.DataArray,
    angle_types: xr.Dataset | None = None,
    mol: Mol | None = None,
    deg: bool = False,
):
    """Identify triples of bonded atoms (using RDKit) and calculate every bond angle for each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        `frame`, `atom` and `direction`
    angle_types, optional
        An :py:class:`xarray.Dataset` containing atom indices, atomic numbers, bond types, angle type
        and angle symbol for each angle to be calculated.
        It may be convenient to use :py:func:`shnitsel.core.geom.identify_angles`
        to create a Dataset in the correct format, and then customize it. If omitted,
        angles are identified automatically based on the `mol` argument.
    mol, optional
        An RDKit `Mol` object, which is generated from `atXYZ` if this argument is omitted.
    deg, optional
        Whether to return angles in degrees (as opposed to radians), by default False

    Returns
    -------
        An :py:class:`xarray.DataArray` of bond angles with dimensions `frame` and `angle`.

    Raises
    ------
    UserWarning
        If both `angle_types` and `mol` are specified.
    """
    if angle_types is None:
        if mol is None:
            mol = default_mol(atXYZ)
        angle_types = identify_angles(mol)
    elif mol is not None:
        raise UserWarning("angle_types passed, so mol will not be used")

    al = atXYZ.sel(atom=angle_types.at_idx[:, 0])
    ac = atXYZ.sel(atom=angle_types.at_idx[:, 1])
    ar = atXYZ.sel(atom=angle_types.at_idx[:, 2])

    def f(var, n):
        xs = var.data
        # res = np.empty(len(xs), dtype=object)
        # res[:] = [tuple(x) for x in xs]
        # return res
        # return var.data.astype('f,f')
        dtype = ','.join(['i'] * n)
        return np.array([tuple(x) for x in xs], dtype=dtype)

    at_idxs = {f'atom{i}': angle_types['at_idx'].isel(atom=i) for i in range(3)}
    angles = (
        angle_(al - ac, ar - ac)
        .assign_coords(
            # at_nums=('angle', angle_types['at_num'].astype('i,i,i').data),
            # bond_types=('angle', angle_types['bond_type'].astype('i,i').data),
            # at_nums=('angle', f(angle_types['at_num'], 3)),
            # bond_types=('angle', f(angle_types['bond_type'], 2)),
            descriptor_type=angle_types['descriptor_type'],
            descriptor_tex=angle_types['descriptor_tex'],
            atom_indices=(
                'descriptor',
                ndarray_of_tuples(angle_types['at_idx'], 'atom'),
            ),
            **at_idxs,
        )
        .set_xindex('descriptor_tex')
    )
    if deg:
        angles *= 180 / np.pi
    return angles


@internal()
def identify_torsions(mol: Mol) -> xr.Dataset:
    """Method to find sets of four contiguous atoms that form a torsion.

    Args:
        mol (Mol): The molecule as RDkit Molecule to derive the angles for.

    Returns:
        xr.Dataset: Dataset containing the full torsion information.
            Will contain variables:
             - `at_idx`: The ids of the atoms in the molecule part of the respective torsion.
             - `at_nums`: The elements of the atoms in the molecule part of the respective torsion.
             - `bond_type`: The types of bonds involved in the respective torsion.
             - `torsion_type`: str representation of the bond types and elements involved in the torsion.
             - `torsion_symbol`: LaTeX-printable torsion label including the atom indices.
            Will have the dimensions:
             - `torsion`: to iterate through all torsions
             - `atom`: To iterate through the atoms as parts of the respective torsion (not all atoms in the molecule). Length 4.
             - `bond`: To iterate through the bonds involved in the respective torsion (not all bonds in the molecule). Length 3.
    """
    quadruples = []
    at_nums = []
    bond_types = []
    torsion_types = []
    torsion_symbols = []
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() > a2.GetAtomicNum():
            a1, a2 = a2, a1
        i1, i2 = a1.GetIdx(), a2.GetIdx()
        bonds1, bonds2 = a1.GetBonds(), a2.GetBonds()
        if len(bonds1) < 2 or len(bonds2) < 2:
            continue
        for b1, b2 in product(bonds1, bonds2):
            atoms = [
                b1.GetBeginAtom(),
                b1.GetEndAtom(),
                b2.GetBeginAtom(),
                b2.GetEndAtom(),
            ]
            idxs = [a.GetIdx() for a in atoms]
            i0 = idxs[0] if i1 == idxs[1] else idxs[1]
            i3 = idxs[3] if i2 == idxs[2] else idxs[2]
            if i0 == i2 or i3 == i1:
                continue
            idxs = [i0, i1, i2, i3]
            atoms = [mol.GetAtomWithIdx(i) for i in idxs]
            an = [a.GetAtomicNum() for a in atoms]
            bt = [int(bond.GetBondType()) for bond in [b1, b, b2]]
            if an[1] == an[2] and (
                an[0] > an[3]  #
                or (an[0] == an[3] and bt[0] > bt[2])  #
            ):
                idxs.reverse()
                atoms.reverse()
                an.reverse()
                b1, b2 = b2, b1
                bt.reverse()
            quadruples.append(idxs)
            at_nums.append(an)
            bond_types.append(bt)
            torsion_types.append(
                '{}'.join([get_symbol_from_atom_number(n) for n in an]).format(*bt)
            )
            torsion_symbols.append(r"$\varphi_{%d,%d,%d,%d}$" % (i0, i1, i2, i3))
    return xr.Dataset(
        {
            'at_idx': (('descriptor', 'atom'), np.array(quadruples)),
            'at_num': (('descriptor', 'atom'), np.array(at_nums)),
            'bond_type': (('descriptor', 'bond'), np.array(bond_types)),
            'descriptor_type': ('descriptor', torsion_types),
            'descriptor_tex': ('descriptor', torsion_symbols),
        }
    )


@API()
@needs(dims={'atom', 'direction'})
def get_bond_torsions(
    atXYZ: xr.DataArray,
    quadruple_types: xr.Dataset | None = None,
    mol: Mol | None = None,
    signed: bool = False,
    deg: bool = False,
):
    """Identify quadruples of bonded atoms (using RDKit) and calculate the corresponding proper bond torsion for each
    frame.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        `frame`, `atom` and `direction`
    quadruple_types, optional
        An :py:class:`xarray.Dataset` containing atom indices, atomic numbers, bond types, torsion type
        and torsion symbol for each angle to be calculated.
        It may be convenient to use :py:func:`shnitsel.core.geom.identify_torsions`
        to create a Dataset in the correct format, and then customize it. If omitted,
        angles are identified automatically based on the `mol` argument.
    mol, optional
        An RDKit `Mol` object, which is generated from `atXYZ` if this argument is omitted.
    signed, optional
        Whether to distinguish between clockwise and anticlockwise rotation, by default False
    deg, optional
        Whether to return angles in degrees (as opposed to radians), by default False

    Returns
    -------
        An :py:class:`xarray.DataArray` of bond torsions with dimensions `frame` and `torsion`.

    Raises
    ------
    UserWarning
        If both `torsion_types` and `mol` are specified.
    """
    if quadruple_types is None:
        if mol is None:
            mol = default_mol(atXYZ)
        quadruple_types = identify_torsions(mol)
    elif mol is not None:
        raise UserWarning("quadruple_types passed, so mol will not be used")
    if 'atNames' in atXYZ.coords or 'atNames' in atXYZ:
        atXYZ = atXYZ.drop_vars('atNames')
    atom_positions = [atXYZ.sel(atom=quadruple_types.at_idx[:, i]) for i in range(4)]
    if signed:
        res = full_dihedral_(*atom_positions)
    else:
        res = dihedral_(*atom_positions)
    if deg:
        res *= 180 / np.pi
    at_idxs = {f'atom{i}': quadruple_types['at_idx'].isel(atom=i) for i in range(4)}
    return res.assign_coords(
        descriptor_type=quadruple_types['descriptor_type'],
        descriptor_tex=quadruple_types['descriptor_tex'],
        atom_indices=(
            'descriptor',
            ndarray_of_tuples(quadruple_types['at_idx'], 'atom'),
        ),
        **at_idxs,
    ).set_xindex('descriptor_tex')


@needs(dims={'atom', 'direction'})
def get_bats(
    atXYZ: xr.DataArray,
    mol: Mol | None = None,
    signed: bool = False,
    deg: bool = False,
    pyr=False,
):
    """Get bond lengths, angles and torsions.

    Parameters
    ----------
    atXYZ
        The coordinates to use.
    mol, optional
        An rdkit Mol object used to determine connectivity; by default this is
        determined automatically based on the first frame of ``atXYZ``.
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
    if mol is None:
        mol = default_mol(atXYZ)
    d = {
        'bond': get_bond_lengths(atXYZ, mol=mol),
        'angle': get_bond_angles(atXYZ, mol=mol, deg=deg),
        'torsion': get_bond_torsions(atXYZ, mol=mol, signed=signed, deg=deg),
    }

    d['bond'] = d['bond'].drop_vars(['atom0', 'atom1'])
    d['angle'] = d['angle'].drop_vars(['atom0', 'atom1', 'atom2'])
    d['torsion'] = d['torsion'].drop_vars(['atom0', 'atom1', 'atom2', 'atom3'])

    if pyr:
        d['pyr'] = get_pyramids(atXYZ, mol=mol, deg=deg, signed=signed)
        if 'atNames' in d['pyr'].coords:
            d['pyr'] = d['pyr'].drop_vars('atNames')

    to_concat = [d['bond'], d['angle'], d['torsion']]

    if pyr:
        to_concat.append(d['pyr'])

    return xr.concat(to_concat, dim='descriptor')


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
def pyramid_(
    a: xr.DataArray, b: xr.DataArray, c: xr.DataArray, x: xr.DataArray
) -> xr.DataArray:
    """Method to calculate the pyramidalization angle of a quadruple of atoms.

    The result will be \\pi/2 minus the angle between the normal of ABC and the vector BX.
    (I.e.: the pyramidalization at atom b)

    Args:
        a (xr.DataArray): The first array to use for the operation
        b (xr.DataArray): The second array to use for the operation
        c (xr.DataArray): The third array to use for the operation
        x (xr.DataArray): The fourth array to use for the operation

    Returns:
        xr.DataArray: The array of pyramidalization angles
    """
    abc = normal(a, b, c)
    return 0.5 * np.pi - angle_(abc, x - b)  # type: ignore # Cannot be a variable if provided with a DataArray


@API()
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


@API()
@needs(dims={'atom'})
def center_geoms(atXYZ: AtXYZ, by_mass: Literal[False] = False) -> AtXYZ:
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


def rotational_procrustes_(A, B, weight=None):
    # TODO: FIXME: Document
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
    # TODO: FIXME: Document
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
    # TODO: FIXME: Document
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
