"""\
    This module will contain two types of closely-related functionality:
    - generic geometry functions
    - functions that use RDKit to identify internal coordinates, and then the above to calculate the values
    Currently, the first category of function is found under postprocess
"""

from itertools import combinations, product

import numpy as np
from rdkit.Chem import Mol
import xarray as xr

from .parse.common import __atnum2symbol__
from .postprocess import (
    distance as distance,
    angle as angle,
    dihedral as dihedral,
    angle_,
    dihedral_,
    full_dihedral_,
    subtract_combinations,
    norm,
    default_mol,
)
from .xrhelpers import expand_midx


def bond_type_to_symbols(e1, e2):
    s1 = __atnum2symbol__[e1]
    s2 = __atnum2symbol__[e2]
    return s1 + s2


def identify_bonds(mol, symbols=True):
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


def get_bond_lengths(atXYZ, bond_types=None, mol=None):
    dists = atXYZ.pipe(subtract_combinations, 'atom', labels=True).pipe(norm)
    if bond_types is None:
        if mol is None:
            mol = default_mol(atXYZ)
        bond_types = identify_bonds(mol, symbols=True)
    elif mol is not None:
        raise UserWarning("bond_types passed, so mol will not be used")
    return (
        # Change how this works! We don't need to calculate them separately!
        xr.concat(
            [
                dists.sel(atomcomb=bonds).pipe(
                    expand_midx, 'atomcomb', 'bond_type', bond_type
                )
                for bond_type, bonds in bond_types.items()
            ],
            dim='atomcomb',
        )
        .rename({'from': 'atom1', 'to': 'atom2', 'atomcomb': 'bond'})
        .transpose('frame', ...)
    )


def identify_angles(mol: Mol):
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
            s = __atnum2symbol__
            angle_types.append(f"{s[n1]}{b10}{s[n0]}{b02}{s[n2]}")
            angle_symbols.append(r"$\theta_{%d,%d,%d}$" % (a1, a0, a2))

    return xr.Dataset(
        {
            'at_idx': (('angle', 'atom'), np.array(triples)),
            'at_num': (('angle', 'atom'), np.array(at_nums)),
            'bond_type': (('angle', 'bond'), np.array(bond_types)),
            'angle_type': ('angle', angle_types),
            'angle_symbol': ('angle', angle_symbols),
        }
    )


def get_bond_angles(atXYZ, angle_types=None, mol=None):
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

    angles = (
        angle_(al - ac, ar - ac)
        .assign_coords(
            # at_nums=('angle', angle_types['at_num'].astype('i,i,i').data),
            # bond_types=('angle', angle_types['bond_type'].astype('i,i').data),
            # at_nums=('angle', f(angle_types['at_num'], 3)),
            # bond_types=('angle', f(angle_types['bond_type'], 2)),
            angle_type=angle_types['angle_type'],
            angle_symbol=angle_types['angle_symbol'],
        )
        .set_xindex('angle_symbol')
    )
    return angles


def identify_torsions(mol: Mol) -> xr.Dataset:
    """Finds sets of four contiguous atoms that form a torsion"""
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
            s = __atnum2symbol__
            torsion_types.append('{}'.join([s[n] for n in an]).format(*bt))
            torsion_symbols.append(r"$\varphi_{%d,%d,%d,%d}$" % (i0, i1, i2, i3))
    return xr.Dataset(
        {
            'at_idx': (('torsion', 'atom'), np.array(quadruples)),
            'at_num': (('torsion', 'atom'), np.array(at_nums)),
            'bond_type': (('torsion', 'bond'), np.array(bond_types)),
            'torsion_type': ('torsion', torsion_types),
            'torsion_symbol': ('torsion', torsion_symbols),
        }
    )


def get_bond_torsions(atXYZ, quadruple_types=None, mol=None, signed=False):
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
    return res.assign_coords(
        torsion_type=quadruple_types['torsion_type'],
        torsion_symbol=quadruple_types['torsion_symbol'],
    )