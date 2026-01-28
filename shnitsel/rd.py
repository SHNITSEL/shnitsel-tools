"""This module contains functions that accept an RDKit.Chem.Mol object;
but *not* necessarily functions that *return* a Mol object."""

from typing import Literal, Mapping, Sequence, TYPE_CHECKING

import rdkit.Chem as rc
import rdkit.Chem.rdDetermineBonds  # noqa: F401
import matplotlib as mpl
import numpy as np

if TYPE_CHECKING:
    from shnitsel.filtering.structure_selection import FeatureDescriptor

#################################################
# Functions for converting RDKit objects to
# SMILES annotated with the original atom indices
# to maintain the order in the `atom` index


def set_atom_props(
    mol: rc.Mol,
    inplace: bool = False,
    **kws: Sequence[str | int]
    | Literal[True]
    | Literal[False]
    | Mapping[int, int | str]
    | None,
) -> rc.Mol | None:
    """Set properties on atoms of an ``rdkit.Chem.Mol`` object

    Parameters
    ----------
    mol : rdkit.Chem.mol
        The ``Mol`` object
    inplace : bool, optional
        Whether to alter ``mol``; , by default False (returns a copy)
    **kws : Sequence[str | int] | Literal[True] | Mapping[int, int | str] | None
        A mapping where parameter names represent the name of a property
        and the arguments are either
            - a dict mapping the atom indices to values that should be assigned. Missing atom indices are ignored
            - a sequence of str or int values the atoms should be set to;
            - ``True``, in which case the atom indices as
              assigned by RDKit will be used as values;
            - ``False``, in which the property will be cleared on every atom.
            - ``None`` values, which are simply ignored

    Returns
    -------
        A copy of the ``Mol`` object, if ``inplace=False``, otherwise the provided mol

    Raises
    ------
    ValueError
        If the amount of values passed to the properties kwargs did not match the amount of
        atoms in the mol.
    """
    if not inplace:
        mol = rc.Mol(mol)
    natoms = mol.GetNumAtoms()
    for prop, vals in kws.items():
        if vals is None:
            # None is not assigned, just ignored
            continue
        elif isinstance(vals, dict):
            # Try and assign the values to each atom identified by the key
            for atom_id, val in vals.items():
                atom = mol.GetAtomWithIdx(atom_id)
                if atom is not None:
                    atom.SetProp(prop, str(val))
        else:
            if vals is True:
                # atom indices are values if `vals=True`
                vals = range(natoms)
            elif vals is False:
                # `False` means, the value should be cleared.
                for atom in mol.GetAtoms():
                    atom.ClearProp(prop)
                continue
            elif natoms != len(vals):
                raise ValueError(
                    f"{len(vals)} values were passed for {prop}, but 'mol' has {natoms} atoms"
                )

            for atom, val in zip(mol.GetAtoms(), vals):
                atom.SetProp(prop, str(val))
    return mol


def mol_to_numbered_smiles(mol: rc.Mol) -> str:
    """Generate a SMILES string containing mapping numbers
    corresponding to the atom indices in the mol object

    Parameters
    ----------
    mol
        An ``rdkit.Chem.Mol`` object

    Returns
    -------
        A SMILES string

    Notes
    -----
        This is intended as a way to store the connectivity
        and order of a matrix of coordinates
    """
    mol = rc.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return rc.MolToSmiles(mol)


def highlight_pairs(mol: rc.Mol, feature_indices: Sequence["FeatureDescriptor"]):
    """Highlight specified pairs of atoms in an image of an ``rdkit.Chem.Mol`` object

    Parameters
    ----------
    mol : rc.Mol
        The ``Mol`` object
    feature_indices : Sequence[FeatureDescriptor]
        A list of tuples of indices for various features.

    Returns
    -------
        Raw PNG data
    """
    d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(320, 240)
    # colors = iter(mpl.colormaps['tab10'](range(10)))
    colors = iter(mpl.colormaps['rainbow'](np.linspace(0, 1, len(feature_indices))))

    acolors: dict[int, list[tuple[float, float, float]]] = {}
    bonds: dict[int, list[tuple[float, float, float]]] = {}
    for feature in feature_indices:
        c = tuple(next(colors))
        if isinstance(feature, int):
            # Position
            a = feature
            if a not in acolors:
                acolors[a] = []
            acolors[a].append(c)
        elif isinstance(feature, tuple):
            flen = len(feature)
            if flen == 2:
                if isinstance(feature[1], int):
                    # Bond
                    a1, a2 = feature
                    if a1 < 0 or a2 < 0:
                        continue
                    if (bond := mol.GetBondBetweenAtoms(a1, a2)) is not None:
                        bondid = bond.GetIdx()
                        if bondid not in bonds:
                            bonds[bondid] = []
                        bonds[bond.GetIdx()].append(c)
                    else:
                        for a in [a1, a2]:
                            if a not in acolors:
                                acolors[a] = []
                            acolors[a].append(c)
                elif isinstance(feature[1], tuple):
                    # pyramid
                    a1, (a2, a3, a4) = feature
                    if a1 < 0 or a2 < 0 or a3 < 0 or a4 < 0:
                        continue
                    # Mark bonds
                    for other in (a2, a3, a4):
                        if (bond := mol.GetBondBetweenAtoms(a1, other)) is not None:
                            bondid = bond.GetIdx()
                            if bondid not in bonds:
                                bonds[bondid] = []
                            bonds[bond.GetIdx()].append(c)
                    # Mark all atoms
                    for a in [a1, a2, a3, a4]:
                        if a not in acolors:
                            acolors[a] = []
                        acolors[a].append(c)
            elif flen == 3 or flen == 4:
                # angle or dihedral
                # Mark bonds
                if not all(x >= 0 for x in feature):
                    continue

                for i in range(flen - 1):
                    a1, a2 = feature[i], feature[i + 1]
                    if (bond := mol.GetBondBetweenAtoms(a1, a2)) is not None:
                        bondid = bond.GetIdx()
                        if bondid not in bonds:
                            bonds[bondid] = []
                        bonds[bond.GetIdx()].append(c)

                # Mark all atoms
                for a in feature:
                    if a not in acolors:
                        acolors[a] = []
                    acolors[a].append(c)

    # d.drawOptions().fillHighlights = False
    d.drawOptions().setBackgroundColour((0.8, 0.8, 0.8, 0.5))
    d.drawOptions().padding = 0

    d.DrawMoleculeWithHighlights(mol, '', acolors, bonds, {}, {}, -1)
    d.FinishDrawing()
    return d.GetDrawingText()
