"""This module contains functions that accept an RDKit.Chem.Mol object;
but *not* necessarily functions that *return* a Mol object."""

from typing import Literal

import rdkit.Chem as rc
import rdkit.Chem.rdDetermineBonds  # noqa: F401
import matplotlib as mpl
import numpy as np

#################################################
# Functions for converting RDKit objects to
# SMILES annotated with the original atom indices
# to maintain the order in the `atom` index


def set_atom_props(
    mol: rc.Mol, inplace: bool = False, **kws: list[str] | Literal[True]
) -> rc.Mol | None:
    """Set properties on atoms of an ``rdkit.Chem.Mol`` object

    Parameters
    ----------
    mol : rdkit.Chem.mol
        The ``Mol`` object
    inplace : bool, optional
        Whether to alter ``mol``; , by default False (returns a copy)
    **kws
        A mapping where parameter names represent the name of a property
        and the arguments are either

            - a list of str values the atoms should be set to;
            - ``True``, in which case the atom indices as
              assigned by RDKit will be used as values;
            - ``False``, in which the property will be cleared on every atom.

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
            continue
        elif vals is True:
            vals = range(natoms)
        elif vals is False:
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


def highlight_pairs(mol: rc.Mol, pairs: list[tuple[int, int]]):
    """Highlight specified pairs of atoms in an image of an ``rdkit.Chem.Mol`` object

    Parameters
    ----------
    mol
        The ``Mol`` object
    pairs
        A list of pairs of atom indices

    Returns
    -------
        Raw PNG data
    """
    d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(320, 240)
    # colors = iter(mpl.colormaps['tab10'](range(10)))
    colors = iter(mpl.colormaps['rainbow'](np.linspace(0, 1, len(pairs))))

    acolors: dict[int, list[tuple[float, float, float]]] = {}
    bonds = {}
    for a1, a2 in pairs:
        if (bond := mol.GetBondBetweenAtoms(a1, a2)) is not None:
            bonds[bond.GetIdx()] = [(1, 0.5, 0.5)]
        else:
            c = tuple(next(colors))
            for a in [a1, a2]:
                if a not in acolors:
                    acolors[a] = []
                acolors[a].append(c)

    # d.drawOptions().fillHighlights = False
    d.drawOptions().setBackgroundColour((0.8, 0.8, 0.8, 0.5))
    d.drawOptions().padding = 0

    d.DrawMoleculeWithHighlights(mol, '', acolors, bonds, {}, {}, -1)
    d.FinishDrawing()
    return d.GetDrawingText()
