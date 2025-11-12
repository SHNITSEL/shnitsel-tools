"""This module contains functions that accept an RDKit.Chem.Mol object;
but *not* necessarily functions that *return* a Mol object."""

import rdkit.Chem as rc
import rdkit.Chem.rdDetermineBonds  # noqa: F401
import matplotlib as mpl
import numpy as np

#################################################
# Functions for converting RDKit objects to
# SMILES annotated with the original atom indices
# to maintain the order in the `atom` index


def set_atom_props(mol, **kws):
    natoms = mol.GetNumAtoms()
    for prop, vals in kws.items():
        if vals is None:
            continue
        elif vals is True:
            vals = range(natoms)
        elif natoms != len(vals):
            raise ValueError(
                f"{len(vals)} values were passed for {prop}, but 'mol' has {natoms} atoms"
            )

        for atom, val in zip(mol.GetAtoms(), vals):
            atom.SetProp(prop, str(val))
    return mol


def mol_to_numbered_smiles(mol: rc.Mol) -> str:
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return rc.MolToSmiles(mol)


def highlight_pairs(mol, pairs):
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