"""This module contains functions that accept an RDKit.Chem.Mol object;
but *not* necessarily functions that *return* a Mol object."""

import rdkit.Chem as rc
import rdkit.Chem.rdDetermineBonds  # noqa: F401


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


