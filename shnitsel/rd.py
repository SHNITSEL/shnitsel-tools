from typing import Literal

import rdkit.Chem as rc
import rdkit.Chem.rdDetermineBonds  # noqa: F401
import xarray as xr

from ._contracts import needs
from .core.postprocess import to_xyz

from .core.datasheet.structure import mol_to_png as mol_to_png

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


@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def to_mol(
    atXYZ_frame: xr.DataArray,
    charge: int | None = None,
    covFactor: float = 1.2,
    to2D: bool = True,
    molAtomMapNumber: list | Literal[True] | None = None,
    atomNote: list | Literal[True] | None = None,
    atomLabel: list | Literal[True] | None = None,
) -> rc.Mol:
    """Convert a single frame's geometry to an RDKit Mol object

    Parameters
    ----------
    atXYZ_frame
        The ``xr.DataArray`` object to be converted; must have 'atom' and 'direction' dims,
        must not have 'frame' dim.
    charge
        Charge of the molecule, used by RDKit to determine bond orders; if ``None`` (the default),
        this function will try ``charge=0`` and leave the bond orders undetermined if that causes
        an error; otherwise failure to determine bond order will raise an error.
    covFactor
        Scales the distance at which atoms are considered bonded, by default 1.2
    to2D
        Discard 3D information and generate 2D conformer (useful for displaying), by default True
    molAtomMapNumber
        Set the ``molAtomMapNumber`` properties to values provided in a list,
        or (if ``True`` is passed) set the properties to the respective atom indices
    atomNote
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomNote`` properties
    atomLabel
        Behaves like the ``molAtomMapNumber`` parameter above, but for the ``atomLabel`` properties

    Returns
    -------
        An RDKit Mol object

    Raises
    ------
    ValueError
        If ``charge`` is not ``None`` and bond order determination fails
    """
    mol = rc.rdmolfiles.MolFromXYZBlock(to_xyz(atXYZ_frame))
    rc.rdDetermineBonds.DetermineConnectivity(mol, useVdw=True, covFactor=covFactor)
    try:
        rc.rdDetermineBonds.DetermineBondOrders(mol, charge=(charge or 0))
    except ValueError as err:
        if charge is not None:
            raise err
    if to2D:
        rc.rdDepictor.Compute2DCoords(mol)  # type: ignore
    return set_atom_props(
        mol, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel
    )


def mol_to_numbered_smiles(mol: rc.Mol) -> str:
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return rc.MolToSmiles(mol)


def numbered_smiles_to_mol(smiles: str) -> rc.Mol:
    mol = rc.MolFromSmiles(smiles, sanitize=False)  # sanitizing would strip hydrogens
    map_new_to_old = [-1 for i in range(mol.GetNumAtoms())]
    for atom in mol.GetAtoms():
        # Renumbering with e.g. [3, 2, 0, 1] means atom 3 gets new index 0, not vice-versa!
        map_new_to_old[int(atom.GetProp("molAtomMapNumber"))] = atom.GetIdx()
    return rc.RenumberAtoms(mol, map_new_to_old)


@needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
def smiles_map(atXYZ_frame, charge=0, covFactor=1.5) -> str:
    mol = to_mol(atXYZ_frame, charge=charge, covFactor=covFactor, to2D=True)
    return mol_to_numbered_smiles(mol)


def default_mol(obj) -> rc.Mol:
    if 'atXYZ' in obj:  # We have a frames Dataset
        atXYZ = obj['atXYZ']
    else:
        atXYZ = obj  # We have an atXYZ DataArray

    if 'smiles_map' in obj.attrs:
        return numbered_smiles_to_mol(obj.attrs['smiles_map'])
    elif 'smiles_map' in atXYZ.attrs:
        return numbered_smiles_to_mol(atXYZ.attrs['smiles_map'])

    try:
        charge = obj.attrs.get('charge', 0)
        return to_mol(atXYZ.isel(frame=0), charge=charge)
    except (KeyError, ValueError):
        raise ValueError(
            "Failed to get default mol, please set a smiles map. "
            "For example, if the compound has charge c and frame i contains a representative geometry, use "
            "frames.attrs['smiles_map'] = frames.atXYZ.isel(frame=i).sh.get_smiles_map(charge=c)"
        )