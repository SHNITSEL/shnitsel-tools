import numpy as np
import xarray as xr

from rdkit import Chem as rkc
from rdkit.Chem import rdDetermineBonds, AllChem

from . import postprocess as P


def find_bonds_by_element(mol, elem1: int, elem2: int):
    def elems_correct(b):
        atnums = {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()}
        return atnums == {elem1, elem2}

    def indices(b):
        return (b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx())

    return [indices(b) for b in mol.GetBonds() if elems_correct(b)]


def is_c_h_bond(b):
    atnums = {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()}
    return atnums == {1, 6}

def find_c_h_bonds(mol):
    def indices(b):
        return (b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx())

    return [indices(b) for b in mol.GetBonds() if is_c_h_bond(b)]

def mol_from_atXYZ(atXYZ_frame, charge=0, covFactor=1.5, to2D=True):
    mol = rc.rdmolfiles.MolFromXYZBlock(P.to_xyz(atXYZ_frame))
    # rdDetermineBonds.DetermineConnectivity(mol) # 2025-02-03 TODO Unify!
    rdDetermineBonds.DetermineBonds(
        mol, charge=charge, useVdw=True, covFactor=covFactor
    )
    if to2D:
        AllChem.Compute2DCoords(mol)
    return mol

def mol_to_numbered_smiles(mol):
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return rc.MolToSmiles(mol)


def numbered_smiles_to_mol(smiles):
    mol = rc.MolFromSmiles(smiles, sanitize=False)  # sanitizing would strip hydrogens
    map_new_to_old = [-1 for i in range(mol.GetNumAtoms())]
    for atom in mol.GetAtoms():
        # Renumbering with e.g. [3, 2, 0, 1] means atom 3 gets new index 0, not vice-versa!
        map_new_to_old[int(atom.GetProp("molAtomMapNumber"))] = atom.GetIdx()
    return rc.RenumberAtoms(mol, map_new_to_old)


def max_bond_lengths(atXYZ, elem1=1, elem2=6):
    def dists(a1, a2):
        return P.norm(atXYZ.isel(atom=a1) - atXYZ.isel(atom=a2), dim='direction')
    mol = mol_from_atXYZ((atXYZ.isel(frame=0)))
    bonds = find_c_h_bonds(mol)
    maxlengths = xr.concat([
      dists(a1, a2).groupby('trajid').map(np.max)
      for a1, a2 in bonds ], dim='bond')
    atoms1, atoms2 = zip(*bonds)
    maxlengths.coords['atom1'] = 'bond',list(atoms1)
    maxlengths.coords['atom2'] = 'bond',list(atoms2)
    maxlengths = maxlengths.set_xindex(['atom1', 'atom2'])
    return maxlengths

def max_ch_lengths(atXYZ):
    return max_bond_lengths(atXYZ, elem1=1, elem2=6)


def lengths_sorted(atXYZ, elem1=1, elem2=6):
    lengths = max_bond_lengths(atXYZ, elem1, elem2)
    return lengths.sortby(lengths.sum(dim='bond'))

def find_overlong(atXYZ, elem1=1, elem2=6, cutoff=2):
    lengths = lengths_sorted(atXYZ, elem1, elem2)
    mask = (lengths>cutoff).any('bond')
    return lengths.trajid.sel(trajid=mask).values

def exclude_trajs(frames, trajids):
    if isinstance(trajids, set):
        trajids = list(trajids)
    return frames.sel(frame=~frames.trajid.isin(trajids))

def exclude_overlong(frames, cutoff=2):
    return exclude_trajs(frames, find_overlong(frames.atXYZ, cutoff=cutoff))   

def find_eccentric(atXYZ, maskfn=None):
    if not isinstance(atXYZ, xr.DataArray):
        raise TypeError()
    noodle = P.pairwise_dists_pca(atXYZ)
    noodle = noodle.to_dataset('PC').rename({0: 'PC1', 1: 'PC2'})
    maskfn = maskfn or (lambda data: data.PC1**2 + data.PC2**2 > 1.5)
    mask = maskfn(noodle)
    return np.unique(noodle.sel(frame=mask).trajid)

def exclude_involving_state(frames, state):
    return exclude_trajs(frames, frames.trajid[frames.astate==3])