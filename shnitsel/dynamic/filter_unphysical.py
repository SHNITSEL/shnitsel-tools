import numpy as np
import scipy.stats as st
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import rdkit as rd
from rdkit import Chem as rkc
from rdkit.Chem import rdDetermineBonds, AllChem

from . import postprocess as P, xrhelpers as xh

def is_c_h_bond(b):
    atnums = {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()}
    return atnums == {1, 6}

def find_c_h_bonds(mol):
    def indices(b):
        return (b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx())

    return [indices(b) for b in mol.GetBonds() if is_c_h_bond(b)]

def mol_from_atXYZ(atXYZ_frame):
    mol = rkc.rdmolfiles.MolFromXYZBlock(P.to_xyz(atXYZ_frame))
    rdDetermineBonds.DetermineConnectivity(mol)
    AllChem.Compute2DCoords(mol)
    return mol

def max_ch_lengths(atXYZ):
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

def lengths_sorted(atXYZ):
    lengths = max_ch_lengths(atXYZ)
    return lengths.sortby(lengths.sum(dim='bond'))

def find_overlong(atXYZ, cutoff=2):
    lengths = lengths_sorted(atXYZ)
    mask = (lengths>cutoff).any('bond')
    return lengths.trajid.sel(trajid=mask).values

def exclude_trajs(frames, trajids):
    if isinstance(trajids, set):
        trajids = list(trajids)
    return frames.sel(frame=~frames.trajid.isin(trajids))

def exclude_overlong(frames, cutoff=2):
    return exclude_trajs(frames, find_overlong(frames.atXYZ, cutoff=cutoff))   

def find_eccentric(atXYZ, maskfn=None):
    noodle = P.pairwise_dists_pca(atXYZ)
    noodle = noodle.to_dataset('PC').rename({0: 'PC1', 1: 'PC2'})
    maskfn = maskfn or (lambda data: data.PC1**2 + data.PC2**2 > 1.5)
    mask = maskfn(noodle)
    return np.unique(noodle.sel(frame=mask).trajid)

def exclude_involving_state(frames, state):
    return exclude_trajs(frames, frames.trajid[frames.astate==3])