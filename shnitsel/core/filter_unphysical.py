from logging import warning

import numpy as np
import xarray as xr

from rdkit import Chem as rc
from rdkit.Chem import rdDetermineBonds, AllChem

from . import postprocess as P, xrhelpers as xh
from .plot import pca_biplot
from .postprocess import (
    mol_to_numbered_smiles as mol_to_numbered_smiles,
    numbered_smiles_to_mol,
)


def find_bonds_by_element(mol, elem1: int, elem2: int):
    def elems_correct(b):
        atnums = {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()}
        return atnums == {elem1, elem2}

    def indices(b):
        return (b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx())

    return [indices(b) for b in mol.GetBonds() if elems_correct(b)]

def mol_from_atXYZ(atXYZ_frame, charge=0, covFactor=1.5, to2D=True):
    mol = rc.rdmolfiles.MolFromXYZBlock(P.to_xyz(atXYZ_frame))
    # rdDetermineBonds.DetermineConnectivity(mol) # 2025-02-03 TODO Unify!
    rdDetermineBonds.DetermineBonds(
        mol, charge=charge, useVdw=True, covFactor=covFactor
    )
    if to2D:
        AllChem.Compute2DCoords(mol)  # type: ignore
    return mol

def max_bond_lengths(atXYZ, elem1=1, elem2=6):
    def dists(a1, a2):
        return P.norm(atXYZ.isel(atom=a1) - atXYZ.isel(atom=a2), dim='direction')

    if 'smiles_map' in atXYZ.attrs:
        mol = numbered_smiles_to_mol(atXYZ.attrs['smiles_map'])
    else:
        warning("`smiles_map` attribute missing; falling back to frame 0")
        mol = mol_from_atXYZ((atXYZ.isel(frame=0)))

    bonds = find_bonds_by_element(mol, elem1, elem2)
    maxlengths = xr.concat([
      dists(a1, a2).groupby('trajid').map(np.max)
      for a1, a2 in bonds ], dim='bond')
    atoms1, atoms2 = zip(*bonds)
    maxlengths.coords['atom1'] = 'bond',list(atoms1)
    maxlengths.coords['atom2'] = 'bond',list(atoms2)
    maxlengths = maxlengths.set_xindex(['atom1', 'atom2'])
    return maxlengths

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
    noodle_da = P.pairwise_dists_pca(atXYZ)
    noodle = noodle_da.to_dataset('PC').rename({0: 'PC1', 1: 'PC2'})
    maskfn = maskfn or (lambda data: data.PC1**2 + data.PC2**2 > 1.5)
    mask = maskfn(noodle)
    return np.unique(noodle.sel(frame=mask).trajid)

def show_bonds_mol(mol, elem1, elem2, to2D):
    pairs = find_bonds_by_element(mol, elem1, elem2)
    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", str(atom.GetIdx()))
    return pca_biplot.highlight_pairs(mol, pairs)


# TODO: refactor -- different way of specifying which bonds
def filter_cleavage(frames, *, CC=False, CH=False, CN=False, NH=False, verbose=2):
    try:
        from IPython.display import display, Image
    except ImportError:
        can_show = False
    else:
        can_show = True

    def show(elem1, elem2):
        mol = numbered_smiles_to_mol(frames.atXYZ.attrs['smiles_map'])
        display(Image(show_bonds_mol(mol, elem1, elem2, to2D=True)))

    def act(descr, elem1, elem2, cutoff):
        nonlocal frames
        overlong = find_overlong(frames.atXYZ, elem1, elem2, cutoff=cutoff)
        if verbose:
            print(f"Found following {descr} bonds:")
            if can_show:
                show(elem1, elem2)
            ntraj = len(np.unique(overlong))
            nframes = xh.sel_trajids(frames, overlong).sizes['frame']
            print(
                f"Remove {ntraj} trajectories ({nframes} frames) containing {descr} cleavage"
            )
            if verbose >= 2:
                print("with IDs:", overlong)
        frames = xh.sel_trajids(frames, overlong, invert=True)

    if CC:
        act('C-C', 6, 6, 2.8)
    if CH:
        act('C-H', 1, 6, 1.7)
    if CN:
        act('C-N', 6, 7, 1.54 * 2)
    if NH:
        act('C-N', 1, 7, 1.056 * 2)

    if verbose:
        ntraj = len(np.unique(frames.trajid))
        nframes = frames.sizes['frame']
        print(f"Keep {ntraj} trajectories ({nframes} frames)")

    return frames

# TODO 2025-06-16: Does this function belong here?
def smiles_map(atXYZ_frame, charge=0, covFactor=1.5) -> str:
    mol = mol_from_atXYZ(atXYZ_frame, charge=charge, covFactor=covFactor, to2D=True)
    return mol_to_numbered_smiles(mol)


#########################
# Non-geometric filtering


def exclude_involving_state(frames, state):
    return exclude_trajs(frames, frames.trajid[frames.astate == 3])


def cutoffs(mask_da):
    return (
        mask_da.groupby('trajid')
        .apply(
            lambda traj: traj.coords['time'][-1]
            if traj.all()
            else traj.sel(frame=~traj).coords['time'][0]
        )
        .rename(trajid='trajid_')
    )


def all_cutoffs(frames):
    e = frames[['e_kin']]
    e['e_pot'] = frames['energy'].sel(state=frames.astate).drop_vars('state')
    e['e_tot'] = e['e_pot'] + e['e_kin']

    feat = xr.Dataset()
    feat['etot_drift'] = (
        e['e_tot'].groupby('trajid').apply(lambda traj: abs(traj - traj.isel(frame=0)))
    )
    feat['ekin_step'] = P.sudi(e['e_kin'])
    feat['epot_step'] = P.sudi(e['e_pot'])
    feat['etot_step'] = P.sudi(e['e_tot'])
    feat['is_hop'] = P.sudi(frames['astate']) != 0

    c = xr.Dataset()
    c['cutoff_length'] = (
        'trajid_',
        [traj.coords['time'][-1] for _, traj in feat.groupby('trajid')],
    )
    c['cutoff_etot_window'] = cutoffs(abs(feat['etot_drift']) < 0.2)
    c['cutoff_etot_step'] = cutoffs(abs(feat['etot_step']) < 0.1)
    c['cutoff_epot_step'] = cutoffs((abs(feat['epot_step']) < 0.7) | feat['is_hop'])
    c['cutoff_ekin_step'] = cutoffs((abs(feat['ekin_step']) < 0.7) | feat['is_hop'])
    c['cutoff_hop_epot'] = cutoffs((abs(feat['epot_step']) < 1.0) | ~feat['is_hop'])

    typenames = list(c.data_vars)
    cda = c.to_dataarray('cutoff')
    c['cutoff_min'] = cda.min('cutoff')
    c['cutoff_type'] = cda.argmin('cutoff')

    c.attrs['types'] = typenames
    return c