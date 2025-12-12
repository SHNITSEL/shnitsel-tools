from functools import reduce
from itertools import chain
from operator import and_

import numpy as np
import rdkit.Chem as rc
import xarray as xr

from shnitsel.data.multi_indices import expand_midx
from shnitsel.bridges import default_mol, set_atom_props
from shnitsel.clean.common import is_stacked  # TODO: move


def atom_path_to_bond_path(mol, atoms):
    from itertools import combinations

    res = []
    for i, j in combinations(atoms, 2):
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is not None:
            res.append(bond.GetIdx())
    return res


def substruct_match_to_submol(mol, substruct_match):
    bond_path = atom_path_to_bond_path(mol, substruct_match)
    return rc.PathToSubmol(mol, bond_path)


def list_analogs(ensembles, smarts='', vis=False):
    if vis:
        from IPython import display

    mols = [default_mol(x) for x in ensembles]
    if not smarts:
        from rdkit.Chem import rdFMCS

        smarts = rdFMCS.FindMCS(mols).smartsString

    search = rc.MolFromSmarts(smarts)

    results = []
    for compound, mol in zip(ensembles, mols):
        set_atom_props(mol, molAtomMapNumber=True)
        idxs = list(mol.GetSubstructMatch(search))

        res_mol = substruct_match_to_submol(mol, idxs)
        set_atom_props(res_mol, atomNote=True)
        # The following ensures that atoms will be stored in canonical order
        # this is as opposed to using the substructure match directly
        idxs = [a.GetAtomMapNum() for a in res_mol.GetAtoms()]

        if vis:
            display(res_mol)

        range_ = range(len(idxs))
        results.append(
            compound.isel(atom=idxs)
            .assign_coords(atom=range_)
            .sortby('atom')
            .assign_attrs(mol=res_mol)
        )
    return results


def combine_compounds_unstacked(compounds, names=None):
    coord_names = [set(x.coords) for x in compounds]
    coords_shared = reduce(and_, coord_names)
    compounds = [
        x.drop_vars(set(x.coords).difference(coords_shared)) for x in compounds
    ]
    if names is None:
        names = range(len(compounds))
    compounds = [
        x.assign_coords(
            {
                'compound': ('trajid', np.full(x.sizes['trajid'], name)),
                'traj': x.trajid,
            }
        )
        .reset_index('trajid')
        .set_xindex(['compound', 'traj'])
        for x, name in zip(compounds, names)
    ]

    return xr.concat(compounds, dim='trajid')


def combine_compounds_stacked(compounds, names=None):
    concat_dim = 'frame'

    coord_names = [set(x.coords) for x in compounds]
    coords_shared = reduce(and_, coord_names)
    compounds = [
        x.drop_vars(set(x.coords).difference(coords_shared)) for x in compounds
    ]

    if names is None:
        names = range(len(compounds))
    compounds = [
        expand_midx(x, 'frame', 'compound', name).drop_dims('trajid_')
        for x, name in zip(compounds, names)
    ]

    # TODO: Should we set `coords` param? to what?
    # Current code displays a `FutureWarning` for xarray==2025.11.0
    res = xr.concat(compounds, dim=concat_dim)

    trajid_only = xr.concat(
        [
            obj.drop_vars(
                [
                    k
                    for k, v in chain(obj.data_vars.items(), obj.coords.items())
                    if 'trajid_' not in v.dims
                ]
            )
            for obj in compounds
        ],
        dim='trajid_',
    )
    res = res.merge(trajid_only, join='exact')

    if any('time_' in x.dims for x in compounds):
        time_only = xr.concat(
            [obj.drop_dims(['frame', 'trajid_'], errors='ignore') for obj in compounds],
            dim='time_',
        )
        res = res.assign(time_only)
    return res


def combine_analogs(ensembles, smarts='', names=None, vis=False):
    analogs = list_analogs(ensembles, smarts=smarts, vis=False)
    if all(is_stacked(x) for x in analogs):
        res = combine_compounds_stacked(analogs, names=names)
    elif not any(is_stacked(x) for x in analogs):
        res = combine_compounds_unstacked(analogs, names=names)
    else:
        raise ValueError("Inconsistent formats")
    # del res.attrs['mol']
    return res.assign_attrs(mols=[x.attrs['mol'] for x in analogs])