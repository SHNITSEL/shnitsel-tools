from functools import reduce
from itertools import combinations
from operator import and_
from typing import Any, Iterable

import numpy as np
import rdkit.Chem as rc
import xarray as xr

from shnitsel.data.multi_indices import expand_midx
from shnitsel.bridges import default_mol, set_atom_props
from shnitsel.clean.common import is_stacked  # TODO: move


def _find_atom_pairs(mol, atoms):
    # TODO: Might we need this elsewhere?
    res = []
    for i, j in combinations(atoms, 2):
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is not None:
            res.append(bond.GetIdx())
    return res


def _substruct_match_to_submol(mol, substruct_match):
    bond_path = _find_atom_pairs(mol, substruct_match)
    return rc.PathToSubmol(mol, bond_path)


def list_analogs(
    ensembles: Iterable[xr.DataArray], smarts: str = '', vis: bool = False
) -> Iterable[xr.DataArray]:
    """Extract a common moiety from a selection of ensembles

    Parameters
    ----------
    ensembles
        An ``Iterable`` of ``xr.DataArray``s, each containing the geometries of an ensemble of
        trajectories for a different compound; they
    smarts, optional
        A SMARTS-string indicating the moiety to cut out of each compound;
        in each case, the match returned by :py:func:`rdkit.Chem.Mol.GetSubstrucMatch`
        (not necessarily the only possible match) will be used;
        if no SMARTS is provided, a minimal common submol will be extracted using
        ``rdFMCS.FindMCS``
    vis, optional
        Whether to display a visual indication of the match

    Returns
    -------
       An ``Iterable`` of ``xr.DataArray``s
    """
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

        res_mol = _substruct_match_to_submol(mol, idxs)
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
            # .drop_vars(['atNames', 'atNums'], errors='ignore') # TODO
        )
    return results


def _combine_compounds_unstacked(compounds, names=None, concat_kws=None):
    if concat_kws is None:
        concat_kws = {}

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

    return xr.concat(compounds, dim='trajid', **concat_kws)


def _combine_compounds_stacked(compounds, names=None, concat_kws=None):
    if concat_kws is None:
        concat_kws = {}

    concat_dim = 'frame'

    coord_names = [set(x.coords) for x in compounds]
    coords_shared = reduce(and_, coord_names)
    compounds = [
        x.drop_vars(set(x.coords).difference(coords_shared)) for x in compounds
    ]

    if names is None:
        names = range(len(compounds))

    # Which coords are unique on a compound level? So far:
    c_per_compound = ['atNames', 'atNums']

    per_compound = {
        crd: xr.concat(
            [obj[crd] for obj in compounds],
            dim='compound_',
        )
        for crd in c_per_compound
        if all(crd in obj.coords for obj in compounds)
    }

    compounds = [
        expand_midx(x, 'frame', 'compound', name)
        .drop_dims('trajid_')
        .drop_vars(c_per_compound, errors='ignore')
        for x, name in zip(compounds, names)
    ]

    res = xr.concat(
        compounds, dim=concat_dim, **({'combine_attrs': 'drop_conflicts'} | concat_kws)
    )
    res = res.assign_coords(per_compound)
    res = res.assign_coords(compound_=names)

    if any('time_' in x.dims for x in compounds):
        time_only = xr.concat(
            [obj.drop_dims(['frame', 'trajid_'], errors='ignore') for obj in compounds],
            dim='time_',
            **concat_kws,
        )
        res = res.assign(time_only)

    # TODO: consider using MajMinIndex
    return res


def combine_analogs(
    ensembles: Iterable[xr.DataArray],
    smarts: str = '',
    names: Iterable[str] | None = None,
    vis: bool = False,
    *,
    concat_kws: dict[str, Any] = None,
) -> xr.DataArray:
    """Combine ensembles for different compounds by finding the
    moieties they have in common

    Parameters
    ----------
    ensembles
        An ``Iterable`` of ``xr.DataArray``s, each containing the geometries of an ensemble of
        trajectories for a different compound; these trajectories should all
        be in the same format, i.e.:

            - all stacked (with 'frames' dimension indexed by'trajid' and 'time' MultiIndex levels)
            - all unstacked (with independent 'trajid' and 'time' dimensions)

    smarts
        A SMARTS-string indicating the moiety to cut out of each compound;
        in each case, the match returned by :py:func:`rdkit.Chem.Mol.GetSubstrucMatch`
        (not necessarily the only possible match) will be used;
        if no SMARTS is provided, a minimal common submol will be extracted using
        ``rdFMCS.FindMCS``
    names
        An ``Iterable`` of ``Hashable`` to identify the compounds;
        these values will end up in the ``compound`` coordinate, by default None
    vis
        Whether to display a visual indication of the match, by default False
    concat_kws
        Keyword arguments for internal calls to ``xr.concat``

    Returns
    -------
        An xr.Dataset of trajectories, with a MultiIndex level identifying each
        trajectory by its compound name (or index, if no names were provided)
        and trajid

    Raises
    ------
    ValueError
        If the ensembles provided are in a mixture of formats (i.e. some have trajectories
        stacked, others unstacked)
    """
    analogs = list_analogs(ensembles, smarts=smarts, vis=False)
    if all(is_stacked(x) for x in analogs):
        res = _combine_compounds_stacked(analogs, names=names, concat_kws=concat_kws)
    elif not any(is_stacked(x) for x in analogs):
        res = _combine_compounds_unstacked(analogs, names=names, concat_kws=concat_kws)
    else:
        raise ValueError("Inconsistent formats")

    mols = [x.attrs['mol'] for x in analogs]
    res = res.assign_attrs(
        mol=mols[0],  # TODO: Try replacing with search pattern object
    ).assign_coords(mols=('compound_', mols))
    return res