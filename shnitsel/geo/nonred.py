import logging

import rdkit.Chem as rc

from shnitsel.bridges import set_atom_props
from shnitsel.geo.geomatch import __get_bond_info as get_bond_info


def get_smiles_order_ignoring_h(mol: rc.Mol) -> list[int]:
    """Returns the order in which atoms would appear in the canonical SMILES of
    ``mol``, ignoring hydrogens

    Parameters
    ----------
    mol
        An ``rdkit.Chem.Mol`` object

    Returns
    -------
        A list of integers representing indices of the original ``mol`` object (as opposed
        to the integers assigned to the copy stripped of hydrogens)
    """
    # Avoid mutating input
    mol = rc.Mol(mol)
    # molAtomMapNumber would interfere with the canonicalization, so use custom property
    set_atom_props(mol, original_index=True)

    mol_no_hs = rc.RemoveHs(mol)
    # The following call causes the _smilesAtomOutputOrder property to be computed and set:
    _ = rc.MolToSmiles(mol_no_hs)
    props = mol_no_hs.GetPropsAsDict(includePrivate=True, includeComputed=True)
    order = list(props['_smilesAtomOutputOrder'])
    return [mol_no_hs.GetAtomWithIdx(i).GetIntProp('original_index') for i in order]


def flag_nonredundant(mol):
    """
    Compute a non-redundant set of bonds, angles, and dihedrals
    sufficient to locate the non-hydrogen atoms of the input.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule under study.

    Returns
    -------
    dict
        {
            'bonds':      [...],
            'angles':     [...],
            'dihedrals':  [...]
        }
    """
    logger = logging.getLogger('flag_nonredundant')
    order = get_smiles_order_ignoring_h(mol)
    bonds = []
    angles = []
    dihedrals = []
    runs = {}
    for i in order:
        logger.info(f'Working on atom {i}')
        if len(runs) == 0:
            logger.info('Nothing to do')
            runs[i] = [i]
            continue
        run_lens = [
            ((j := neighbor.GetIdx()), len(runs.get(j, [])))
            for neighbor in mol.GetAtomWithIdx(i).GetNeighbors()
        ]
        j = max(run_lens, key=lambda tup: tup[1])[0]
        logger.info(f'Maximum run len was {len(runs[j])} starting back from atom {j}')
        run = runs[j] + [i]
        if len(run) >= 2:
            logger.info(f'Add bond: {run[-2:]!r}')
            bonds.append((1, tuple(run[-2:])))
        if len(run) >= 3:
            logger.info(f'Add angle: {run[-3:]!r}')
            angles.append((1, tuple(run[-3:])))
        if len(run) >= 4:
            logger.info(f'Add dihedral: {run[-4:]!r}')
            dihedrals.append((1, tuple(run[-4:])))
        runs[i] = run

    # TODO: deal with the hydrogens here

    return {
        'bonds': get_bond_info(mol, bonds),
        'angles': get_bond_info(mol, angles),
        'dihedrals': get_bond_info(mol, dihedrals),
    }