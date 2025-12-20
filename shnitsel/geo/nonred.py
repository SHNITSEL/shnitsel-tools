import logging

import rdkit.Chem as rc

from shnitsel.bridges import set_atom_props
from shnitsel.geo.geomatch import _get_bond_info as get_bond_info


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
    set_atom_props(mol, inplace=True, original_index=True)

    mol_no_hs = rc.RemoveHs(mol)
    # The following call causes the _smilesAtomOutputOrder property to be computed and set:
    _ = rc.MolToSmiles(mol_no_hs)
    props = mol_no_hs.GetPropsAsDict(includePrivate=True, includeComputed=True)
    order = list(props['_smilesAtomOutputOrder'])
    return [mol_no_hs.GetAtomWithIdx(i).GetIntProp('original_index') for i in order]



def flag_nonredundant(mol: rc.Mol, include_h: bool = True):
    """
    Compute a non-redundant set of bonds, angles, and dihedrals
    sufficient to uniquely determine the atoms of the input,
    given a fixed centre and whole-molecular orientation.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule under study.
    include_h : bool
        Whether to include internal coordinates for hydrogen atoms

    Returns
    -------
    dict
        {
            'bonds':      [...],
            'angles':     [...],
            'dihedrals':  [...]
        }
    """

    def f(s):
        return ' '.join(str(x) for x in s)

    logger = logging.getLogger('flag_nonredundant')
    order = get_smiles_order_ignoring_h(mol)
    if include_h:
        order.extend(set(range(mol.GetNumAtoms())).difference(order))

    bonds = []
    angles = []
    dihedrals = []
    runs = {}
    min_run_len = 0
    for i in order:
        if len(runs) == 0:
            logger.info(f'Atom {i}: Nothing to do')
            runs[i] = [i]
            continue

        neigh_runs = (
            runs.get(neighbor.GetIdx(), [])
            for neighbor in mol.GetAtomWithIdx(i).GetNeighbors()
        )
        runs[i] = run = max(neigh_runs, key=lambda x: len(x)) + [i]

        assert len(run) >= min_run_len

        if len(run) > 4:
            logger.info(f"Atom {i}: Using run ({f(run[:-4])}) {f(run[-4:])}")
        else:
            logger.info(f"Atom {i}: Using run {f(run)}")

        for n, k in enumerate(run[:-1]):
            if len(runs.get(k, [])) < 4 <= len(run) - n:
                new_run = run[n:][::-1]
                logger.info(f"Overwriting run for {k} with {f(new_run)}")
                runs[k] = new_run

        if min_run_len < 4 and len(run) > min_run_len:
            min_run_len = len(run)
            logger.info(f'{min_run_len=}')

        if len(run) >= 2:
            bonds.append((1, tuple(run[-2:])))
        if len(run) >= 3:
            angles.append((1, tuple(run[-3:])))
        if len(run) >= 4:
            dihedrals.append((1, tuple(run[-4:])))

    return {
        'bonds': get_bond_info(mol, bonds),
        'angles': get_bond_info(mol, angles),
        'dihedrals': get_bond_info(mol, dihedrals),
    }