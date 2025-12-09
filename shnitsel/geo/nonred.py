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


# def runs_for_order(mol, order):


def flag_nonredundant(mol, include_h=True):
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
    expectation = 0
    for i in order:
        if len(runs) == 0:
            logger.info(f'Atom {i}: Nothing to do')
            runs[i] = [i]
            continue

        run_lens = [
            ((j := neighbor.GetIdx()), len(runs.get(j, [])))
            for neighbor in mol.GetAtomWithIdx(i).GetNeighbors()
        ]
        j, rl = max(run_lens, key=lambda tup: tup[1])
        assert rl + 1 >= expectation

        run = runs[j] + [i]
        runs[i] = run

        if len(run) > 4:
            logger.info(f"Atom {i}: Using run ({f(run[:-4])}) {f(run[-4:])}")
        else:
            logger.info(f"Atom {i}: Using run {f(run)}")

        for n, k in enumerate(run[:-1]):
            if len(runs.get(k, [])) < 4 <= len(run) - n:
                new_run = run[n:][::-1]
                logger.info(f"Overwriting run for {k} with {f(new_run)}")
                runs[k] = new_run

        if expectation < 4 and len(run) > expectation:
            expectation = len(run)
            logger.info(f'{expectation=}')

        if len(run) >= 2:
            bonds.append((1, tuple(run[-2:])))
        if len(run) >= 3:
            angles.append((1, tuple(run[-3:])))
        if len(run) >= 4:
            dihedrals.append((1, tuple(run[-4:])))
    # for run in runs.values():

    return {
        'bonds': get_bond_info(mol, bonds),
        'angles': get_bond_info(mol, angles),
        'dihedrals': get_bond_info(mol, dihedrals),
    }