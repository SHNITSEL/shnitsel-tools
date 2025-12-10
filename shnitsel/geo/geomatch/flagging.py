import logging
from logging import warning, info

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

from shnitsel.filtering.structure_selection import StructureSelection
from .features.angles import flag_angles
from .features.atoms import flag_atoms
from .features.bonds import flag_bonds
from .features.dihedrals import flag_dihedrals
from .presentation import __get_highlight_molimg, __get_img_multiple_mols
from shnitsel.vis.colormaps import st_yellow


def __match_pattern(mol: Mol, smarts: str) -> list[tuple] | None:
    """
    Find all substructure matches of a SMARTS pattern in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.
    smarts : str
        SMARTS pattern to search for.

    Returns
    -------
    list of tuples
        Each tuple contains atom indices corresponding to one match of the SMARTS pattern.
        Returns an empty list if no match is found.
    None
        if the provided SMARTS string was invalid.
    """
    pattern = Chem.MolFromSmarts(smarts)

    if pattern is None:
        info(f"Invalid SMARTS '{smarts}'. Falling back to full reference molecule.")
        matches = None
    else:
        matches = list(mol.GetSubstructMatches(pattern))

    return matches


def __get_bond_info(mol: Mol, flagged_tuples: StructureSelection):
    """
    Extend flagged tuple of bonds, angles or dihedrals by a tuple of bond indices and
    an information of their respective bond types as double (1.0: SINGLE, 1.5: AROMATIC,
    2.0: DOUBLE, and 3.0: TRIPLE)

    Parameters
    ----------
    mol : RDKit Mol
        Molecule object.
    flagged_tuples : list of tuples
        Each entry: (flag, atom_tuple)
        Example: [(1, (0,1,2,3)), (0, (1,2,3,4))]

    Returns
    -------
    flagged_tuples_bond_info: list of tuples containing additional bond information
        Each entry: (flag, atom_tuple, bond_tuple, bondtype_tuple,
        rdkit.Chem.rdchem.Mol of submolecule/subgraph )
        The bond tuple is the tuple of all bond indices in the RDKit Mol.
        The bondtype tuple is the tuple of all bond types aligned with the indices
        in `bond_tuple`.
        Example: [(0, (5, 0, 1), (6, 0), (1.0, 2.0), rdkit.Chem.rdchem.Mol object)]
    """

    l_flags = []
    l_atom_idxs = []
    l_bond_idxs = []
    l_bond_types = []
    l_submols = []

    for flag, t in flagged_tuples:
        if len(t) < 2:
            continue
        l_flags.append(flag)
        l_atom_idxs.append(t)

        # consecutive pairs (divide angles, torsions back to bonds
        atom_pairs = [(t[i], t[i + 1]) for i in range(len(t) - 1)]

        inner_bond_idxs = []
        inner_bond_types = []
        for i, j in atom_pairs:
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond:
                inner_bond_idxs.append(bond.GetIdx())
                inner_bond_types.append(bond.GetBondTypeAsDouble())

        l_bond_idxs.append(tuple(inner_bond_idxs))
        l_bond_types.append(tuple(inner_bond_types))
        l_submols.append(Chem.PathToSubmol(mol, inner_bond_idxs))

    flagged_tuples_bond_info = list(
        zip(l_flags, l_atom_idxs, l_bond_idxs, l_bond_types, l_submols)
    )

    return flagged_tuples_bond_info


def __collect_tuples(selection: dict):
    """
    Select the appropriate list of tuples based on hierarchy:
    dihedrals -> angles -> bonds.

    Each value is a list of tuples of the form:
        (flag, atom_tuple, bond_tuple, bondtype_tuple, Mol)
    """

    hierarchy = ["dihedrals", "angles", "bonds"]
    flagged_tuples = {}
    for key in hierarchy:
        if key in selection:
            tuples = selection[key]
            if any(entry[0] != 0 for entry in tuples):
                flagged_tuples = tuples
            else:
                warning(
                    f"All flags in {key} are zero. " "Will collect all atoms instead."
                )
                if 'atoms' in selection:
                    if any(entry[0] != 0 for entry in tuples):
                        flagged_tuples = selection['atoms']
                    else:
                        flagged_tuples = [
                            (1, *info) for (_, *info) in selection['atoms']
                        ]
                else:
                    warning("There are no flags provided. Plot impossible.")
                    flagged_tuples = []

    return flagged_tuples


# -----------------------------------------------------------------------------
# --- Flag all bonds, angles, dihedrals ---------------------------------------
# -----------------------------------------------------------------------------


def flag_bats(mol: Mol, smarts: str = None, t_idxs: tuple = (), draw=False) -> dict:
    """
    Compute and flag bonds, angles, and dihedrals in a single call,
    automatically determining which interactions can be filtered
    based on the size of the SMARTS pattern and/or the number of
    atom indices supplied in t_idxs.

    Rules
    -----
    - If SMARTS has:
        2 atoms: only bonds can be filtered
        3 atoms: bonds + angles can be filtered
        >=4 atoms: bonds + angles + dihedrals can be filtered

    - If t_idxs has:
        len=2: only bonds
        len=3: bonds + angles
        len>=4: bonds + angles + dihedrals

    If both SMARTS and t_idxs are provided, the *maximal allowed degree*
    is the minimum of the two limits.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule under study.
    smarts : str, optional
        SMARTS pattern for filtering interactions.
    t_idxs : tuple[int], optional
        Atom indices for filtering interactions.
    draw: True or False flag for returning flagged
        rdkit.Chem.rdchem.Mol object

    Returns
    -------
    dict
        {
            'bonds':      [...],
            'angles':     [...],
            'dihedrals':  [...]
        }

    Interaction types that cannot be filtered due to SMARTS/t_idx size
    are returned fully active.

    if draw: rdkit.Chem.rdchem.Mol object of filtered features
    """

    # ------------------------------------------------------------
    # 1) Determine allowed level from SMARTS size
    # ------------------------------------------------------------
    smarts_level = 3  # bonds + angles + dihedrals

    if smarts:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            info(f"Invalid SMARTS '{smarts}' .\
                    Falling back to full reference molecule.")
        else:
            n_atoms = patt.GetNumAtoms()

            if n_atoms <= 1:
                info(
                    f"SMARTS '{smarts}' contains <2 atoms. No filtering "
                    f"is possible. All interactions returned active."
                )
                smarts_level = 0
                smarts = None
                t_idxs = ()

            elif n_atoms == 2:
                info(
                    f"SMARTS {smarts} contains 2 atoms. All angles and "
                    f"dihedrals will be flagged as inactive."
                )
                smarts_level = 1  # bonds

            elif n_atoms == 3:
                info(
                    f"SMARTS {smarts} contains 3 atoms. All dihedrals "
                    f"will be flagged as inactive."
                )
                smarts_level = 2  # bonds + angles

    # ------------------------------------------------------------
    # 2) Determine allowed level from t_idxs size
    # ------------------------------------------------------------
    t_level = 3

    if t_idxs:
        n = len(t_idxs)
        if n == 1:
            info(
                f"t_idxs has only 1 atom, filtering not possible. "
                f"All interactions returned active."
            )
            t_level = 0
            smarts = None
            t_idxs = ()

        elif n == 2:
            info(
                f"SMARTS {smarts} contains 2 atoms. All angles and "
                f"dihedrals will be flagged as inactive."
            )
            t_level = 1

        elif n == 3:
            info(
                f"SMARTS {smarts} contains 3 atoms. All dihedrals "
                f"will be flagged as inactive."
            )
            t_level = 2

    # ------------------------------------------------------------
    # 3) Effective filtering level
    # ------------------------------------------------------------
    # The most restrictive level applies
    filter_level = min(smarts_level, t_level)

    # ------------------------------------------------------------
    # 4) Call the individual flaggers
    # ------------------------------------------------------------
    out_atoms = flag_atoms(mol, smarts=smarts, t_idxs=t_idxs, draw=False)['atoms']
    out_bonds = flag_bonds(mol, smarts=smarts, t_idxs=t_idxs, draw=False)['bonds']
    out_angles = flag_angles(mol, smarts=smarts, t_idxs=t_idxs, draw=False)['angles']
    out_dihedrals = flag_dihedrals(mol, smarts=smarts, t_idxs=t_idxs, draw=False)[
        'dihedrals'
    ]

    # ------------------------------------------------------------
    # 5) Apply level truncation
    # ------------------------------------------------------------
    if filter_level == 0:
        info(
            f"SMARTS/t_idxs only point to a single atom; no filtering possible. "
            f"All bonds, angles, torsion will be returned (all active)."
        )

        d_flag = {
            'bonds': [(1, *info) for (_, *info) in out_bonds],  # force active
            'angles': [(1, *info) for (_, *info) in out_angles],  # force active
            'dihedrals': [(1, *info) for (_, *info) in out_dihedrals],  # force active
        }

    if filter_level == 1:
        d_flag = {
            'bonds': out_bonds,
            'angles': [(0, *ang) for (_, *ang) in out_angles],  # force inactive
            'dihedrals': [(0, *dih) for (_, *dih) in out_dihedrals],  # force inactive
        }

    if filter_level == 2:
        d_flag = {
            'bonds': out_bonds,
            'angles': out_angles,
            'dihedrals': [(0, *dih) for (_, *dih) in out_dihedrals],  # force inactive
        }

    else:
        d_flag = {'bonds': out_bonds, 'angles': out_angles, 'dihedrals': out_dihedrals}

    if draw:
        img = __get_highlight_molimg(mol, d_flag)
        return d_flag, img
    else:
        return d_flag, filter_level


def flag_bats_multiple(
    mol: Mol, l_smarts: list[str] = None, l_t_idxs: list[tuple] = (), draw=False
) -> dict:
    """
    Compute and flag bonds, angles, and dihedrals in a single call,
    ifor multiple structural features

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule under study.
    l_smarts : list[str], optional
        SMARTS patterns for filtering interactions.
    t_idxs : list[tuple[int]], optional
        list of atom indices tuples for filtering interactions.
    draw: True or False flag for returning flagged
        rdkit.Chem.rdchem.Mol object

    Returns
    -------
    dict
        {
            'bonds':      [...],
            'angles':     [...],
            'dihedrals':  [...]
        }

    Interaction types that cannot be filtered due to SMARTS/t_idx size
    are returned fully active.

    if draw: rdkit.Chem.rdchem.Mol object of filtered features
    """

    if l_smarts and not l_t_idxs:
        res = {}
        l_img = []
        l_levels = []
        for smarts in l_smarts:
            d_flag, flag_level = flag_bats(mol=mol, smarts=smarts, draw=False)
            res[smarts] = d_flag
            l_levels.append(flag_level)

        if draw:
            img = __get_img_multiple_mols(mol, res, l_smarts, l_levels)
            return res, img
        else:
            return res

    elif l_t_idxs and not l_smarts:
        res = {}
        res_img = []
        for t in l_t_idxs:
            d_flag = flag_bats(mol=mol, t_idxs=t, draw=False)
            res[t] = d_flag

        if draw:
            img = __get_img_multiple_mols(mol, res, l_t_idxs, l_levels)
            return res, img
        else:
            return res

    elif smarts and l_t_idxs:
        info(
            f"Info text with important data {var}"
        )  # where var is a variable you want printed
        warning(
            f"Either define a list of smarts or a list of tuples with atom indices!"
            f"Here, the list of smarts are ignored and only the"
            f"l_t_idxs: {l_t_idxs} are used for further processing."
        )

        res = {}
        res_img = []
        for t in l_t_idxs:
            d_flag = flag_bats(mol=mol, t_idxs=t, draw=False)
            res[t] = d_flag

        if draw:
            img = __get_img_multiple_mols(mol, res, l_t_idxs, l_levels)
            return res, img
        else:
            return res
