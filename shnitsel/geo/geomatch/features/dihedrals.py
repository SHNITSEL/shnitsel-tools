import logging
from logging import warning, info

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

from shnitsel.filtering.structure_selection import StructureSelection
from .atoms import flag_atoms
from ..flagging import __get_bond_info, __match_pattern
from ..presentation import __get_highlight_molimg
from shnitsel.vis.colormaps import st_yellow


# -------------------------------------------------------------------
# -- Torsion specific functions -------------------------------------
# -------------------------------------------------------------------


def __get_all_dihedrals(mol: Mol) -> dict:
    """
    Return all dihedrals in the molecule as a dictionary with flags.

    Dihedral quadruples are of the form (i, j, k, l) where:
        i-j, j-k, and k-l are all bonds.

    All dihedrals are initially flagged as active (1).

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.

    Returns
    -------
    dict
        Dictionary with key 'dihedrals' mapping to a list of:
        (flag, (i, j, k, l))
    """
    dihedrals = []

    for bond_jk in mol.GetBonds():
        j = bond_jk.GetBeginAtomIdx()
        k = bond_jk.GetEndAtomIdx()

        # atoms bonded to j (except k)
        neighbors_j = [
            nbr.GetIdx()
            for nbr in mol.GetAtomWithIdx(j).GetNeighbors()
            if nbr.GetIdx() != k
        ]

        # atoms bonded to k (except j)
        neighbors_k = [
            nbr.GetIdx()
            for nbr in mol.GetAtomWithIdx(k).GetNeighbors()
            if nbr.GetIdx() != j
        ]

        # form dihedrals (i, j, k, l)
        for i in neighbors_j:
            for l in neighbors_k:
                dihedrals.append((1, (i, j, k, l)))

        # also handle reversed central bond direction (k–j)
        # giving quadruples (i, k, j, l)
        for i in neighbors_k:
            for l in neighbors_j:
                dihedrals.append((1, (i, k, j, l)))

    return {'dihedrals': dihedrals}


def __get_dihedrals_by_indices(match_list: list[tuple], d_dihedrals: dict) -> dict:
    """
    Flag dihedrals as active (1) if all four atoms (i, j, k, l)
    belong to the same SMARTS match.

    Parameters
    ----------
    match_list : list of tuples
        SMARTS matches from __match_pattern.
    d_dihedrals : dict
        Output of __get_all_dihedrals().

    Returns
    -------
    dict
        Updated dihedral dictionary with flags.
    """
    updated = []

    for flag, (i, j, k, l) in d_dihedrals['dihedrals']:
        active = any(
            (i in match and j in match and k in match and l in match)
            for match in match_list
        )
        updated.append((1 if active else 0, (i, j, k, l)))

    return {'dihedrals': updated}


def flag_dihedrals(
    mol: Mol, smarts: str = None, t_idxs: tuple = (), draw=False
) -> dict:
    """
    Flag dihedrals in a molecule based on SMARTS and/or atom indices.

    Modes
    -----
    1) No SMARTS + no t_idxs: return all dihedrals active
    2) SMARTS only: dihedrals part of SMARTS matches are active
    3) t_idxs only: dihedrals fully inside t_idxs are active
    4) SMARTS + t_idxs: Find intersection behavior:
            - No SMARTS match: return t_idxs only.
            - No overlap: return t_idxs only.
            - Overlap: return only intersecting dihedrals.

    Parameters
    ----------
    mol : RDKit Mol
        Molecule under study.
    smarts : str, optional
        SMARTS pattern.
    t_idxs : tuple, optional
        Atom index tuple for filtering.

    Returns
    -------
    dict
        {'dihedrals': [(flag, (i,j,k,l)), ...]}
    """

    out_atoms = flag_atoms(mol=mol, smarts=smarts, t_idxs=t_idxs, draw=False)['atoms']
    d_all = __get_all_dihedrals(mol)

    # CASE 1: no filtering
    if not smarts and not t_idxs:
        d_flag = d_all

    # CASE 2: SMARTS only
    if smarts and not t_idxs:
        matches = __match_pattern(mol, smarts)
        if not matches:
            info(
                f"SMARTS pattern '{smarts}' not found. "
                "Returning all dihedrals active."
            )
            d_flag = d_all
        else:
            d_flag = __get_dihedrals_by_indices(matches, d_all)

    # CASE 3: t_idxs only
    if not smarts and t_idxs:
        d_flag = __get_dihedrals_by_indices([t_idxs], d_all)

    # CASE 4: SMARTS and t_idxs
    if smarts and t_idxs:
        matches = __match_pattern(mol, smarts)
        t_set = set(t_idxs)

        if not matches:
            info(f"SMARTS '{smarts}' not found. Returning dihedrals for t_idxs only.")

            d_flag = __get_dihedrals_by_indices([t_idxs], d_all)

        # flatten match atoms
        matched_atoms = set().union(*matches)
        inter = matched_atoms & t_set

        if not inter:
            info(
                f"No overlap between SMARTS match and t_idxs. "
                f"Returning dihedrals for atom indices {t_idxs} only."
            )

            d_flag = __get_dihedrals_by_indices([t_idxs], d_all)

        else:
            info(
                f"Partial overlap between SMARTS and t_idxs. "
                f"Using atoms {sorted(inter)} for dihedral filtering."
            )

            d_flag = __get_dihedrals_by_indices([tuple(inter)], d_all)

    d_flag_binfo = {}
    d_flag_binfo['atoms'] = out_atoms
    d_flag_binfo['dihedrals'] = __get_bond_info(mol, d_flag['dihedrals'])

    if draw:
        img = __get_highlight_molimg(mol, d_flag_binfo)
        return d_flag_binfo, img
    else:
        return d_flag_binfo
