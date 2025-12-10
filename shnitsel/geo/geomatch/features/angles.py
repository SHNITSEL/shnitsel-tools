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
# -- Angle specific functions ---------------------------------------
# -------------------------------------------------------------------


def __get_all_angles(mol: Mol) -> dict:
    """
    Return all angles in the molecule as a dictionary with flags.

    Angles are defined as atom triples (i, j, k) where
    i-j and j-k are both bonded.

    All angles are initially flagged as active (1).

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.

    Returns
    -------
    dict
        Dictionary with key 'angles' mapping to a list of tuples:
        (flag, (i, j, k))
    """
    angles = []
    for bond_j in mol.GetBonds():
        j = bond_j.GetBeginAtomIdx()
        k = bond_j.GetEndAtomIdx()

        # find atoms bonded to j (except k)
        neighbors_j = [
            nbr.GetIdx()
            for nbr in mol.GetAtomWithIdx(j).GetNeighbors()
            if nbr.GetIdx() != k
        ]

        # angles are (i, j, k)
        for i in neighbors_j:
            angles.append((1, (i, j, k)))

        # also angles (k, j, i) by symmetry
        neighbors_k = [
            nbr.GetIdx()
            for nbr in mol.GetAtomWithIdx(k).GetNeighbors()
            if nbr.GetIdx() != j
        ]
        for i in neighbors_k:
            angles.append((1, (i, k, j)))

    return {'angles': angles}


def __get_angles_by_indices(match_list: list[tuple], d_angles: dict) -> dict:
    """
    Flag angles as active (1) or inactive (0) based on substructure matches.

    An angle (i, j, k) is active if all three indices belong to
    the *same* match tuple.

    Parameters
    ----------
    match_list : list of tuples
        Each tuple contains atom indices of a SMARTS match.
    d_angles : dict
        Angle dictionary from __get_all_angles().

    Returns
    -------
    dict
        Updated angle dictionary with active/inactive flags.
    """
    updated = []

    for flag, (i, j, k) in d_angles['angles']:
        active = any((i in match and j in match and k in match) for match in match_list)
        updated.append((1 if active else 0, (i, j, k)))

    return {'angles': updated}


def flag_angles(mol: Mol, smarts: str = None, t_idxs: tuple = (), draw=False) -> dict:
    """
    Flag molecule angles based on SMARTS patterns and/or atom indices.

    Modes of operation
    ------------------
    1) No SMARTS and no t_idxs: return all angles as active.
    2) SMARTS only: angles part of any SMARTS match are active.
    3) t_idxs only: only angles fully inside t_idxs are active.
    4) SMARTS + t_idxs:
          - If SMARTS fails: warn, return angles from t_idxs only.
          - If SMARTS matches but no overlap: warn, return angles from t_idxs.
          - If partial overlap: warn, return only angles in the intersection.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.
    smarts : str, optional
        SMARTS pattern representing the angle substructure.
    t_idxs : tuple, optional
        Atom indices to filter angles by.

    Returns
    -------
    dict
        Dictionary with key 'angles' mapping to list of (flag, (i, j, k)).
        Active = 1, inactive = 0.
    """

    out_atoms = flag_atoms(mol=mol, smarts=smarts, t_idxs=t_idxs, draw=False)['atoms']

    d_all_angles = __get_all_angles(mol)

    # ---- CASE 1: no filters ----
    if not smarts and not t_idxs:
        d_flag = d_all_angles

    # ---- CASE 2: SMARTS only ----
    if smarts and not t_idxs:
        matches = __match_pattern(mol, smarts)
        if not matches:
            info(
                f"SMARTS pattern '{smarts}' was not found. "
                f"Returning all angles as active."
            )
            d_flag = d_all_angles
        else:
            d_flag = __get_angles_by_indices(matches, d_all_angles)

    # ---- CASE 3: t_idxs only ----
    if not smarts and t_idxs:
        match_list = [t_idxs]
        d_flag = __get_angles_by_indices(match_list, d_all_angles)

    # ---- CASE 4: SMARTS + t_idxs ----
    if smarts and t_idxs:
        matches = __match_pattern(mol, smarts)
        t_set = set(t_idxs)

        # SMARTS fails: fallback
        if not matches:
            info(
                f"SMARTS '{smarts}' not found. "
                f"Returning angles for provided atom indices only."
            )
            d_flag = __get_angles_by_indices([t_idxs], d_all_angles)

        # Collect all atoms from matches
        matched_atoms = set()
        for match in matches:
            matched_atoms.update(match)

        # Intersection
        inter = matched_atoms & t_set

        if not inter:
            info(
                "No overlap between SMARTS matches and target indices. "
                f"Returning angles for {t_idxs} only."
            )
            d_flag = __get_angles_by_indices([t_idxs], d_all_angles)
        else:
            info(
                f"Partial overlap between SMARTS and t_idxs. "
                f"Using only atoms {sorted(inter)}."
            )
            d_flag = __get_angles_by_indices([tuple(inter)], d_all_angles)

    d_flag_binfo = {}
    d_flag_binfo['atoms'] = out_atoms
    d_flag_binfo['angles'] = __get_bond_info(mol, d_flag['angles'])

    if draw:
        img = __get_highlight_molimg(mol, d_flag_binfo)
        return d_flag_binfo, img
    else:
        return d_flag_binfo
