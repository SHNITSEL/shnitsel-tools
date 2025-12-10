import sys, os
import logging
from logging import warning, info

import pandas as pd
import numpy as np

import rdkit
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

st_yellow = (196/255, 160/255, 0/255)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    force=True
)


def __match_pattern(mol: Mol, smarts: str) -> list[tuple]:
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
    """
    pattern = Chem.MolFromSmarts(smarts)

    if pattern is None:
        info(f"Invalid SMARTS '{smarts}'. Falling back to full reference molecule.")
        matches = ()
    else:
        matches = mol.GetSubstructMatches(pattern)

    return matches

def __get_bond_info(
        mol: Mol,
        flagged_tuples
        ):
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
    flagged_tuples_binfo: list of tuples
        Each entry: (flag, atom_tuple, bond_tuple, bondtype_tuple, 
        rdkit.Chem.rdchem.Mol of submoli/subgraph )
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
        atom_pairs = [(t[i], t[i+1]) for i in range(len(t)-1)]
        
        inner_bond_idxs = []
        inner_bond_types = []
        for i,j in atom_pairs:
            bond = mol.GetBondBetweenAtoms(i,j)
            if bond:
                inner_bond_idxs.append(bond.GetIdx())
                inner_bond_types.append(bond.GetBondTypeAsDouble())

        l_bond_idxs.append(tuple(inner_bond_idxs))
        l_bond_types.append(tuple(inner_bond_types))
        l_submols.append(Chem.PathToSubmol(mol, inner_bond_idxs))
    
    flagged_tuples_binfo = list(zip(l_flags, l_atom_idxs, l_bond_idxs, l_bond_types, l_submols))

    return flagged_tuples_binfo


def __level_to_geo(flag_level):

    level_to_prop = {
            0: 'bonds',
            1: 'bonds',
            2: 'angles',
            3: 'dihedrals',
            4: 'dihedrals'
            }

    return level_to_prop[flag_level]

def __get_color_atoms(d_flag, flag_level):

    flagged_tuples = d_flag[__level_to_geo(flag_level)]
    color_atoms = [a for flag, at, *_ in flagged_tuples if flag==1 for a in at]

    return color_atoms

def __get_color_bonds(d_flag, flag_level):

    flagged_tuples = d_flag[__level_to_geo(flag_level)]
    color_bonds = [b for flag, at, bt, *_ in flagged_tuples if flag==1 for b in bt]

    return color_bonds

def __collect_tuples(d_flag):
    """
    Select the appropriate list of tuples based on hierarchy:
    dihedrals -> angles -> bonds.

    Each value is a list of tuples of the form:
        (flag, atom_tuple, bond_tuple, bondtype_tuple, Mol)
    """

    hierarchy = ["dihedrals", "angles", "bonds"]
    for key in hierarchy:
        if key in d_flag:
            tuples = d_flag[key]
            if any(entry[0] != 0 for entry in tuples):
                flagged_tuples = tuples
            else:
                warning("All flags in dihedrals, angles, and bonds are zero. "\
                        "Will collect all atoms instead.")
                if 'atoms' in d_flag:
                    if any(entry[0] != 0 for entry in tuples):
                        flagged_tuples = d_flag['atoms']
                    else:
                        flagged_tuples = [(1, *info) for (_, *info) in d_flag['atoms']]
                else:
                    warning("There are no flags provided. Plot impossible.")
                    flagged_tuples = []

    return flagged_tuples


def __get_highlight_molimg(
        mol: Mol, 
        d_flag,
        highlight_color=st_yellow,
        width=300,
        height=300
        ):
    """
    Convert a list of flagged atom tuples into bonds and atoms for highlighting.

    Parameters
    ----------
    mol : RDKit Mol
        Molecule object.
    d_flag : dictionary
        keys: 'atoms', 'bonds', 'angles' or 'dihedrals'
        values: list of tuples
        Each entry: (flag, atom_tuple, bond_tuple, bondtype_tuple, Mol)
        Example: [(1, (0,1), (0), (1.0)), (0, (1,2), (1), (2.0))]
    highlight_color : tuple, optional
        RGB tuple for highlighting bonds/atoms.

    Returns
    -------
    img : PIL.Image
        RDKit molecule image with highlighted atoms and bonds.
    """

    flagged_tuples = __collect_tuples(d_flag)

    # draw molecule with highlights
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().fillHighlights=True
    drawer.drawOptions().addAtomIndices = True
    drawer.drawOptions().setHighlightColour(highlight_color)
    drawer.drawOptions().clearBackground = False

    has_binfo = all(len(t) >= 5 for t in flagged_tuples)

    if has_binfo:
        highlights = [(at, bt) for flag, at, bt, bot, m in flagged_tuples if flag == 1]
        highlight_bonds = [bidx for at, bt in highlights for bidx in bt]
        highlight_atoms = [aidx for at, bt in highlights for aidx in at]
        
        rdMolDraw2D.PrepareAndDrawMolecule(
                drawer, 
                mol, 
                highlightAtoms=highlight_atoms, 
                highlightBonds=highlight_bonds)

    else:
        highlight_atoms = [idx for flag, (idx, symbol) in flagged_tuples if flag == 1]

        rdMolDraw2D.PrepareAndDrawMolecule(
                drawer, 
                mol, 
                highlightAtoms=highlight_atoms
                ) 

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img_text = drawer.GetDrawingText()
    img = SVG(img_text)

    return img


# -------------------------------------------------------------------
# -- Atom specific functions ----------------------------------------
# -------------------------------------------------------------------

def __get_all_atoms(mol: Mol) -> dict:
    """
    Return all atoms in the molecule as a dictionary with flags.

    All atoms are initially flagged as active (1).

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.

    Returns
    -------
    dict
        Dictionary with key 'atoms' mapping to a list of tuples.
        Each tuple has the form `(flag, (atom_idx, atom_symbol))`, where:
            - flag : int
                1 indicates active atom.
            - (atom_idx, atom_symbol) : tuple
                The atom index in the molecule and its atomic symbol
                (e.g., (0, 'C'), (1, 'H')).
    """
    atom_list = [
        (1, (atom.GetIdx(), atom.GetSymbol()))
        for atom in mol.GetAtoms()
    ]

    return {'atoms': atom_list}


def __get_atoms_by_indices(
    match_indices_list: list[tuple],
    d_atoms: dict) -> dict:
    """
    Set flags for atoms based on whether they are contained within substructure matches.

    Atoms are active (1) if they represent also nodes in the matched substructures.
    All other atoms are flagged as inactive (0).

    Parameters
    ----------
    match_indices_list : list of tuples
        Each tuple contains atom indices (e.g. obtained form substructure match).
    d_bonds : dict
        Dictionary of all bonds in the molecule (output of `__get_all_bonds`).

    Returns
    -------
    dict
        Updated atoms dictionary with flags updated based on the substructure matches.
    """
    updated_atoms = []
    for flag, (idx, symbol) in d_atoms['atoms']:
        atom_active = any(idx in match for match in match_indices_list)
        updated_atoms.append((1 if atom_active else 0, (idx, symbol)))

    return {'atoms': updated_atoms}


def flag_atoms(
    mol: Mol,
    smarts: str = None,
    t_idxs: tuple = (),
    draw=False) -> dict:
    """
    Flag atoms in a molecule based on substructure patterns or atom indices.

    Atoms can be flagged in four ways:
        1. No SMARTS or target indices provided: all atoms are returned as active.
        2. SMARTS provided: atoms belonging to the substructure matches are active; others are inactive.
        3. Target atom indices provided: only atoms in the provided tuple are active; others are inactive.
        4. Both SMARTS and target indices provided: only the intersection of SMARTS matches and t_idxs are considered active.
           Warnings are issued if there is no overlap or partial overlap.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.
    smarts : str, optional
        SMARTS pattern to filter atoms. Default is None.
    t_idxs : tuple, optional
        Tuple of atom indices to filter atoms. Default is empty tuple.
    draw : bool, optional
        If True, returns an image highlighting the active atoms. Default is True.

    Returns
    -------
    dict or tuple
        If draw=False: dictionary with key 'atoms' mapping to a list of tuples `(flag, (atom_idx, atom_symbol))`.
        If draw=True: tuple of (atom dictionary, highlighted molecule image).
    """

    d_all_atoms = __get_all_atoms(mol)

    # Case 1: no filtering
    if not smarts and not t_idxs:
        d_flag = d_all_atoms

    # Case 2: SMARTS filtering
    elif smarts and not t_idxs:
        matches = __match_pattern(mol, smarts)

        if not matches:
            info(f"SMARTS pattern '{smarts}' was not found in the molecule. "
                 "Returning all atoms as active.")
            d_flag = d_all_atoms
        else:
            d_flag = __get_atoms_by_indices(matches, d_all_atoms)

    # Case 3: target indices filtering
    elif not smarts and t_idxs:
        t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
        d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)

    # Case 4: both SMARTS and target indices
    elif smarts and t_idxs:
        matches = __match_pattern(mol, smarts)
        t_idxs_set = set(t_idxs)

        if not matches:
            info(f"SMARTS pattern '{smarts}' was not found in the molecule. "
                 "Returning atoms for the provided target indices only.")
            t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
            d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)
        else:
            # flatten all matched indices into a set
            matched_atoms = set()
            for match in matches:
                matched_atoms.update(match)

            intersection = t_idxs_set & matched_atoms

            if not intersection:
                info(f"No overlap between SMARTS matches {matches} and "
                     f"target indices {t_idxs}. Returning atoms for the target indices only.")
                t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
                d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)
            else:
                info(f"Partial overlap found between SMARTS matches {matches} and "
                     f"target indices {t_idxs}. Only indices {sorted(intersection)} "
                     "are considered active.")
                t_idxs_list = [tuple(intersection)]
                d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)

    if draw:
        img = __get_highlight_molimg(mol, d_flag)
        return d_flag, img
    else:
        return d_flag 


# -------------------------------------------------------------------
# -- Bond specific functions ----------------------------------------
# -------------------------------------------------------------------

def __get_all_bonds(mol: Mol) -> dict:
    """
    Return all bonds in the molecule as a dictionary with flags.

    All bonds are initially flagged as active (1).

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.

    Returns
    -------
    dict
        Dictionary with key 'bonds' mapping to a list of tuples.
        Each tuple has the form `(flag, (atom_idx1, atom_idx2))`, where:
            - flag : int
                0 indicates active bond.
            - (atom_idx1, atom_idx2) : tuple
                Pair of atom indices defining the bond.
    """
    bond_list = [(1, (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                 for bond in mol.GetBonds()]

    return {'bonds': bond_list}


def __get_bonds_by_indices(
    match_indices_list: list[tuple],
    d_bonds: dict) -> dict:
    """
    Set flags for bonds based on whether they are contained within substructure matches.

    Bonds are active (1) if both atoms belong to any of the matched substructures.
    All other bonds are flagged as inactive (0).

    Parameters
    ----------
    match_indices_list : list of tuples
        Each tuple contains atom indices (e.g. obtained form substructure match).
    d_bonds : dict
        Dictionary of all bonds in the molecule (output of `__get_all_bonds`).

    Returns
    -------
    dict
        Updated bond dictionary with the same structure as `d_bonds`,
        but flags updated based on the substructure matches.
    """
    updated_bonds = []
    for flag, (i, j) in d_bonds['bonds']:
        bond_active = any(i in match and j in match for match in match_indices_list)
        updated_bonds.append((1 if bond_active else 0, (i, j)))

    return {'bonds': updated_bonds}


def flag_bonds(
    mol: Mol,
    smarts: str = None,
    t_idxs: tuple = (),
    draw=True) -> dict:
    """
    Flag bonds in a molecule based on substructure patterns or atom indices.

    Bonds can be flagged in three ways:
        1. No SMARTS or target indices provided: all bonds are returned as active.
        2. SMARTS provided: bonds belonging to the substructure matches are active; others are inactive.
        3. Target atom indices provided: bonds entirely contained in the atom index tuple are active; others are inactive.
        4. Both SMARTS and target indices provided: only the intersection of SMARTS matches and t_idxs are considered active.
           Warnings are issued if there is no overlap or partial overlap.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object.
    smarts : str, optional
        SMARTS pattern to filter bonds. Default is None.
    t_idxs : tuple, optional
        Tuple of atom indices to filter bonds. Default is empty tuple.

    Returns
    -------
    dict
        Dictionary with key 'bonds' mapping to a list of tuples `(flag, (atom_idx1, atom_idx2))`.
    """

    out_atoms = flag_atoms(mol=mol, smarts=smarts, t_idxs=t_idxs, draw=False)['atoms']

    d_all_bonds = __get_all_bonds(mol)

    # case 1: no filtering
    if not smarts and not t_idxs:
        d_flag = d_all_bonds

    # case 2: SMARTS filtering
    elif smarts and not t_idxs:
        matches = __match_pattern(mol, smarts)

        if not matches:
            info(f"SMARTS pattern '{smarts}' was not found in the molecule. "\
                    f"Returning all bonds as active.")
            d_flag = d_all_bonds

        else:
            d_flag = __get_bonds_by_indices(matches, d_all_bonds)

    # case 3: target indices filtering
    elif not smarts and t_idxs:
        if isinstance(t_idxs, tuple):
            t_idxs_list = [t_idxs]
        else:
            t_idxs_list = list(t_idxs)
        d_flag = __get_bonds_by_indices(t_idxs_list, d_all_bonds)

    # case 4: both SMARTS and target indices
    elif smarts and t_idxs:
        matches = __match_pattern(mol, smarts)
        t_idxs_set = set(t_idxs)

        if not matches:
            info(f"SMARTS pattern '{smarts}' was not found in the molecule. "\
                    f"Returning bonds for the provided target indices only.")
            t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
            d_flag = __get_bonds_by_indices(t_idxs_list, d_all_bonds)

        else:
            # flatten all matched indices into a set
            matched_atoms = set()
            for match in matches:
                matched_atoms.update(match)

            intersection = t_idxs_set & matched_atoms

            if not intersection:
                info(f"There is no overlap between the SMARTS macthes {matches} "\
                        f"and the target indices {t_idxs}. Returning bonds for "\
                        f"the provided target indices only.")
                t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
                d_flag = __get_bonds_by_indices(t_idxs_list, d_all_bonds)

            else:
                info(f"Partial overlap found between SMARTS matches {matches} and "\
                    f"the target indices {t_idxs}. Only indices {sorted(intersection)} "\
                    f"are considered for active bonds.")
                t_idxs_list = [tuple(intersection)]
                d_flag = __get_bonds_by_indices(t_idxs_list, d_all_bonds)

    d_flag_binfo = {}
    d_flag_binfo['atoms'] = out_atoms
    d_flag_binfo['bonds'] = __get_bond_info(mol, d_flag['bonds'])

    if draw:
        img = __get_highlight_molimg(mol, d_flag_binfo) 
        return d_flag_binfo, img
    else:
        return d_flag_binfo


# -------------------------------------------------------------------
# -- Angle specific functions ---------------------------------------
# -------------------------------------------------------------------

def __get_all_angles(mol: Mol) -> dict:
    """
    Return all angles in the molecule as a dictionary with flags.

    Angles are defined as atom triples (i, j, k) where
    i–j and j–k are both bonded.

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
        active = any(
            (i in match and j in match and k in match)
            for match in match_list
        )
        updated.append((1 if active else 0, (i, j, k)))

    return {'angles': updated}


def flag_angles(
    mol: Mol,
    smarts: str = None,
    t_idxs: tuple = (),
    draw=False) -> dict:
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
            info(f"SMARTS pattern '{smarts}' was not found. "\
                    f"Returning all angles as active.")
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
            info(f"SMARTS '{smarts}' not found. "\
                    f"Returning angles for provided atom indices only.")
            d_flag = __get_angles_by_indices([t_idxs], d_all_angles)

        # Collect all atoms from matches
        matched_atoms = set()
        for match in matches:
            matched_atoms.update(match)

        # Intersection
        inter = matched_atoms & t_set

        if not inter:
            info("No overlap between SMARTS matches and target indices. "\
                    f"Returning angles for t_idxs only.")
            d_flag = __get_angles_by_indices([t_idxs], d_all_angles)
        else:
            info(f"Partial overlap between SMARTS and t_idxs. "\
                    f"Using only atoms {sorted(inter)}.")
            d_flag = __get_angles_by_indices([tuple(inter)], d_all_angles)

    d_flag_binfo = {}
    d_flag_binfo['atoms'] = out_atoms
    d_flag_binfo['angles'] = __get_bond_info(mol, d_flag['angles'])

    if draw:
        img = __get_highlight_molimg(mol, d_flag_binfo) 
        return d_flag_binfo, img
    else:
        return d_flag_binfo


# -------------------------------------------------------------------
# -- Torsion specific functions -------------------------------------
# -------------------------------------------------------------------

def __get_all_dihedrals(mol: Mol) -> dict:
    """
    Return all dihedrals in the molecule as a dictionary with flags.

    Dihedral quadruples are of the form (i, j, k, l) where:
        i–j, j–k, and k–l are all bonds.

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
        mol: Mol, 
        smarts: str = None, 
        t_idxs: tuple = (),
        draw=False) -> dict:
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
            info(f"SMARTS pattern '{smarts}' not found. "\
                    "Returning all dihedrals active.")
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
            info(f"No overlap between SMARTS match and t_idxs. "\
                    f"Returning dihedrals for atom indices {t_idxs} only.")

            d_flag =  __get_dihedrals_by_indices([t_idxs], d_all)

        else:
            info(f"Partial overlap between SMARTS and t_idxs. "\
                    f"Using atoms {sorted(inter)} for dihedral filtering.")

            d_flag = __get_dihedrals_by_indices([tuple(inter)], d_all)

    d_flag_binfo = {}
    d_flag_binfo['atoms'] = out_atoms
    d_flag_binfo['dihedrals'] = __get_bond_info(mol, d_flag['dihedrals'])

    if draw:
        img = __get_highlight_molimg(mol, d_flag_binfo)
        return d_flag_binfo, img
    else:
        return d_flag_binfo

# -----------------------------------------------------------------------------
# --- Flag all bonds, angles, dihedrals ---------------------------------------
# -----------------------------------------------------------------------------

def flag_bats(
    mol: Mol,
    smarts: str = None,
    t_idxs: tuple = (),
    draw=False) -> dict:
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
    smarts_level = 3 # bonds + angles + dihedrals

    if smarts:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            info(f"Invalid SMARTS '{smarts}' .\
                    Falling back to full reference molecule.")
        else:
            n_atoms = patt.GetNumAtoms()

            if n_atoms <= 1:
                info(f"SMARTS '{smarts}' contains <2 atoms. No filtering "\
                        f"is possible. All interactions returned active.")
                smarts_level = 0
                smarts = None
                t_idxs = ()

            elif n_atoms == 2:
                info(f"SMARTS {smarts} contains 2 atoms. All angles and "\
                        f"dihedrals will be flagged as inactive.")
                smarts_level = 1 # bonds

            elif n_atoms == 3:
                info(f"SMARTS {smarts} contains 3 atoms. All dihedrals "\
                        f"will be flagged as inactive.")
                smarts_level = 2 # bonds + angles

    # ------------------------------------------------------------
    # 2) Determine allowed level from t_idxs size
    # ------------------------------------------------------------
    t_level = 3

    if t_idxs:
        n = len(t_idxs)
        if n == 1:
            info(f"t_idxs has only 1 atom, filtering not possible. "\
                    f"All interactions returned active.")
            t_level = 0
            smarts = None
            t_idxs = ()

        elif n == 2:
            info(f"SMARTS {smarts} contains 2 atoms. All angles and "\
                    f"dihedrals will be flagged as inactive.")
            t_level = 1

        elif n == 3:
            info(f"SMARTS {smarts} contains 3 atoms. All dihedrals "\
                    f"will be flagged as inactive.")
            t_level = 2

    # ------------------------------------------------------------
    # 3) Effective filtering level
    # ------------------------------------------------------------
    # The most restrictive level applies
    filter_level = min(smarts_level, t_level)

    # ------------------------------------------------------------
    # 4) Call the individual flaggers
    # ------------------------------------------------------------
    out_atoms     = flag_atoms(mol, smarts=smarts, t_idxs=t_idxs, draw=False)['atoms']
    out_bonds     = flag_bonds(mol, smarts=smarts, t_idxs=t_idxs, draw=False)['bonds']
    out_angles    = flag_angles(mol, smarts=smarts, t_idxs=t_idxs, draw=False)['angles']
    out_dihedrals = flag_dihedrals(mol, smarts=smarts, t_idxs=t_idxs, draw=False)['dihedrals']

    # ------------------------------------------------------------
    # 5) Apply level truncation
    # ------------------------------------------------------------
    if filter_level == 0:
        info(f"SMARTS/t_idxs only point to a single atom; no filtering possible. "\
                f"All bonds, angles, torsion will be returned (all active).")

        d_flag = {
            'bonds':     [(1, *info) for (_, *info) in out_bonds],  # force active
            'angles':    [(1, *info) for (_, *info) in out_angles],   # force active
            'dihedrals': [(1, *info) for (_, *info) in out_dihedrals] # force active
        }

    if filter_level == 1:
        d_flag = {
            'bonds':     out_bonds,
            'angles':    [(0, *ang) for (_, *ang) in out_angles],   # force inactive
            'dihedrals': [(0, *dih) for (_, *dih) in out_dihedrals] # force inactive
        }

    if filter_level == 2:
        d_flag = {
            'bonds':     out_bonds,
            'angles':    out_angles,
            'dihedrals': [(0, *dih) for (_, *dih) in out_dihedrals] # force inactive
        }

    else:
        d_flag = {
            'bonds':     out_bonds,
            'angles':    out_angles,
            'dihedrals': out_dihedrals
        }

    if draw:
        img = __get_highlight_molimg(mol, d_flag)
        return d_flag, img
    else:
        return d_flag, filter_level


def flag_bats_multiple(
    mol: Mol,
    l_smarts: list[str] = None,
    l_t_idxs: list[tuple] = (),
    draw=False) -> dict:
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
        info(f"Info text with important data {var}") # where var is a variable you want printed
        warning(f"Either define a list of smarts or a list of tuples with atom indices!"
                f"Here, the list of smarts are ignored and only the"
                f"l_t_idxs: {l_t_idxs} are used for further processing.")
        
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


def __get_img_multiple_mols(
        mol: Mol, 
        d_multi_flag: dict, 
        l_patterns: list, 
        l_levels: list) -> SVG:
    """
    Generate a single SVG image containing multiple copies of a molecule,
    each highlighted according to supplied atom/bond pattern levels.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object. The same molecule is drawn once for each
        entry in `l_patterns` / `l_levels`.
    d_multi_flag : dict
        Dictionary with mapping information on bonds, angles torsion
        Each entry: {'bonds': [(0, (5, 0, 1), (6, 0), (1.0, 2.0), rdkit.Chem.rdchem.Mol object)]}

        Each entry contains information on the atom indexes (2nd) and bond indexes (3rd)
        element needed by the helper functions `__get_color_atoms()` and `__get_color_bonds()` 
        to extract the atoms and bonds to be highlighted.
    l_patterns : Iterable
        A list of patterns (smarts or tuples of atom indices) indexing into `d_multi_flag`.  
    l_levels : Iterable
        List of highlight “levels” corresponding to `l_patterns`.
        Each level is passed to the highlight extraction helpers to control
        the highlight color intensity or style.

    Returns
    -------
    PIL.Image.Image or IPython.display.SVG
        A grid image produced by `rdkit.Chem.Draw.MolsToGridImage`, 
        containing all molecule renderings arranged in a single row. 
        Each copy of the molecule is highlighted with its own atom and bond sets determined 
        by the input patterns (smarts or atom indices).
    """

    highlight_atoms = []
    highlight_bonds = []
    mols = []
    for i, level in zip(l_patterns, l_levels):
        mols.append(mol)
        highlight_atoms.append(__get_color_atoms(d_multi_flag[i], level))
        highlight_bonds.append(__get_color_bonds(d_multi_flag[i], level))

    opts = Draw.MolDrawOptions()
    opts.addAtomIndices = True
    opts.fillHighlights = True
    opts.setHighlightColour(st_yellow)

    img = Draw.MolsToGridImage(
            mols,
            highlightAtomLists=highlight_atoms,
            highlightBondLists=highlight_bonds,
            molsPerRow=len(mols),
            subImgSize=(300, 300),
            useSVG=True,
            drawOptions=opts
            )

    return img

# -----------------------------------------------------------------------------
# --- BLA ---------------------------------------------------------------------
# -----------------------------------------------------------------------------

def __build_conjugated_smarts(n_double: int,
                             elems: str = "#6,#7,#8,#15,#16") -> str:
    """
    Build a SMARTS pattern for a linear conjugated system with `n_double`
    alternating double bonds.

    Example (n_double=2):
    [#6,#7]=[#6,#7]-[#6,#7]=[#6,#7]

    Parameters
    ----------
    n_double : int
        Number of C=C-like double bonds.
    elems : str
        SMARTS atomic specification (default: C,N,O,P,S).

    Returns
    -------
    str
        SMARTS string encoding the conjugated system.
    """

    if n_double < 1:
        raise ValueError("n_double must be >= 1")

    unit = f"[{elems}]=[{elems}]"
    return "-".join([unit] * n_double)


def __match_bla_chromophor(mol,
                          smarts: str | None = None,
                          n_double: int | None = None,
                          elems: str = "#6,#7,#8,#15,#16"):
    """
    Detect conjugated chromophores defined either by SMARTS, number of
    alternating double bonds, or automatically by maximum extension.

    Decision logic
    --------------
    1. If `smarts` is given → use SMARTS
    2. If both `smarts` and `n_double` are given → validate consistency
    3. If only `n_double` is given → generate SMARTS
    4. If neither is given → search for maximum conjugated chromophore

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule to search.
    smarts : str, optional
        Explicit SMARTS pattern defining the chromophore.
    n_double : int, optional
        Number of alternating double bonds.
    elems : str
        Allowed atoms in conjugation.

    Returns
    -------
    dict
        {
            "smarts": SMARTS used,
            "n_double": number of double bonds,
            "matches": substructure matches
        }
    """

    take = None
    info_msg = None

    # --- case 1: SMARTS only ---
    if smarts and not n_double:
        take = "smarts"

    # --- case 2: SMARTS and n_double ---
    elif smarts and n_double:
        patt = Chem.MolFromSmarts(smarts)
        dbl = sum(
            1 for b in patt.GetBonds()
            if b.GetBondType() == Chem.rdchem.BondType.DOUBLE
        )

        if dbl >= n_double:
            take = "smarts"
        else:
            take = "n_double"

        info(f"SMARTS encodes {dbl} double bonds, "\
             f"n_double={n_double}. Using {take} for pattern matching.")

    # --- case 3: n_double only ---
    elif n_double and not smarts:
        print('entered here..')
        take = "n_double"

    # --- case 4: neither given: auto-discovery ---
    else:
        max_matches = ()
        best_smarts = None
        best_n = 0

        n = 2
        while True:
            smi = __build_conjugated_smarts(n, elems)
            matches = __match_pattern(mol, smi)
            if not matches:
                break
            best_smarts = smi
            best_n = n
            max_matches = matches
            n += 1

        # overwrite smarts with best smarts
        smarts = best_smarts

    # --- execute selected mode ---
    if take == "n_double":
        smarts = __build_conjugated_smarts(n_double, elems)
        print(smarts)

    matches = __match_pattern(mol, smarts)
    if not matches:
        info(f"Given {take} parameter doesn't match the structure.")

    return {
        "smarts": smarts,
        "n_double": n_double,
        "matches": matches
    }


def flag_bla_chromophor(mol,
                        smarts: str | None = None,
                        n_double: int | None = None,
                        elems: str = "#6,#7,#8,#15,#16",
                        draw=True,
                        width: int = 500,
                        height: int = 300):

    match_chromo = __match_bla_chromophor(mol, smarts, n_double, elems)

    if not match_chromo['matches']:
        info(f"SMARTS ({match_chromo['smarts']}) and/or n_double "\
             f"({match_chromo['n_double']}) not found in the molecule. "\
             f"No bonds will be flagged.")
        d_flag = {}
        return d_flag

    else:
        out_atoms = flag_atoms(mol=mol, 
                smarts=match_chromo['smarts'], 
                t_idxs=None, 
                draw=False)['atoms']

        d_all_bonds = __get_all_bonds(mol)
        d_flag = __get_bonds_by_indices(match_chromo['matches'], d_all_bonds)

        d_flag_binfo = {}
        d_flag_binfo['atoms'] = out_atoms
        d_flag_binfo['bonds'] = __get_bond_info(mol, d_flag['bonds'])

        if draw:
            img = __get_highlight_molimg(mol, d_flag_binfo, width=width, height=height)
            return d_flag_binfo, img
        else:
            return d_flag_binfo


