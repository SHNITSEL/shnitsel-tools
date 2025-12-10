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
    bond_list = [
        (1, (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) for bond in mol.GetBonds()
    ]

    return {'bonds': bond_list}


def __get_bonds_by_indices(match_indices_list: list[tuple], d_bonds: dict) -> dict:
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


def flag_bonds(mol: Mol, smarts: str = None, t_idxs: tuple = (), draw=True) -> dict:
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
            info(
                f"SMARTS pattern '{smarts}' was not found in the molecule. "
                f"Returning all bonds as active."
            )
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
            info(
                f"SMARTS pattern '{smarts}' was not found in the molecule. "
                f"Returning bonds for the provided target indices only."
            )
            t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
            d_flag = __get_bonds_by_indices(t_idxs_list, d_all_bonds)

        else:
            # flatten all matched indices into a set
            matched_atoms = set()
            for match in matches:
                matched_atoms.update(match)

            intersection = t_idxs_set & matched_atoms

            if not intersection:
                info(
                    f"There is no overlap between the SMARTS macthes {matches} "
                    f"and the target indices {t_idxs}. Returning bonds for "
                    f"the provided target indices only."
                )
                t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
                d_flag = __get_bonds_by_indices(t_idxs_list, d_all_bonds)

            else:
                info(
                    f"Partial overlap found between SMARTS matches {matches} and "
                    f"the target indices {t_idxs}. Only indices {sorted(intersection)} "
                    f"are considered for active bonds."
                )
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

