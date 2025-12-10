import logging
from logging import warning, info

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

from shnitsel.filtering.structure_selection import StructureSelection
from shnitsel.geo.geomatch.flagging import __match_pattern
from shnitsel.geo.geomatch.presentation import __get_highlight_molimg
from shnitsel.vis.colormaps import st_yellow

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
    atom_list = [(1, (atom.GetIdx(), atom.GetSymbol())) for atom in mol.GetAtoms()]

    return {'atoms': atom_list}


def __get_atoms_by_indices(match_indices_list: list[tuple], d_atoms: dict) -> dict:
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


def flag_atoms(mol: Mol, smarts: str = None, t_idxs: tuple = (), draw=False) -> dict:
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

    # TODO: FIXME: Currently broken and will not load due to syntax issues

    d_all_atoms = __get_all_atoms(mol)

    # Case 1: no filtering
    if not smarts and not t_idxs:
        d_flag = d_all_atoms

    # Case 2: SMARTS filtering
    elif smarts and not t_idxs:
        matches = __match_pattern(mol, smarts)

        if not matches:
            info(
                f"SMARTS pattern '{smarts}' was not found in the molecule. "
                "Returning all atoms as active."
            )
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
            info(
                f"SMARTS pattern '{smarts}' was not found in the molecule. "
                "Returning atoms for the provided target indices only."
            )
            t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
            d_flag = __get_atomed indices into a set
            matched_atoms = set()
            for match in matches:
                matched_atoms.update(match)

            intersection = t_idxs_set & matched_atoms

            if not intersection:
                info(
                    f"No overlap between SMARTS matches {matches} and "
                    f"target indices {t_idxs}. Returning atoms for the target indices only."
                )
                t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
                d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)
            else:
                info(
                    f"Partial overlap found between SMARTS matches {matches} and "
                    f"target indices {t_idxs}. Only indices {sorted(intersection)} "
                    "are considered active."
                )
                t_idxs_list = [tuple(intersection)]
                d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)

    if draw:
        img = __get_highlight_molimg(mol, d_flag)
        return d_flag, img
    else:
        return d_flags_by_indices(t_idxs_list, d_all_atoms)
        else:
            # flatten all matched indices into a set
            matched_atoms = set()
            for match in matches:
                matched_atoms.update(match)

            intersection = t_idxs_set & matched_atoms

            if not intersection:
                info(
                    f"No overlap between SMARTS matches {matches} and "
                    f"target indices {t_idxs}. Returning atoms for the target indices only."
                )
                t_idxs_list = [t_idxs] if isinstance(t_idxs, tuple) else list(t_idxs)
                d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)
            else:
                info(
                    f"Partial overlap found between SMARTS matches {matches} and "
                    f"target indices {t_idxs}. Only indices {sorted(intersection)} "
                    "are considered active."
                )
                t_idxs_list = [tuple(intersection)]
                d_flag = __get_atoms_by_indices(t_idxs_list, d_all_atoms)

    if draw:
        img = __get_highlight_molimg(mol, d_flag)
        return d_flag, img
    else:
        return d_flag

