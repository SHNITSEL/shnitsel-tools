import sys, os
import logging
from logging import warning, info
import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
import geomatch

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    force=True
)

def __flag_alternating_double_single(d_flag: dict) -> dict:
    """
    Identify double–single–double patterns in bonds from a flagged bond dictionary.

    Only considers bonds with flag==1. Flags the bonds as:
        2 -> double bonds in chain (outer bonds and inner)
        1 -> single bond strictly between two double bonds
        0 -> all other bonds

    If there are fewer than three bonds flagged with 1, the function
    cannot find alternating double–single–double patterns and returns
    all flags set to 0 with an info message.

    Parameters
    ----------
    d_flag : dict
        Dictionary with key 'bonds' containing tuples:
        (flag, (atom1, atom2), (bond_idx,), (bond_order,), mol)

    Returns
    -------
    dict
        Updated bonds dictionary with first element in each tuple set according to rules.
    """
    atoms = d_flag['atoms']
    bonds = d_flag['bonds']

    # Only consider bonds with flag == 1
    active_indices = [i for i, (flag, *_) in enumerate(bonds) if flag == 1]

    if len(active_indices) < 3:
        info(f"Only {len(active_indices)} found, i.e., less than three bonds are flagged. "\
                f"Cannot find alternating double–single–double patterns.")
        return {
                'atoms': atoms,
                'bonds': [(0,) + t[1:] for t in bonds]
                }

    # Build atom -> bonds mapping
    atom_to_bonds = {}
    for i in active_indices:
        _, (a1, a2), _, _, _ = bonds[i]
        atom_to_bonds.setdefault(a1, []).append(i)
        atom_to_bonds.setdefault(a2, []).append(i)

    # Initialize result flags
    result_flags = [0] * len(bonds)

    # Find all active double bonds
    double_bonds = [i for i in active_indices if bonds[i][3][0] == 2.0]

    # Process each double bond as potential start
    for start_idx in double_bonds:
        if result_flags[start_idx] != 0:
            continue  # already processed

        result_flags[start_idx] = 2  # flag first double bond
        chain_queue = [(start_idx, 2.0)]  # bond_idx, last bond order

        while chain_queue:
            bond_idx, last_order = chain_queue.pop(0)
            _, (a1, a2), _, (order,), _ = bonds[bond_idx]
            neighbor_atoms = [a1, a2]
            next_order = 1.0 if last_order == 2.0 else 2.0

            for atom in neighbor_atoms:
                for nb_idx in atom_to_bonds.get(atom, []):
                    if result_flags[nb_idx] != 0:
                        continue
                    nb_order = bonds[nb_idx][3][0]
                    if nb_order == next_order:
                        chain_queue.append((nb_idx, nb_order))
                        # flag single bonds only if between two doubles
                        if nb_order == 1.0:
                            a, b = bonds[nb_idx][1]
                            neighbors = []
                            for neighbor_atom in (a, b):
                                neighbors.extend([
                                    i for i in atom_to_bonds.get(neighbor_atom, [])
                                    if bonds[i][3][0] == 2.0
                                ])
                            if len(set(neighbors)) >= 2:
                                result_flags[nb_idx] = 1
                        elif nb_order == 2.0:
                            result_flags[nb_idx] = 2

    # Build updated bonds dictionary
    updated_bonds = [(result_flags[i],) + bonds[i][1:] for i in range(len(bonds))]

    return {
            'atoms': atoms,
            'bonds': updated_bonds
            }

def flag_bla(
        mol: Mol,
        smarts=None,
        t_idxs=[]) -> dict:

    d_flag_bonds = geomatch.flag_bonds(mol, smarts=smarts, t_idxs=t_idxs, draw=False)
    d_flag_bla = __flag_alternating_double_single(d_flag_bonds)

    return d_flag_bla



