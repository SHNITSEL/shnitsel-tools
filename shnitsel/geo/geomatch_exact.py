from logging import info
from typing import Literal

import rdkit.Chem as rc
import xarray as xr

from shnitsel._contracts import needs
from shnitsel.bridges import default_mol
from shnitsel.geo.geocalc import get_bats
from shnitsel.geo import geomatch


def flag_exact(mol: rc.Mol, l_smarts: list[str] | str) -> dict:
    """
    Compute and flag bonds, angles, and dihedrals in a single call,
    for multiple structural features

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule under study.
    l_smarts : list[str], optional
        SMARTS patterns for filtering interactions.

    Returns
    -------
    dict
        {
            'bonds':      [...],
            'angles':     [...],
            'dihedrals':  [...]
        }
    Note that even if no BATs of a certain type (e.g. no bonds) are
    requested, the returned dictionary will still contain the corresponding
    key, associated with an empty list.
    """
    if isinstance(l_smarts, str):
        l_smarts = [l_smarts]
    res = {'bonds': [], 'angles': [], 'dihedrals': []}
    for smarts in l_smarts:
        patt = rc.MolFromSmarts(smarts)
        if patt is None:
            raise ValueError(f"Invalid SMARTS '{smarts}'")
        if not mol.HasSubstructMatch(patt):
            info(f"No matches for pattern '{smarts}', skipping.")
            continue

        n_atoms = patt.GetNumAtoms()

        d_flag = geomatch.flag_bats(mol, smarts)
        if n_atoms == 2:
            res['bonds'].extend(d_flag[0]['bonds'])
        elif n_atoms == 3:
            res['angles'].extend(d_flag[0]['angles'])
        elif n_atoms == 4:
            res['dihedrals'].extend(d_flag[0]['dihedrals'])

    return res


@needs(dims={'atom', 'direction'})
def get_bats_matching(
    atXYZ: xr.DataArray,
    l_smarts: list[str],
    signed: bool | None = None,
    ang: Literal[False, 'deg', 'rad'] = False,
) -> xr.DataArray:
    """Get bond lengths, angles and torsions according to a list of SMARTS searches.

    Parameters
    ----------
    atXYZ
        The coordinates to use.
    l_smarts : list[str], optional
        SMARTS patterns to search for.
    signed, optional
        Whether to distinguish between clockwise and anticlockwise rotation,
        when returning angles as opposed to cosine & sine values;
        by default, do not distinguish.
        NB. This applies only to the dihedrals, not to the three-center angles.
        The latter are always unsigned.
    ang, optional
        If False (the default), returns sines and cosines;
        if set to 'deg', returns angles in degrees
        if set to 'rad', returns angles in radians

    Returns
    -------
        An :py:class:`xarray.DataArray` containing bond lengths, angles and tensions as specified.

    Examples
    --------
        >>> geom.get_bats_matching(frames['atXYZ'], ['C~C=C~C', '[#7]~[#6]'])
        # Finds and calculates all-carbon torsions with a central double bond
        # and bond-lengths between carbon and nitrogen.
    """
    mol = default_mol(atXYZ)
    matches = flag_exact(mol, l_smarts)
    return get_bats(atXYZ, matches_or_mol=matches, signed=signed, ang=ang, pyr=False)