import logging
from logging import warning, info

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

from shnitsel.filtering.structure_selection import StructureSelection
from .geomatch.features.atoms import flag_atoms
from .geomatch.features.bonds import __get_all_bonds, __get_bonds_by_indices
from .geomatch.flagging import __get_bond_info, __match_pattern
from .geomatch.presentation import __get_highlight_molimg
from shnitsel.vis.colormaps import st_yellow

# -----------------------------------------------------------------------------
# --- BLA ---------------------------------------------------------------------
# -----------------------------------------------------------------------------


def __build_conjugated_smarts(n_double: int, elems: str = "#6,#7,#8,#15,#16") -> str:
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


def __match_bla_chromophor(
    mol,
    smarts: str | None = None,
    n_double: int | None = None,
    elems: str = "#6,#7,#8,#15,#16",
):
    """
    Detect conjugated chromophores defined either by SMARTS, number of
    alternating double bonds, or automatically by maximum extension.

    Decision logic
    --------------
    1. If `smarts` is given -> use SMARTS
    2. If both `smarts` and `n_double` are given -> validate consistency
    3. If only `n_double` is given -> generate SMARTS
    4. If neither is given -> search for maximum conjugated chromophore

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
            1 for b in patt.GetBonds() if b.GetBondType() == Chem.rdchem.BondType.DOUBLE
        )

        if dbl >= n_double:
            take = "smarts"
        else:
            take = "n_double"

        info(
            f"SMARTS encodes {dbl} double bonds, "
            f"n_double={n_double}. Using {take} for pattern matching."
        )

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

    return {"smarts": smarts, "n_double": n_double, "matches": matches}


def flag_bla_chromophor(
    mol,
    smarts: str | None = None,
    n_double: int | None = None,
    elems: str = "#6,#7,#8,#15,#16",
    draw=True,
    width: int = 500,
    height: int = 300,
):
    match_chromo = __match_bla_chromophor(mol, smarts, n_double, elems)

    if not match_chromo['matches']:
        info(
            f"SMARTS ({match_chromo['smarts']}) and/or n_double "
            f"({match_chromo['n_double']}) not found in the molecule. "
            f"No bonds will be flagged."
        )
        d_flag = {}
        return d_flag

    else:
        out_atoms = flag_atoms(
            mol=mol, smarts=match_chromo['smarts'], t_idxs=None, draw=False
        )['atoms']

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
