import logging
from logging import warning, info

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

from shnitsel.filtering.structure_selection import StructureSelection
from .flagging import __collect_tuples
from .helpers import __level_to_geo
from shnitsel.vis.colormaps import hex2rgb, st_yellow


def __get_color_atoms(d_flag, flag_level):
    flagged_tuples = d_flag[__level_to_geo(flag_level)]
    color_atoms = [a for flag, at, *_ in flagged_tuples if flag == 1 for a in at]

    return color_atoms


def __get_color_bonds(d_flag, flag_level):
    flagged_tuples = d_flag[__level_to_geo(flag_level)]
    color_bonds = [b for flag, at, bt, *_ in flagged_tuples if flag == 1 for b in bt]

    return color_bonds


def __get_highlight_molimg(
    mol: Mol,
    d_flag,
    highlight_color: tuple[float, float, float] | str = st_yellow,
    width=300,
    height=300,
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

    if isinstance(highlight_color, str):
        highlight_color = tuple(*hex2rgb(highlight_color))

    flagged_tuples = __collect_tuples(d_flag)

    # draw molecule with highlights
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().fillHighlights = True
    drawer.drawOptions().addAtomIndices = True
    drawer.drawOptions().setHighlightColour(highlight_color)
    drawer.drawOptions().clearBackground = False

    has_binfo = all(len(t) >= 5 for t in flagged_tuples)

    if has_binfo:
        highlights = [(at, bt) for flag, at, bt, bot, m in flagged_tuples if flag == 1]
        highlight_bonds = [bidx for at, bt in highlights for bidx in bt]
        highlight_atoms = [aidx for at, bt in highlights for aidx in at]

        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer, mol, highlightAtoms=highlight_atoms, highlightBonds=highlight_bonds
        )

    else:
        highlight_atoms = [idx for flag, (idx, symbol) in flagged_tuples if flag == 1]

        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=highlight_atoms)

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img_text = drawer.GetDrawingText()
    img = SVG(img_text)

    return img


def __get_img_multiple_mols(
    mol: Mol, d_multi_flag: dict, l_patterns: list, l_levels: list
) -> SVG:
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
    opts.setHighlightColour(tuple(*hex2rgb(st_yellow)))

    img = Draw.MolsToGridImage(
        mols,
        highlightAtomLists=highlight_atoms,
        highlightBondLists=highlight_bonds,
        molsPerRow=len(mols),
        subImgSize=(300, 300),
        useSVG=True,
        drawOptions=opts,
    )

    return img
