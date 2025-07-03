from logging import warning

import rdkit
import rdkit.Chem.Draw
import rdkit.Chem as rc

import matplotlib as mpl

from ..plot import pca_biplot
from .. import postprocess as P
from .common import figax

def xyz_to_mol(atXYZ, charge=0, covFactor=1.5) -> rc.Mol:
    mol = rc.rdmolfiles.MolFromXYZBlock(P.to_xyz(atXYZ))
    rc.rdDetermineBonds.DetermineBonds(
        mol, charge=charge, useVdw=True, covFactor=covFactor
    )
    rc.rdDepictor.Compute2DCoords(mol)
    return mol


def mol_to_png(mol, width=320, height=240):
    d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(width, height)

    d.drawOptions().setBackgroundColour((1, 1, 1, 0))
    d.drawOptions().padding = 0.05

    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()

# TODO DEPRECATE
def show_atXYZ(
    atXYZ, charge=0, name='', smiles=None, inchi=None, skeletal=True, ax=None
) -> mpl.axes.Axes:
    fig, ax = pca_biplot.figax(ax)

    mol = pca_biplot.xyz_to_mol(atXYZ, charge=charge)
    smol = rdkit.Chem.RemoveHs(mol)
    rdkit.Chem.RemoveStereochemistry(smol)
    smiles = rdkit.Chem.MolToSmiles(smol) if smiles is None else smiles
    inchi = rdkit.Chem.MolToInchi(smol) if inchi is None else inchi

    png = mol_to_png(rdkit.Chem.RemoveHs(mol) if skeletal else mol)
    pca_biplot.mpl_imshow_png(ax, png)
    ax.set_title(name)
    ax.axis('on')
    ax.get_yaxis().set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.set_xlabel(f"SMILES={smiles}\n{inchi}", wrap=True)
    print(smiles, inchi)
    # axy.tick_params(axis="y", labelleft=False)
    return ax

def format_inchi(inchi: str) -> str:
    if len(inchi) < 30:
        return inchi
    else:
        split = inchi.split('/')
        if len(split) not in {4, 5}:
            warning(f"Unexpected InChi: {split=}")
        lens = [len(s) for s in split]
        split[2] = '\n' + split[2]
        if sum(lens[2:]) > 30:
            split[3] = '\n' + split[3]
        return '/'.join(split)


def plot_structure(
    mol, name='', smiles=None, inchi=None, fig=None, ax=None
) -> mpl.axes.Axes:
    fig, ax = figax(fig, ax)
    png = mol_to_png(mol)
    pca_biplot.mpl_imshow_png(ax, png)
    ax.set_title(name)
    ax.axis('on')
    ax.get_yaxis().set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    inchi = format_inchi(inchi)
    ax.set_xlabel(f"SMILES={smiles}\n{inchi}", fontsize='small')
    print(smiles, inchi)
    # axy.tick_params(axis="y", labelleft=False)
    return ax
