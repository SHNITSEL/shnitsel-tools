import rdkit
import rdkit.Chem as rc

import matplotlib as mpl

from ... import pca_biplot, postprocess as P
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
    ax.set_xlabel(f"SMILES={smiles}\n{inchi}", wrap=True)
    print(smiles, inchi)
    # axy.tick_params(axis="y", labelleft=False)
    return ax
