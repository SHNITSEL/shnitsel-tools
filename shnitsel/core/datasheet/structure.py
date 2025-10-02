from logging import warning

import rdkit
import rdkit.Chem.Draw

import matplotlib as mpl

from ..plot import pca_biplot
from .common import figax

def mol_to_png(mol, width=320, height=240):
    d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(width, height)

    d.drawOptions().setBackgroundColour((1, 1, 1, 0))
    d.drawOptions().padding = 0.05

    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()


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
