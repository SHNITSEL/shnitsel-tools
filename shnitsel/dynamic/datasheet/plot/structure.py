import rdkit
from ... import pca_biplot


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
):
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

def plot_structure(mol, name='', smiles=None, inchi=None, ax=None):
    fig, ax = pca_biplot.figax(ax)
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
