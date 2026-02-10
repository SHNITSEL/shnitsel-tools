"""This module contains functions that accept an RDKit.Chem.Mol object;
but *not* necessarily functions that *return* a Mol object."""

from typing import Iterable, Literal, Mapping, Sequence, TYPE_CHECKING

from matplotlib.colors import Colormap
import rdkit.Chem as rc
import rdkit.Chem.rdDetermineBonds  # noqa: F401
import matplotlib as mpl
import numpy as np

from shnitsel.vis.colormaps import hex2rgb

if TYPE_CHECKING:
    from shnitsel.filtering.structure_selection import FeatureDescriptor

#################################################
# Functions for converting RDKit objects to
# SMILES annotated with the original atom indices
# to maintain the order in the `atom` index


def set_atom_props(
    mol: rc.Mol,
    inplace: bool = False,
    **kws: Sequence[str | int]
    | Literal[True]
    | Literal[False]
    | Mapping[int, int | str]
    | None,
) -> rc.Mol | None:
    """Set properties on atoms of an ``rdkit.Chem.Mol`` object

    Parameters
    ----------
    mol : rdkit.Chem.mol
        The ``Mol`` object
    inplace : bool, optional
        Whether to alter ``mol``; , by default False (returns a copy)
    **kws : Sequence[str | int] | Literal[True] | Mapping[int, int | str] | None
        A mapping where parameter names represent the name of a property
        and the arguments are either
            - a dict mapping the atom indices to values that should be assigned. Missing atom indices are ignored
            - a sequence of str or int values the atoms should be set to;
            - ``True``, in which case the atom indices as
              assigned by RDKit will be used as values;
            - ``False``, in which the property will be cleared on every atom.
            - ``None`` values, which are simply ignored

    Returns
    -------
        A copy of the ``Mol`` object, if ``inplace=False``, otherwise the provided mol

    Raises
    ------
    ValueError
        If the amount of values passed to the properties kwargs did not match the amount of
        atoms in the mol.
    """
    if not inplace:
        mol = rc.Mol(mol)
    natoms = mol.GetNumAtoms()
    for prop, vals in kws.items():
        if vals is None:
            # None is not assigned, just ignored
            continue
        elif isinstance(vals, dict):
            # Try and assign the values to each atom identified by the key
            for atom_id, val in vals.items():
                atom = mol.GetAtomWithIdx(atom_id)
                if atom is not None:
                    atom.SetProp(prop, str(val))
        else:
            if vals is True:
                # atom indices are values if `vals=True`
                vals = range(natoms)
            elif vals is False:
                # `False` means, the value should be cleared.
                for atom in mol.GetAtoms():
                    atom.ClearProp(prop)
                continue
            elif natoms != len(vals):
                raise ValueError(
                    f"{len(vals)} values were passed for {prop}, but 'mol' has {natoms} atoms"
                )

            for atom, val in zip(mol.GetAtoms(), vals):
                atom.SetProp(prop, str(val))
    return mol


def mol_to_numbered_smiles(mol: rc.Mol) -> str:
    """Generate a SMILES string containing mapping numbers
    corresponding to the atom indices in the mol object

    Parameters
    ----------
    mol
        An ``rdkit.Chem.Mol`` object

    Returns
    -------
        A SMILES string

    Notes
    -----
        This is intended as a way to store the connectivity
        and order of a matrix of coordinates
    """
    mol = rc.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return rc.MolToSmiles(mol)


def highlight_features(
    mol: rc.Mol,
    feature_indices: Sequence["FeatureDescriptor | None"],
    fmt: Literal['png', 'svg'] = 'png',
    colors: Sequence[tuple] | None = None,
    cmap: str | Colormap | None = None,
    width: int | None = None,
    height: int | None = None,
):
    """Highlight specified pairs of atoms in an image of an ``rdkit.Chem.Mol`` object

    Parameters
    ----------
    mol : rc.Mol
        The ``Mol`` object
    feature_indices : Sequence[FeatureDescriptor | None]
        A list of tuples of indices for various features.
    fmt : {'png', 'svg'}, default='png'
        Optional format specifier for the desired image output format.
    colors : Sequence[tuple], optional
        Optionally a sequence of colors to plot the features in. If set, does not consider the value in `cmap`.
        If not provided, will be generated from `cmap`.
    cmap : str | Colormap, optional
        Either the name of a `Colormap` or a `Colormap` to use to assign colors to the features if `colors` is not set.
        If required but not provided defaults to the `rainbow` colormap.
    width: int, default=320
        The width of the target graphic in pixels
    height: int, default=320
        The height of the target graphic in pixels

    Returns
    -------
        Raw PNG data or the SVG string depending on the choice of `fmt`.
    """
    if width is None:
        width = 320
    if height is None:
        height = 320

    if fmt == 'png':
        d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(width, height)
    elif fmt == 'svg':
        d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    else:
        raise ValueError(
            f"Unsupported format for molecular highlight graph: {fmt=}. Supported are 'png' and 'svg'"
        )

    plot_colors: Iterable[tuple]
    if colors is None:
        if cmap is None:
            # colors = iter(mpl.colormaps['tab10'](range(10)))
            plot_colors = iter(
                mpl.colormaps['rainbow'](np.linspace(0, 1, len(feature_indices)))
            )
        else:
            if isinstance(cmap, str):
                cmap = mpl.colormaps[cmap]
            plot_colors = cmap(np.linspace(0, 1, len(feature_indices)))
    else:
        assert len(colors) >= len(feature_indices), (
            "Insufficient number of colors for the number of features provided."
        )
        plot_colors = colors
    atom_colors: dict[int, list[tuple[float, float, float]]] = {}
    bond_colors: dict[int, list[tuple[float, float, float]]] = {}

    # We expect here that there are sufficient colors for our features
    for feature, c in zip(feature_indices, plot_colors):
        if feature is None:
            continue

        atoms_to_color = []

        if isinstance(c, str):
            c = tuple(hex2rgb(c))

        if not isinstance(c, tuple):
            c = tuple(c)

        if isinstance(feature, int):
            # Position
            a = feature
            atoms_to_color = [a]
        elif isinstance(feature, tuple):
            match feature:
                case a1, (a2, a3, a4):
                    # pyramid
                    if a1 < 0 or a2 < 0 or a3 < 0 or a4 < 0:
                        continue
                    # Mark bonds
                    for other in (a2, a3, a4):
                        if (bond := mol.GetBondBetweenAtoms(a1, other)) is not None:
                            bondid = bond.GetIdx()
                            if bondid not in bond_colors:
                                bond_colors[bondid] = []
                            bond_colors[bond.GetIdx()].append(c)
                    # Mark all atoms
                    atoms_to_color = [a1, a2, a3, a4]
                case (a1, a2):
                    # Bond
                    a1, a2 = feature
                    if a1 < 0 or a2 < 0:
                        continue

                    if (bond := mol.GetBondBetweenAtoms(a1, a2)) is not None:
                        bondid = bond.GetIdx()
                        if bondid not in bond_colors:
                            bond_colors[bondid] = []
                        bond_colors[bond.GetIdx()].append(c)
                    else:
                        atoms_to_color = [a1, a2]
                case _:
                    # angle or dihedral
                    # Mark bonds
                    if not all(x >= 0 for x in feature):
                        continue

                    flen = len(feature)

                    for i in range(flen - 1):
                        a1, a2 = feature[i], feature[i + 1]
                        if isinstance(a1, int) and isinstance(a2, int):
                            if (bond := mol.GetBondBetweenAtoms(a1, a2)) is not None:
                                bondid = bond.GetIdx()
                                if bondid not in bond_colors:
                                    bond_colors[bondid] = []
                                bond_colors[bond.GetIdx()].append(c)

                    # Mark all atoms
                    atoms_to_color = [x for x in feature if isinstance(x, int)]

        # Mark all atoms
        for a in atoms_to_color:
            if a < 0:
                # Ignore negative entries
                continue

            if a not in atom_colors:
                atom_colors[a] = []
            atom_colors[a].append(c)

    # d.drawOptions().fillHighlights = False
    d.drawOptions().setBackgroundColour((0.8, 0.8, 0.8, 0.5))
    d.drawOptions().padding = 0
    # TODO: FIXME: We seem to have double indices because of this
    d.drawOptions().addAtomIndices = True

    d.DrawMoleculeWithHighlights(mol, '', atom_colors, bond_colors, {}, {}, -1)
    d.FinishDrawing()
    return d.GetDrawingText()
