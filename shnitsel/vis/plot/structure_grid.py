from math import ceil
from typing import Iterable, Sequence

from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from shnitsel.filtering.structure_selection import FeatureDescriptor, StructureSelection
from rdkit.Chem import Mol

from shnitsel.rd import highlight_features
from shnitsel.vis.plot.common import mpl_imshow_png


def _get_axs_grid(
    num_axs: int, ax_labels: Sequence[str]
) -> tuple[Figure, dict[str, Axes]]:
    """Helper function to get a mosaic grid of axes for a number of desired axes
    that is as quadratic as possible.

    Returns the created figure together with a dict mapping the provided labels for the axes to the
    constructed ``Axes`` objects.

    Parameters
    ----------
    num_axs : int
        The number of axes to be created in the grid. Will be padded to an appropriate number.
    ax_labels : Sequence[str]
        The labels to associate the axes with. Must be unique and collision-free, i.e. no two labels must be the same.

    Returns
    -------
    Figure
        The created figure holding all the axes.
    dict[str, Axes]
        The mapping from axis labels to ``Axes`` objects.
    """
    ncols = ceil(num_axs**0.5)
    nblanks = ncols - (num_axs % ncols)
    flat = ax_labels[:num_axs] + [None] * nblanks
    mosaic = np.array(flat).reshape(-1, ncols)
    fig, axs = plt.subplot_mosaic(mosaic.tolist())
    return fig, axs


def feature_highlight_grid(
    mol: Mol | StructureSelection,
    features: Sequence[Sequence[FeatureDescriptor | None]],
    labels: Sequence[str] | None = None,
    colors: dict[str, Sequence[tuple]] | None = None,
    cmaps: str | Colormap | dict[str, str | Colormap] | None = None,
    axs: dict[str, Axes] | None = None,
    width: int | None = None,
    height: int | None = None,
) -> tuple[Figure | None, dict[str, Axes]]:
    """Plot a grid highlighting the respective features in different axes
    using calls to RDKit functions to highlight the features on a molecular
    representation.

    Uses the `labels` object to set the title of each highlight plot
    and uses the `mol` object to draw the features onto.

    Parameters
    ----------
    mol : Mol
        An RDKit ``Mol`` object to be used for structure display
    features : Sequence[Iterable[FeatureDescriptor | None]]
        The tuples describing the various features on a molecule as defined in the
        `shnitsel.filtering.structure_selection` module.
    labels : Sequence[str], optional
        Labels for the plots of the feature graphs
    colors : dict[str, Sequence[tuple]], optional
        Optional parameter to directly set the colors of each feature in the `features` sequence associated with the `label` key.
        If not set, colors may be picked from `cmaps` or from any default setting of `highlight_features()`.
    cmaps : str | Colormap |  dict[str, str | Colormap], optional
        Optional parameter to set a colormap for all features, or per axis.
        If not set, `highlight_features()` will apply default colors.
    axs : dict[str, Axes], optional
        A dictionary mapping from plot labels to :py:class:`matplotlib.pyplot.axes.Axes`
        objects
        (If not provided, one will be created.)
    width: int, optional
        Optional argument to specify the width of a molecular structure plot in pixels. Defaults to 320 on a lower level
    height: int, optional
        Optional argument to specify the height of a molecular structure plot. Defaults to 240 on a lower level
    """
    assert mol is not None, (
        "No feature plotting can be performed without a `Mol` object describing the structure."
    )
    plot_mol: Mol
    if isinstance(mol, StructureSelection):
        assert mol.mol is not None, (
            "The provided structure selection did not contain a `Mol` object to use for feature highlights."
        )
        plot_mol = mol.mol
    else:
        plot_mol = mol

    if labels is None:
        labels = list('abcdefghijklmnopqrstuvwxyz')

    # Get axes grid if not provided
    if axs is None:
        num_plots = min(len(features), len(labels))
        fig, axs = _get_axs_grid(num_plots, labels)
    else:
        fig = None

    if cmaps is not None:
        if isinstance(cmaps, (str, Colormap)):
            # turn into a dict
            cmaps = {label: cmaps for label in labels}

    # Disable axes labels and other decoration
    for mol_ax in axs.values():
        mol_ax.axis('off')

    if axs is not None and plot_mol is not None:
        # Plot the highlighted features
        for label, feature_set in zip(labels, features):
            if label in axs:
                # Acquire optional settings for colors from parameters
                f_colors = None
                f_cmap = None
                if colors is not None and label in colors:
                    f_colors = colors[label]
                elif cmaps is not None and label in cmaps:
                    f_cmap = cmaps[label]

                png = highlight_features(
                    plot_mol,
                    feature_set,
                    colors=f_colors,
                    cmap=f_cmap,
                    fmt='png',
                    width=width,
                    height=height,
                )
                mpl_imshow_png(axs[label], png)
                axs[label].set_title(label)

    return fig, axs
