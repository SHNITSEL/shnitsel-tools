from collections.abc import Sequence
from functools import wraps
from typing import Hashable

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.typing import HashableList
from matplotlib.text import Text

from ...plot.common import figax as figax

symbols = dict(energy=r"$E_i$", force=r"$\mathbf{F}_i$", dip_perm=r"$\mathbf{\mu}_i$")


def figaxs_defaults(
    mosaic: list[HashableList[Hashable]],
    scale_factors: Sequence[float] | None = None,
    height_ratios: Sequence[float] | None = None,
):
    """Decorator to automatically create a mosaic of subfigures and provide the axes to the decorated function if only a figure is provided.

    Args:
        mosaic (list[HashableList[Hashable]]): Matrix of keys, where the individual subplots should go
        scale_factors (Sequence[float] , optional): Sequence of scale factors for the individual plots. Defaults to None.
        height_ratios (Sequence[float] , optional): Height ratios of the individual plots. Defaults to None.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(
            *args, fig: Figure | None = None, axs: dict[str, Axes] | None = None, **kws
        ):
            nonlocal func, scale_factors, mosaic, height_ratios
            if scale_factors is None:
                scale_factors = (1, 1)
            set_size = fig is None and axs is None
            if fig is None:
                if len(plt.get_fignums()):
                    fig = plt.gcf()
                else:
                    fig = plt.figure(layout='constrained')
            if axs is None:
                axs = fig.subplot_mosaic(mosaic=mosaic, height_ratios=height_ratios)
            if set_size:
                fig.set_size_inches(8.27 * scale_factors[0], 11.69 * scale_factors[1])
            return func(*args, fig=fig, axs=axs, **kws)

        return wrapper

    return decorator


def centertext(
    text: str,
    ax: Axes,
    clearticks='y',
    background_color: ColorType | None = None,
    color: ColorType | None = None,
) -> Text:
    """Helper method to center the text within the axes.

    Optionally removes ticks in the dimensions `x` or `y`.

    Args:
        text (str): Message to center in the frame
        ax (Axes): Axes to plot the text into
        clearticks (str, optional): String of all dimensions to clear the ticks for (may contain `x` and/or `y`). Defaults to 'y'.
        background_color (ColorType|None, optional): Color argument to set for the background of the plot
        color (ColorType|None, optional): Color argument to set for font on the plot

    Returns:
        Text: The Text object created by a call to `.text()` on the `ax` object.
    """
    if 'x' in clearticks:
        ax.tick_params(axis='x', labelbottom=False)
    if 'y' in clearticks:
        ax.tick_params(axis='y', labelleft=False)
    if background_color is not None:
        ax.set_facecolor(background_color)
    return ax.text(
        0.5, 0.5, text, transform=ax.transAxes, ha='center', va='center', color=color
    )
