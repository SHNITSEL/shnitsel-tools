from functools import wraps

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from ...pca_biplot import figax as figax


# TODO: DEPRECATE
def figaxs(
    fig=Figure | None,
    axs=list[Axes] | dict[str, Axes] | None,
    mosaic=list | str | None,
    scale_factors=tuple[float, float],
):
    if scale_factors is None:
        scale_factors = (1, 1)
    set_size = fig is None and axs is None
    if fig is None:
        if len(plt.get_fignums()):
            fig = plt.gcf()
        else:
            fig = plt.figure(layout='constrained')
    if axs is None:
        axs = fig.subplot_mosaic(mosaic=mosaic)
    if set_size:
        fig.set_size_inches(8.27 * scale_factors[0], 11.69 * scale_factors[1])
    return fig, axs


def figaxs_defaults(mosaic, scale_factors=None, height_ratios=None):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, fig=None, axs=None, **kws):
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