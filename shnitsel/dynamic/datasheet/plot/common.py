import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


def figaxs(
    fig=Figure | None,
    axs=Axes | dict | None,
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


# Could make a decorator later... could cause more problems than it solves
def figaxs_wrapper(func, *args, fig=None, axs=None, **kws):
    return NotImplemented