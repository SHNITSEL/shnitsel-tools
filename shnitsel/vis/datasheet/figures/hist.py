from matplotlib.axes import Axes
import numpy as np
import numpy.typing as npt

from shnitsel.core._api_info import internal


# it will be useful to truncate some of the histograms
# this should be noted in the text
# a logarithmic histogram could show outliers... probably not worth it
@internal()
def calc_truncation_maximum(data, rel_cutoff: float = 0.01, bins: int = 1000) -> float:
    """Function to calculate the upper cutoff-threshold of data such that
    the frequency in the last bin is at least `rel_cutoff` times the maximum frequency.

    Helps to limit outliers.

    Parameters
    ----------
    data
        The data that should be histogrammed and filtered.
    rel_cutoff : float, optional
        Factor relative to the frequency maximum that should be used for determining the cutoff. Defaults to 0.01.
    bins : int, optional
        How many bins should be used for getting the correct threshold location. Defaults to 1000.

    Returns
    -------
    float
        Upper cutoff position to stay above the `rel_cutoff` relative threshold.
    """
    freqs, edges = np.histogram(
        data, bins=bins, range=(np.nanmin(data), np.nanmax(data))
    )
    cutoff = freqs.max() * rel_cutoff
    relevant = edges[:-1][freqs > cutoff]
    return relevant.max()


@internal()
def truncate_from_above(
    data: npt.NDArray, rel_cutoff: float = 0.01, bins: int = 1000
) -> npt.NDArray:
    """Helper function to truncate the `data` array on the upper end by a threshold
    such that the histogram frequency stays above `rel_cutoff*max(frequency)` and is below this relative cutoff above the cutoff.

    Parameters
    ----------
    data
        The data that should be histogrammed and filtered.
    rel_cutoff : float, optional
        Factor relative to the frequency maximum that should be used for determining the cutoff. Defaults to 0.01.
    bins : int, optional
        How many bins should be used for getting the correct threshold location. Defaults to 1000.

    Returns
    -------
    npt.NDArray
        The filtered data array
    """
    sup = calc_truncation_maximum(data, rel_cutoff=rel_cutoff, bins=bins)
    plot_data = data[data <= sup]
    return plot_data


@internal()
def create_marginals(ax: Axes) -> tuple[Axes, Axes]:
    """Function to create a pair of axes on top of and beside the axes passed as an argument to plot additional
    data into.

    Generally used for plotting dimension-specific histograms next to xy-plots.

    Parameters
    ----------
    ax : Axes
        Axes to create marginal plots inside of

    Returns
    -------
    tuple[Axes, Axes
        Resulting pair of outset axes.
    """
    axx = ax.inset_axes((0.0, 1.05, 1.0, 0.25), sharex=ax)
    axy = ax.inset_axes((1.05, 0.0, 0.25, 1.0), sharey=ax)
    # no labels next to main plot, no axis at all on other side
    axx.tick_params(axis="x", labelbottom=False)
    axy.tick_params(axis="y", labelleft=False)
    axx.get_yaxis().set_visible(False)
    axy.get_xaxis().set_visible(False)
    return axx, axy


def create_marginals_dict(axs: dict[str, Axes], label: str) -> dict[str, Axes]:
    """Function to add a set of marginal axes with `create_marginals()` for the axes at key `label` in `axs`
    and add the marginal axes back into the dict with appended `x` and `y` suffixes.

    Parameters
    ----------
    axs : dict[str, Axes]
        Dict of axes from which to pick the axis object and amend with marginal axes.
    label : str
        Key in `axs` for which the marginal axes should be created.

    Returns
    -------
    dict[str, Axes]
        `axs` but with the new marginal axes inserted at `{label}x` and `{label}y.
    """
    ax = axs[label]
    axx, axy = create_marginals(ax)
    axs[f'{label}x'], axs[f'{label}y'] = axx, axy
    return axs
