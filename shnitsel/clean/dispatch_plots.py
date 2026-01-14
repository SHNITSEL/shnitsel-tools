from typing import Literal, Sequence
import xarray as xr



# This function is in a separate file to avoid a circular dependency
# between `shnitsel.clean` and `shnitsel.vis.plot.filtration`


def dispatch_plots(
    filtranda: xr.DataArray,
    plot_thresholds: bool | Sequence[float],
    plot_populations: Literal['independent', 'intersections', False],
):
    """Call filtration-related plotting functions depending on arguments

    Parameters
    ----------
    filtranda
        Data according to which to filter
    plot_thresholds
        - If ``True``, will plot using ``check_thresholds`` with
        default quantiles
        - If a ``Sequence``, will plot using ``check_thresholds``
        with specified quantiles
        - If ``False``, will not plot threshold plot
    plot_populations
        - If ``'intersections'``, will plot populations of
        trajectories satisfying intersecting conditions
        - If ``'independent'``, will plot populations of
        trajectories satisfying conditions taken independently
        - If ``False``, will not plot populations plot


    Raises
    ------
    ValueError
        If ``plot_populations`` is an invalid value
    """
    from shnitsel.vis.plot.filtration import check_thresholds, validity_populations
    if plot_thresholds:
        if isinstance(plot_thresholds, Sequence):
            quantiles = plot_thresholds
        else:
            quantiles = None
        check_thresholds(filtranda, quantiles)
    if plot_populations == 'intersections':
        validity_populations(filtranda, intersections=True)
    elif plot_populations == 'independent':
        validity_populations(filtranda, intersections=False)
    elif plot_populations is not False:
        raise ValueError(f"Invalid value: {plot_populations=}")
