# from .filter_energy import (
#     energy_filtranda as energy_filtranda,
#     sanity_check as sanity_check,
# )
# from .filter_geo import (
#     bond_length_filtranda as bond_length_filtranda,
#     filter_by_length as filter_by_length,
# )
# from .common import (
#     omit as omit,
#     truncate as truncate,
#     transect as transect,
#     cum_max_quantiles as cum_max_quantiles,
#     true_upto as true_upto,
#     cum_mask_from_dataset as cum_mask_from_dataset,
#     cum_mask_from_filtranda as cum_mask_from_filtranda,
# )

from numbers import Number
from typing import Sequence
from typing_extensions import Literal

from shnitsel.clean.filter_energy import EnergyFiltrationThresholds, filter_by_energy
from shnitsel.clean.filter_geo import GeometryFiltrationThresholds, filter_by_length
from rdkit.Chem import Mol


def sanity_check(
    frames,
    filter_method: Literal["truncate", "omit", "annotate"] | Number = "truncate",
    *,
    energy_thresholds: EnergyFiltrationThresholds | None = None,
    geometry_thresholds: GeometryFiltrationThresholds | None = None,
    plot_thresholds: bool | Sequence[float] = False,
    plot_populations: bool | Literal["independent", "intersections"] = False,
    mol: Mol | None = None,
):
    """Filter trajectories according to energy to exclude unphysical (insane) behaviour

    Parameters
    ----------
    frames
        A xr.Dataset with an ``atXYZ`` variable as well as ``astate``, ``energy``, and ideally ``e_kin`` variables
    filter_method, optional
        Specifies the manner in which to remove data;

            - if 'omit', drop trajectories unless all frames meet criteria (:py:func:`shnitsel.clean.omit`)
            - if 'truncate', cut each trajectory off just before the first frame that doesn't meet criteria
                (:py:func:`shnitsel.clean.truncate`)
            - if 'annotate', merely annotate the data;
            - if a number, interpret this number as a time, and cut all trajectories off at this time,
                discarding those which violate criteria before reaching the given limit,
                (:py:func:`shnitsel.clean.transect`)
        see :py:func:`shnitsel.clean.dispatch_filter`.
    energy_thresholds, optional
        Threshold for total, potential and kinetic energy of the system.
        Can specify thresholds for overall drift and individual time step changes.
        Can also specify thresholds for energy steps at hops.
        Unit should be specified as a member variable.
        If not provided will default to some reasonable default values as seen in `EnergyThresholds` definition.
    geometry_thresholds, optional
        A mapping from SMARTS-strings to length-thresholds.

            - The SMARTS-strings describe bonds which are searched
                for in an RDKit Mol object obtained via :py:func:`shnitsel.bridges.default_mol`
            - The thresholds describe maximal tolerable bond-lengths; if there are multiple matches
                for a given search, the longest bond-length will be considered for each frame
            - The unit for the maximum length is provided in the member variable `length_unit` which defaults to `angstrom`.
            - If not provided will be initialized with thresholds for H-(C/N) bonds and one for all bonds.
    plot_thresholds
        See :py:func:`shnitsel.vis.plot.filtration.check_thresholds`.

        - If ``True``, will plot using ``check_thresholds`` with
        default quantiles
        - If a ``Sequence``, will plot using ``check_thresholds``
        with specified quantiles
        - If ``False``, will not plot threshold plot
    plot_populations
        See :py:func:`shnitsel.vis.plot.filtration.validity_populations`.

        - If ``True`` or ``'intersections'``, will plot populations of
        trajectories satisfying intersecting conditions
        - If ``'independent'``, will plot populations of
        trajectories satisfying conditions taken independently
        - If ``False``, will not plot populations plot
    Returns
    -------
        The sanitized xr.Dataset

    Notes
    -----
    The resulting object has a ``energy_filtranda`` and a ``geometry_filtranda`` data_var, representing the values by which the data were filtered.
    If the input has a ``filtranda`` data_var, it is overwritten.
    If the input has a `criterion` dimension, it will be dropped.
    """

    # Perform energy filtering
    ds_energy = filter_by_energy(
        frames,
        filter_method,
        energy_thresholds=energy_thresholds,
        plot_thresholds=plot_thresholds,
        plot_populations=plot_populations,
    )

    # Rename to filter-method prefixed names
    ds_tmp = ds_energy.rename_dims({"criterion": "energy_criterion"}).rename(
        {"filtranda": "energy_filtranda", "thresholds": "energy_thresholds"}
    )

    # Perform length filtering
    ds_lengths = filter_by_length(
        ds_tmp,
        filter_method,
        geometry_thresholds=geometry_thresholds,
        mol=mol,
        plot_thresholds=plot_thresholds,
        plot_populations=plot_populations,
    )

    # Rename to filter-method prefixed names
    ds_tmp = ds_lengths.rename_dims({"criterion": "length_criterion"}).rename(
        {"filtranda": "length_filtranda", "thresholds": "lengths_thresholds"}
    )

    return ds_tmp
