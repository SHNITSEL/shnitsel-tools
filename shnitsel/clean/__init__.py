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

import logging
from typing import Any, Sequence, TypeVar
from typing_extensions import Literal

from xarray import Dataset

from shnitsel.bridges import construct_default_mol
from shnitsel.core._api_info import API, internal
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.structure_selection import SMARTSstring


from .filter_energy import EnergyFiltrationThresholds, filter_by_energy
from .filter_geo import GeometryFiltrationThresholds, filter_by_length
from rdkit.Chem import Mol

TrajectoryOrFrames = TypeVar("TrajectoryOrFrames", bound=Trajectory | Frames)


@API()
def sanity_check(
    trajectory_or_frames: TreeNode[Any, TrajectoryOrFrames] | TrajectoryOrFrames,
    filter_method: Literal["truncate", "omit", "annotate"] | float = "truncate",
    *,
    energy_thresholds: dict[str, float] | EnergyFiltrationThresholds | None = None,
    geometry_thresholds: dict[SMARTSstring, float]
    | GeometryFiltrationThresholds
    | None = None,
    plot_thresholds: bool | Sequence[float] = False,
    plot_populations: Literal["independent", "intersections", False] = False,
    mol: Mol | None = None,
    drop_empty_trajectories: bool = False,
) -> TreeNode[Any, TrajectoryOrFrames] | TrajectoryOrFrames | None:
    """Filter trajectories according to energy to exclude unphysical (insane) behaviour

    Parameters
    ----------
    trajectory_or_frames : Trajectory | Frames | TreeNode[Any, Trajectory|Frames]
        A Trajectory or Frames object (or a ShnitselDB structure holding such objects) with an ``atXYZ`` variable as well as ``astate``, ``energy``, and ideally ``e_kin`` variables
    filter_method : Literal["truncate", "omit", "annotate"] | float, optional
        Specifies the manner in which to remove data;
            - if 'omit', drop trajectories unless all frames meet criteria (:py:func:`shnitsel.clean.omit`)
            - if 'truncate', cut each trajectory off just before the first frame that doesn't meet criteria
                (:py:func:`shnitsel.clean.truncate`)
            - if 'annotate', merely annotate the data;
            - if a `float` number, interpret this number as a time, and cut all trajectories off at this time,
                discarding those which violate criteria before reaching the given limit,
                (:py:func:`shnitsel.clean.transect`)
        see :py:func:`shnitsel.clean.dispatch_filter`.
    energy_thresholds : EnergyFiltrationThresholds, optional
        Threshold for total, potential and kinetic energy of the system.
        Can specify thresholds for overall drift and individual time step changes.
        Can also specify thresholds for energy steps at hops.
        Unit should be specified as a member variable.
        If not provided will default to some reasonable default values as seen in `EnergyThresholds` definition.
    geometry_thresholds : GeometryFiltrationThresholds, optional
        A mapping from SMARTS-strings to length-thresholds.

            - The SMARTS-strings describe bonds which are searched
                for in an RDKit Mol object obtained via :py:func:`shnitsel.bridges.default_mol`
            - The thresholds describe maximal tolerable bond-lengths; if there are multiple matches
                for a given search, the longest bond-length will be considered for each frame
            - The unit for the maximum length is provided in the member variable `length_unit` which defaults to `angstrom`.
            - If not provided will be initialized with thresholds for H-(C/N) bonds and one for all bonds.
    plot_thresholds : bool, optional
        See :py:func:`shnitsel.vis.plot.filtration.check_thresholds`.

        - If ``True``, will plot using ``check_thresholds`` with
        default quantiles
        - If a ``Sequence``, will plot using ``check_thresholds``
        with specified quantiles
        - If ``False``, will not plot threshold plot
    plot_populations : Literal ['intersections', 'independent', False], optional
        See :py:func:`shnitsel.vis.plot.filtration.validity_populations`.

        - If ``'intersections'``, will plot populations of
        trajectories satisfying intersecting conditions
        - If ``'independent'``, will plot populations of
        trajectories satisfying conditions taken independently
        - If ``False``, will not plot populations plot
    mol : rdkit.Chem.Mol, optional
        Optional parameter to provide a mol object to base structure analysis on, by default generated from the first frame in the trajectory or frameset.
    drop_empty_trajectories : bool, optional
        Flag to not include trajectories for which the sanity check result was empty in the final result tree, by default False.
        Only used for tree-structure inputs.

    Returns
    -------
        The sanitized trajectory, frames or tree.
        A tree is sanitized by applying the sanitization function to all individual data points in the tree.

    Notes
    -----
    The resulting object has a ``energy_filtranda`` and a ``geometry_filtranda`` data_var, representing the values by which the data were filtered.
    If the input has a ``filtranda`` data_var, it is overwritten.
    If the input has a `criterion` dimension, it will be dropped.
    """
    kws = dict(
        filter_method=filter_method,
        energy_thresholds=energy_thresholds,
        geometry_thresholds=geometry_thresholds,
        plot_populations=plot_populations,
        plot_thresholds=plot_thresholds,
        mol=mol,
    )
    if isinstance(trajectory_or_frames, TreeNode):
        return trajectory_or_frames.map_data(
            lambda x: _sanity_check_per_trajectory(x, **kws),
            keep_empty_branches=not drop_empty_trajectories,
        )
    else:
        return _sanity_check_per_trajectory(trajectory_or_frames, **kws)


@internal()
def _sanity_check_per_trajectory(
    trajectory_or_frames: TrajectoryOrFrames,
    filter_method: Literal["truncate", "omit", "annotate"] | float = "truncate",
    *,
    energy_thresholds: dict[str, float] | EnergyFiltrationThresholds | None = None,
    geometry_thresholds: dict[SMARTSstring, float]
    | GeometryFiltrationThresholds
    | None = None,
    plot_thresholds: bool | Sequence[float] = False,
    plot_populations: Literal["independent", "intersections", False] = False,
    mol: Mol | None = None,
) -> TrajectoryOrFrames | None:
    """Internal function to filter a single Trajectory or Frames object

    Parameters
    ----------
    trajectory_or_frames
        A xr.Dataset with an ``atXYZ`` variable as well as ``astate``, ``energy``, and ideally ``e_kin`` variables
    filter_method, optional
        Specifies the manner in which to remove data;

            - if 'omit', drop trajectories unless all frames meet criteria (:py:func:`shnitsel.clean.omit`)
            - if 'truncate', cut each trajectory off just before the first frame that doesn't meet criteria
                (:py:func:`shnitsel.clean.truncate`)
            - if 'annotate', merely annotate the data;
            - if a `float` number, interpret this number as a time, and cut all trajectories off at this time,
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

        - If ``'intersections'``, will plot populations of
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

    wrapped_ds = wrap_dataset(trajectory_or_frames, Trajectory | Frames)

    assert isinstance(wrapped_ds, (Trajectory, Frames)), (
        "Data provided to `sanity_check()` is neither Trajectory nor Frame data."
    )

    if mol is None:
        mol = construct_default_mol(wrapped_ds)

    # Perform energy filtering
    ds_energy = filter_by_energy(
        wrapped_ds,
        filter_method,
        energy_thresholds=energy_thresholds,
        plot_thresholds=plot_thresholds,
        plot_populations=plot_populations,
    )

    if ds_energy is None:
        logging.info("Rejected trajectory because of energy constraints")
        return None

    rename_keys = [
        "filtranda",
        "thresholds",
        "filter_mask",
        "good_upto",
        "good_throughout",
        "criterion",
    ]
    prefix = "energy"
    # Rename to filter-method prefixed names
    ds_tmp = type(ds_energy)(
        ds_energy.dataset.rename_dims({"criterion": prefix + "_criterion"}).rename_vars(
            {key: prefix + "_" + key for key in rename_keys if key in ds_energy.dataset}
        )
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

    if ds_lengths is None:
        logging.info("Rejected trajectory because of length constraints")
        return None
    prefix = "length"

    # Rename to filter-method prefixed names
    ds_tmp = type(ds_lengths)(
        ds_lengths.dataset.rename_dims(
            {"criterion": prefix + "_criterion"}
        ).rename_vars(
            {
                key: prefix + "_" + key
                for key in rename_keys
                if key in ds_lengths.dataset
            }
        )
    )

    # If input was unwrapped, then unwrap result
    if isinstance(trajectory_or_frames, Dataset):
        return ds_tmp.dataset

    return ds_tmp


__all__ = ['sanity_check', 'filter_by_energy', 'filter_by_length']
