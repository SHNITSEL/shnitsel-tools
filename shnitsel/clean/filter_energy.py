from dataclasses import asdict, dataclass
import logging
from typing import Literal, Sequence, TypeVar

import numpy as np
from shnitsel.data.dataset_containers import Frames, Trajectory, wrap_dataset
import xarray as xr

from shnitsel.data.dataset_containers.data_series import DataSeries
from shnitsel.data.multi_indices import mdiff
from shnitsel.clean.common import dispatch_filter
from shnitsel.clean.dispatch_plots import dispatch_plots
from shnitsel.units.conversion import convert_energy
from shnitsel.units.definitions import energy

TrajectoryOrFrames = TypeVar(
    "TrajectoryOrFrames", bound=Trajectory | Frames | xr.Dataset
)


@dataclass
class EnergyFiltrationThresholds:
    """Helper class to keep an extensible set of threshold values for various physical
    filtration properties. The keys/names of the fields are the properties that can be checked
    and the float values are the energy thresholds that should be applied.
    Additionally, this class has a `energy_unit` field which specifies the unit that all other
    fields are given in (by default: `eV`)
    """

    epot_active_step: float = 0.7
    # TODO:FIXME: The hop step should be smaller than the epot active step.
    epot_hop_step: float = 1.0
    etot_step: float = 0.1
    etot_drift: float = 0.2
    ekin_step: float = 0.7
    energy_unit: str = energy.eV

    def __init__(self, kv_map: dict[str, float] | None = None):
        if kv_map is not None:
            kv_map = dict(kv_map)
            if "epot_active_step" in kv_map:
                self.epot_active_step = kv_map["epot_active_step"]
                del kv_map["epot_active_step"]
            if "epot_hop_step" in kv_map:
                self.epot_hop_step = kv_map["epot_hop_step"]
                del kv_map["epot_hop_step"]
            if "etot_step" in kv_map:
                self.etot_step = kv_map["etot_step"]
                del kv_map["etot_step"]
            if "etot_drift" in kv_map:
                self.etot_drift = kv_map["etot_drift"]
                del kv_map["etot_drift"]
            if "ekin_step" in kv_map:
                self.ekin_step = kv_map["ekin_step"]
                del kv_map["ekin_step"]

            if kv_map:
                logging.warning(
                    f"Unsupported energy filtration criteria: {kv_map.keys}"
                )

    def to_dataarray(
        self, selected_criteria: Sequence[str] | None = None
    ) -> xr.DataArray:
        """Helper function to convert this dataclass object into an xarray.DataArray to be
        assigned to the coordinate of a Filtration dataset.

        Args:
            selected_criteria (Sequence[str] | None, optional):
                The sequence of criteria keys to be a result of the conversion. Defaults to None.
                Must be a subset of the keys of fields on this object if not set to None.
                If set to None, all but the unit field in this object will be used to set the
                list of criterion keys.


        Returns:
            xr.DataArray:
                A DataArray with a 'criterion' coordinate with either all available or all
                selected properties and the threshold values in the array cells.

        """
        dict_repr = asdict(self)
        if selected_criteria is None:
            criteria = list(dict_repr.keys())
            # Don't include the unit property
            criteria.remove("energy_unit")
        else:
            # Make sure all criteria exist on this object
            assert all(x in dict_repr.keys() for x in selected_criteria)
            criteria = list(selected_criteria)

        threshold_values = [dict_repr[c] for c in criteria]
        # Build DataArray from thresholds, criterion names and energy unit
        res = xr.DataArray(
            list(threshold_values),
            coords={"criterion": criteria},
            attrs={"units": self.energy_unit},
        )
        return res.astype(float)


def calculate_energy_filtranda(
    frames_or_trajectory: xr.Dataset | DataSeries,
    *,
    energy_thresholds: dict[str, float] | EnergyFiltrationThresholds | None = None,
) -> xr.DataArray:
    """Derive energetic filtration targets from an xr.Dataset

    Parameters
    ----------
    frames
        A xr.Dataset with ``astate``, ``energy``, and ideally ``e_kin`` variables
    energy_thresholds, optional
        Threshold for total, potential and kinetic energy of the system.
        Can specify thresholds for overall drift and individual time step changes.
        Can also specify thresholds for energy steps at hops.
        Unit should be specified as a member variable.
        If not provided will default to some reasonable default values as seen in `EnergyFiltrationThresholds` definition.

    Returns
    -------
        An xr.DataArray of filtration targets stacked along the ``criterion`` dimension;
        criteria comprise epot_step and hop_epot_step, as well as
        etot_drift, etot_step and ekin_step if the input contains an e_kin variable
    """
    if isinstance(frames_or_trajectory, xr.Dataset):
        frames_or_trajectory = wrap_dataset(frames_or_trajectory, DataSeries)

    elif not isinstance(frames_or_trajectory, DataSeries):
        message: str = 'Filtered dataset object is of type %s instead of the required types Frames or Trajectory'
        logging.warning(message, type(frames_or_trajectory))
        raise ValueError(message % type(frames_or_trajectory))

    if energy_thresholds is None or not isinstance(
        energy_thresholds, EnergyFiltrationThresholds
    ):
        energy_thresholds = EnergyFiltrationThresholds(energy_thresholds)

    if not frames_or_trajectory.has_data('astate'):
        message: str = 'Skipping active energy filtering because of missing variable `astate` in the trajectory'
        logging.warning(message)
        # print(frames_or_trajectory)
        raise ValueError(message)

    if energy_thresholds is None:
        energy_thresholds = EnergyFiltrationThresholds()

    filter_energy_unit = energy_thresholds.energy_unit

    is_hop = mdiff(frames_or_trajectory.active_state) != 0

    res = xr.Dataset()
    if frames_or_trajectory.has_data('energy'):
        # TODO: FIXME: Shouldn't we drop coords instead?
        e_pot_active = frames_or_trajectory.energy.sel(
            state=frames_or_trajectory.active_state
        ).drop_vars("state")
        e_pot_active.attrs["units"] = frames_or_trajectory.energy.attrs["units"]
        e_pot_active = convert_energy(e_pot_active, to=filter_energy_unit)

        res["epot_active_step"] = mdiff(e_pot_active).where(~is_hop, 0)
        res["epot_hop_step"] = mdiff(e_pot_active).where(is_hop, 0)
    else:
        logging.warning(
            'Skipping active state energy filtering because of missing variable `energy` in the trajectory'
        )
        e_pot_active = None

    if frames_or_trajectory.has_data("e_kin"):
        e_kin = frames_or_trajectory.e_kin
        e_kin.attrs["units"] = frames_or_trajectory.e_kin.attrs["units"]
        e_kin = convert_energy(e_kin, to=filter_energy_unit)
        # TODO: FIXME: Do we really only care about the ekin difference at hopping points?
        res["ekin_step"] = mdiff(e_kin).where(~is_hop, 0)

        if e_pot_active is not None:
            e_tot = e_pot_active + e_kin
            res["etot_step"] = mdiff(e_tot)
            # FIXME (thevro): Use more general way to determine correct groupby spec, if any
            if 'trajid' in e_tot.coords:
                res["etot_drift"] = e_tot.groupby('trajid').map(lambda x: x - x.item(0))
            elif 'atrajectory' in e_tot.coords:
                res["etot_drift"] = e_tot.groupby('atrajectory').map(
                    lambda x: x - x.item(0)
                )
            else:
                res["etot_drift"] = e_tot - e_tot.item(0)
    else:
        e_kin = None
        logging.warning("data does not contain kinetic energy variable ('e_kin')")

    abs_criterion = np.abs(res.to_dataarray("criterion"))

    da = abs_criterion.assign_attrs(units=filter_energy_unit)

    # Make threshold coordinates
    da = da.assign_coords(
        thresholds=energy_thresholds.to_dataarray(
            selected_criteria=da.coords["criterion"].values
        )
    )
    return da


def filter_by_energy(
    frames_or_trajectory: TrajectoryOrFrames,
    filter_method: Literal["truncate", "omit", "annotate"] | float = "truncate",
    *,
    energy_thresholds: dict[str, float] | EnergyFiltrationThresholds | None = None,
    plot_thresholds: bool | Sequence[float] = False,
    plot_populations: Literal["independent", "intersections", False] = False,
) -> TrajectoryOrFrames | None:
    """Filter trajectories according to energy to exclude unphysical (insane) behaviour

    Parameters
    ----------
    frames_or_trajectory
        A Frames or Trajectory object with ``astate``, ``energy``, and ideally ``e_kin`` variables.
        If ``astate`` is not set, no filtering will be performed and no filtranda assigned.
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
    plot_thresholds
        See :py:func:`shnitsel.vis.plot.filtration.check_thresholds`.

        - If ``True``, will plot using ``check_thresholds`` with
        default quantiles
        - If a ``Sequence``, will plot using ``check_thresholds``
        with specified quantiles
        - If ``False`` (the default), will not plot threshold plot
    plot_populations
        See :py:func:`shnitsel.vis.plot.filtration.validity_populations`.

        - If ``'intersections'``, will plot populations of
        trajectories satisfying intersecting conditions
        - If ``'independent'``, will plot populations of
        trajectories satisfying conditions taken independently
        - If ``False`` (the default), will not plot populations plot
    Returns
    -------
        The sanitized xr.Dataset

    Notes
    -----
    The resulting object has a ``filtranda`` data_var, representing the values by which the data were filtered.
    If the input has a ``filtranda`` data_var, it is overwritten.
    """

    analysis_data: Trajectory | Frames = wrap_dataset(
        frames_or_trajectory, Trajectory | Frames
    )

    filtranda = calculate_energy_filtranda(
        analysis_data, energy_thresholds=energy_thresholds
    )
    dispatch_plots(filtranda, plot_thresholds, plot_populations)
    # Here we need to build a new object with the criteria assigned.

    filtered_frames = type(analysis_data)(
        analysis_data.dataset.drop_dims(["criterion"], errors="ignore").assign(
            filtranda=filtranda
        )
    )
    filter_res = dispatch_filter(filtered_frames, filter_method)

    if isinstance(frames_or_trajectory, xr.Dataset):
        return filter_res.dataset
    else:
        return filter_res
