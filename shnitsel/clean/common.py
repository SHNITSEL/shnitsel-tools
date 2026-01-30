import logging
from typing import Literal, Sequence, TypeVar

import numpy as np
import xarray as xr

from shnitsel.core.typedefs import DimName
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory

from shnitsel.data.multi_indices import ensure_unstacked, stack_trajs

########################
# Formats we will need #
########################


def true_upto(mask: xr.DataArray, dim: str) -> xr.DataArray:
    """Helper function to assess whether a mask has only `true` entries up until a certain point
    along dimension `dim`.
    Used to check if criterion validity is maintained along the `time` dimension.

    Returns array with values of `dim` coordinate up to which the values are all `true` or `-np.inf` if no
    frame is valid.

    Parameters
    ----------
    mask : xr.DataArray
        The mask holding boolean flags whether criteria are fulfilled.
    dim : str
        The dimension along which to check continuous validity of criteria.

    Returns
    -------
    xr.DataArray
        The point in time up to which the criterion is fulfilled.
    """
    mask, was_stacked = ensure_unstacked(mask.astype(bool), fill_value=False)
    assert dim in mask.dims, "Mask array is missing specified dimension %s" % dim
    shifted_coord = np.concat([[-np.inf], mask.coords[dim].data])
    num_cum_valid_indices = mask.cumprod(dim).sum(dim).astype(int)
    res = np.take(shifted_coord, num_cum_valid_indices.data)
    # We only deal with individual trajectories
    return num_cum_valid_indices.copy(data=res)


_true_upto = true_upto


def _filter_mask_from_criterion_mask(mask: xr.DataArray) -> xr.DataArray:
    """Generate cutoff array from the mask, specifying for each criterion, up to which point
    the criterion is fulfilled.

    Either holds a boolean `filter_mask` or a time values variable `good_upto` depending on
    whether a `time` dimension/coordinate is present.

    Parameters
    ----------
    mask : xr.DataArray
        The xarray holding the boolean flags whether a frame contains valid data for various criteria

    Returns
    -------
    xr.DataArray
        With name `filter_mask` which holds true boolean flags, whether a frame should be kept according to the respective criterion.
        If a time dimension is present, also holds a `good_upto` coordinate, which maps
        criteria to the time value at which the criterion is last fulfilled.
        If the time dimension is missing, will just have boolean flags.
        Also has a coordinate `good_throughout`, which indicates, whether the entire trajectory/frameset satisfies the criterion.
    """
    leading_dim: str
    if "time" in mask.dims:
        leading_dim = "time"
        good_upto = _true_upto(mask, leading_dim)
        good_upto.name = "good_upto"
        filter_mask = (mask.time < good_upto).astype(bool)
        filter_mask.name = "filter_mask"
        filter_mask = filter_mask.assign_coords(good_upto=good_upto)
    else:
        # We have independent frames. Don't try and conceive a `good up to this point` property
        leading_dim = "frame"
        filter_mask = mask.copy()
        filter_mask.name = "filter_mask"

    good_throughout = mask.all(leading_dim)

    return filter_mask.assign_coords(good_throughout=good_throughout)


def _filter_mask_from_filtranda(filtranda: xr.DataArray) -> xr.DataArray:
    """
    Calculates first a filter mask and then the cutoffs from that mask using
    `_cutoffs_from_mask`
    """
    thresholds = filtranda.coords["thresholds"]
    is_good_frame = (filtranda < thresholds).astype(bool)
    return _filter_mask_from_criterion_mask(is_good_frame)


def _filter_mask_from_dataset(ds: xr.Dataset) -> xr.DataArray:
    """
    Returns a da containing cutoff times (the same as the good_upto data_var)
    and with a coord called good_throughout

    The returned object has dimension {'criterion'}.

    Parameters
    ----------
    ds : xr.Dataset
        A Dataset containing either:
            - a 'good_upto' data_var and a 'good_throughout' coordinate
            - a 'filtranda' data_var with a 'threshold' coordinate

    Returns
    -------
        A DataArray containing cutoff times (the same as the good_upto data_var)
        and with a coord called good_throughout
        The returned object has dimension {'criterion'}.

    Raises
    ------
    ValueError
        If there is no filtration information in the Dataset
    """
    if "good_upto" in ds.data_vars:
        mask = (ds.coords["time"] <= ds["good_upto"]).astype(bool)
        mask.name = "filter_mask"
        return mask
    elif "filter_mask" in ds.data_vars:
        mask = ds.data_vars["filter_mask"]
        if "good_throughout" not in mask.coords:
            raise ValueError(
                "data_var 'good_upto' is missing expected coord "
                "'good_throughout'; will recalculate."
            )
        else:
            return mask
    elif "filtranda" in ds.data_vars:
        return _filter_mask_from_filtranda(ds.data_vars["filtranda"])
    else:
        raise ValueError(
            "Please set data_vars 'filtranda' and 'thresholds', "
            "or alternatively supply a filter mask directly using data_var 'filter_mask'"
        )


####################
# Action functions #
####################

# All the action functions take a dataset
# They can use the functions above to get the info they need

TrajectoryOrFrames = TypeVar("TrajectoryOrFrames", bound=Trajectory | Frames)


def omit(frames_or_trajectory: TrajectoryOrFrames) -> TrajectoryOrFrames | None:
    """If all filter criteria are fulfilled throughout, keep the trajectory.
    Otherwise return None to omit it.

    Parameters
    ----------
    frames_or_trajectory : Frames | Trajectory
        Either the Frameset or the trajectory to filter

    Returns
    -------
    Frames | Trajectory | None
        The Frameset or Trajectory if all filter conditions are fulfilled or None if it should be omitted.
    """
    wrapped_dataset = wrap_dataset(frames_or_trajectory)
    try:
        filter_mask = _filter_mask_from_dataset(frames_or_trajectory.dataset)
        good_throughout = filter_mask["good_throughout"]
        all_critera_fulfilled = good_throughout.all("criterion").item()
        if all_critera_fulfilled:
            return frames_or_trajectory
    except:
        pass
    # FIXME (thevro): This doesn't work for MultiFrames aka. stacked
    return None


_omit = omit


def _log_omit(before, after):
    kept = set(after.trajid.values.tolist())
    omitted = set(before.trajid.values.tolist()).difference(kept)
    logging.info(
        f"Kept {len(kept)} trajectories, IDs: {kept}; \n"
        f"Dropped {len(omitted)} trajectories, IDs: {omitted}"
    )


def truncate(
    frames_or_trajectory: TrajectoryOrFrames | xr.Dataset,
) -> TrajectoryOrFrames | Trajectory | Frames:
    """Perform a truncation on the trajectory or frameset, i.e. cut off the trajectory
    after the last frame that fulfils all filtration conditions.

    Parameters
    ----------
    frames_or_trajectory : TrajectoryOrFrames | xr.Dataset
        The dataset to truncate

    Returns
    -------
    TrajectoryOrFrames | Trajectory | Frames
        The truncated dataset.
    """

    wrapped_dataset = wrap_dataset(frames_or_trajectory, Trajectory | Frames)

    filter_mask_all_criteria = _filter_mask_from_dataset(wrapped_dataset.dataset).all(
        "criterion"
    )

    tmp_res = wrapped_dataset.dataset.isel(
        {frames_or_trajectory.leading_dim: filter_mask_all_criteria}
    )
    # TODO: FIXME: Test whether this works. May be wrong shape
    if not isinstance(frames_or_trajectory, xr.Dataset):
        return type(frames_or_trajectory)(tmp_res)
    else:
        return wrap_dataset(tmp_res, Trajectory | Frames)


_truncate = truncate


def transect(
    trajectory: Trajectory | xr.Dataset, cutoff_time: float
) -> Trajectory | None:
    """Perform a transect, i.e. cut off the trajetory at time `cutoff_time` if it is valid until then
    or omit it, if it is not valid for long enough.

    Trajectory must be a trajectory with `time` dimension.

    Parameters
    ----------
    trajectory : Trajectory | xr.Dataset
        The trajectory to transect
    cutoff_time : float
        Time at which the trajectory should be cut off or discarded entirely if conditions are not satisfied until this time.

    Returns
    -------
    Trajectory | None
        Either the filtered trajectory with all frames being valid up until `cutoff_time` or None if the trajectory is not valid for long enough.
    """
    wrapped_dataset = wrap_dataset(trajectory, Trajectory)

    assert "time" in wrapped_dataset.dims, (
        "Dataset has no coordinate `time` but time-based truncation has been requested, which cannot be performed!"
    )

    time_sliced_dataset = wrapped_dataset.loc[{"time": slice(float(cutoff_time))}]
    good_upto = _filter_mask_from_dataset(time_sliced_dataset)
    assert good_upto.name == "good_upto", (
        "Despite a `time` dimension being present, the filter mask returned for the dataset was not a `good_upto` value."
    )
    # TODO: FIXME: We may want to accept the last time before `cutoff_time` to be true.
    is_trajectory_good = (good_upto >= cutoff_time).all("criterion").item()
    if is_trajectory_good:
        return Trajectory(time_sliced_dataset)
    else:
        return None


_transect = transect


def dispatch_filter(
    frames_or_trajectory: TrajectoryOrFrames,
    filter_method: Literal["truncate", "omit", "annotate"] | float = "truncate",
) -> TrajectoryOrFrames | None:
    """Filter trajectories according to energy to exclude unphysical (insane) behaviour

    Parameters
    ----------
    frames_or_trajectory
        A Frames or Trajectory object with a `filtranda` variable set and a `thresholds` coordinate both along a `criterion` dimension.
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

    Returns
    ----------
        The modified dataset with either data violating the

    Raises
    ----------
    ValueError
        If an unsupported value for the `cut` parameter was provided.
    """
    if not frames_or_trajectory.has_variable(
        "filtranda"
    ) or not frames_or_trajectory.has_coordinate("thresholds"):
        logging.warning(
            "Trajectory is missing the required variable `filtranda` (possibly because they could not be calculated) or the required coordinate `thresholds` for those filtranda. No filtering is performed."
        )
        return frames_or_trajectory

    if filter_method == "annotate":
        filter_mask = _filter_mask_from_dataset(frames_or_trajectory.dataset)
        return type(frames_or_trajectory)(
            frames_or_trajectory.dataset.assign({filter_mask.name: filter_mask})
        )
    elif filter_method == "truncate":
        return truncate(frames_or_trajectory)
    elif filter_method == "omit":
        return omit(frames_or_trajectory)
    elif isinstance(filter_method, float):
        transect_position: float = filter_method
        # assert isinstance(frames_or_trajectory, Trajectory), (
        #     "Cannot provide a `Frames` object to a `transect()` call. Unsupported operation."
        # )
        return transect(frames_or_trajectory, transect_position)  # type: ignore # Must be Trajectory here or the _transect() call will error out anyway.

    else:
        raise ValueError(
            "`filter_method` should be one of {'truncate', 'omit', 'annotate'}, or a number, "
            f"not {filter_method}"
        )


###########################################
# Formats directly prerequisite for plots #
###########################################


def cum_max_quantiles(
    filtranda_array: xr.DataArray,
    quantiles: Sequence[float] | None = None,
    cum_dim: DimName = "time",
    group_dim: DimName = "trajectory",
) -> xr.DataArray:
    """Quantiles of cumulative maxima

    Parameters
    ----------
    filtranda_array : xr.DataArray
        A DataArray, or a Dataset with a data_var 'filtranda';
        either way, the Variable should have dimensions and
        coordinates corresponding to a
        (stacked or unstacked) ensemble of trajectories.
    quantiles : Sequence[float], optional
        Which quantiles to calculate,
        by default ``[0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]``
    cum_dim : DimName, optional
        The dimension along which to accumulate the maxima, by default `time`
    group_dim : DimName, optional
        The key/dimension along which to calculate the quantiles of the maxima, by default `atrajectory`.

    Returns
    -------
    xr.DataArray
        A DataArray with 'quantile' and 'cum_dim' dimensions;
        'group_dim' dimension will have been removed to calculate quantiles;
        other dimensions remain unaffected.

    See also
    --------
    The data returned by this function is intended for consumption by
    :py:func:`shnitsel.vis.plot.filtration.check_thresholds`
    """
    # # NOTE (thevro): This function doesn't accept stacked, only unstacked (layered)
    # filtranda_array, _ = ensure_unstacked(filtranda_array)

    if quantiles is None:
        quantiles = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]

    filtranda_array = filtranda_array.fillna(0)
    time_axis = filtranda_array.get_axis_num(cum_dim)

    cum_max = filtranda_array.copy(
        data=np.maximum.accumulate(filtranda_array.data, axis=time_axis)
    )
    # TODO: FIXME: Rewrite this to allow for accumulation across collection of trajectories, not only stacked.
    return cum_max.quantile(quantiles, group_dim)
