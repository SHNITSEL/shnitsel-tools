import logging
from typing import Any, Literal, TypeVar, overload
import re

import xarray as xr
import numpy as np

from shnitsel.core.typedefs import DimName
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.data_series import DataSeries
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.multi_indices import mdiff, sel_trajs
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.state_selection import StateSelection, StateSelectionDescriptor
from shnitsel.filtering.helpers import _get_default_state_selection


TrajectoryOrFrames = TypeVar("TrajectoryOrFrames", bound=DataSeries)


# TODO: Finish documentation
@overload
def hops_mask_from_active_state(
    active_state_source: xr.Dataset | xr.DataArray | DataSeries,
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
    dim: DimName | None = None,
) -> xr.DataArray:
    """Overload to specify simple return type for simple (flat) input types.

    See `hops_mask_from_active_state()` for details
    """
    ...


@overload
def hops_mask_from_active_state(
    active_state_source: TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
    dim: DimName | None = None,
) -> TreeNode[Any, xr.DataArray]:
    """Overload to specify hierarchical return type for hierarchical input types.

    See other overloads of `hops_mask_from_active_state` for details of how the individual data entries are
    mapped.
    """
    ...


def hops_mask_from_active_state(
    active_state_source: xr.Dataset
    | xr.DataArray
    | DataSeries
    | TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
    dim: DimName | None = None,
) -> xr.DataArray | TreeNode[Any, xr.DataArray]:
    """Generate boolean masks marking hopping points by identifying changes in the active state of provided
    data source.

    Needs to be fed either with (hierarchical) trajectory data that has `active_state` (`astate`) information
    or directly with the xr.DataArray holding `astate` information.

    Parameters
    ----------
    active_state_source : xr.Dataset | Trajectory | Frames | xr.DataArray | TreeNode[Any, xr.Dataset | Trajectory  |  Frames  |  xr.DataArray]
        A potential source for extracting the active state along a leading dimension and the leading dimension name.
    hop_type_selection: StateSelection | StateSelectionDescriptor, optional
        A state selection holding state transitions that should be used in hop filtering.
    dim : DimName, optional,
        The dimension along which the hops should be detected. For most cases, this should be `frame` or `time`.

    Returns
    -------
    xr.DataArray | TreeNode[Any, xr.DataArray]
        Either the flat boolean mask of leading dimension instances where a hop happens or a hierarchical structure holding
        such a flat mask for every original data entry in the hierarchical input data

    Raises
    ------
    ValueError
        If an unsupported input type was provided
    """
    if isinstance(active_state_source, TreeNode):
        return active_state_source.map_data(
            hops_mask_from_active_state,
            hop_type_selection=hop_type_selection,
            dtype=xr.DataArray,
        )
    else:
        active_state_data: xr.DataArray
        leading_dim: str | None
        if isinstance(active_state_source, xr.Dataset):
            tmp = wrap_dataset(active_state_source, DataSeries)
            active_state_data = tmp.active_state
            leading_dim = tmp.leading_dim
        elif isinstance(active_state_source, ShnitselDataset):
            active_state_data = active_state_source.active_state
            leading_dim = active_state_source.leading_dim
        elif isinstance(active_state_source, xr.DataArray):
            if 'astate' in active_state_source.coords:
                active_state_data = active_state_source.coords['astate']
            else:
                active_state_data = active_state_source
            leading_dim = str(active_state_data.dims[0])
        else:
            raise ValueError(
                "Unknown type of provided source for `active_state` data: %s"
                % type(active_state_source)
            )

        # Overwrite the leading dim detection
        if dim is not None:
            leading_dim = dim

        is_hop_mask = mdiff(active_state_data, dim=leading_dim) != 0
        # Add prior and current state back as coordinates
        is_hop_mask = is_hop_mask.assign_coords(
            hop_from=(active_state_data.shift({leading_dim: 1}, -1)),
            hop_to=active_state_data,
        )

        if hop_type_selection is not None and leading_dim is not None:
            state_selection = _get_default_state_selection(
                hop_type_selection, active_state_source
            )
            type_filter = np.full(is_hop_mask.sizes[leading_dim], False)
            # TODO: FIXME: We need to make sure that the state combinations returned are actually bidirectional is not directed.
            for hop_from, hop_to in state_selection.state_combinations:
                type_filter |= (is_hop_mask.hop_from == hop_from) & (
                    is_hop_mask.hop_to == hop_to
                )
            if not state_selection.is_directed:
                for hop_to, hop_from in state_selection.state_combinations:
                    type_filter |= (is_hop_mask.hop_from == hop_from) & (
                        is_hop_mask.hop_to == hop_to
                    )
            is_hop_mask &= type_filter
        return is_hop_mask


# TODO: FIXME: Make compatible with trees and wrapper datasets
def hops(
    frames, hop_type_selection: StateSelection | StateSelectionDescriptor | None = None
):
    """Select hops

    Parameters
    ----------
    frames
        An Xarray object (Dataset or DataArray) with a ``frames`` dimension
    hop_type_selection
        A list of pairs of states, e.g.:
        ``[(1, 2), (2, 1), (3, 1)]``
        to select only hops between states 1 and 2 as well as from
        3 to 1 (but not from 1 to 3).
        Alternatively, hops may be specified as a single string
        in the following style: ``'1<>2, 3->1'`` -- this specification
        selects the same hops as in the previous example, with ``<>``
        selecting hops in either direction and ``->`` being one-
        directional.

    Returns
    -------
    An indexed version of ``frames``, where each entry in the
    ``frames`` dimension represents a hop.
    The following coordinates are added along ``frames``:

        - ``tidx``: the time-step index of the hop in its trajectory
        - ``hop_from``: the active state before the hop
        - ``hop_to``: the active state after the hop
    """
    is_hop_mask = hops_mask_from_active_state(
        frames, hop_type_selection=hop_type_selection
    )

    res = frames.isel(frame=is_hop_mask)
    tidxs = np.concat(
        [np.arange(traj.sizes['frame']) for _, traj in frames.groupby('atrajectory')]
    )
    hop_tidx = tidxs[is_hop_mask]
    res = res.assign_coords(
        tidx=('frame', hop_tidx),
        hop_from=(frames['astate'].shift({'frame': 1}, -1).isel(frame=is_hop_mask)),
        hop_to=res['astate'],
    )
    if hasattr(res, 'drop_dims'):
        res = res.drop_dims(['trajectory'], errors='ignore')
    return res


@overload
def filter_data_at_hops(
    active_state_and_data_source: xr.DataArray,
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
) -> xr.DataArray:
    """Overload of `filter_data_at_hops` for non-hierarchical input types to indicate non-hierarchical returns.

    Converts raw datasets into a Shnitsel-Style wrapper
    """
    ...


@overload
def filter_data_at_hops(
    active_state_and_data_source: xr.Dataset | DataSeries,
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
) -> DataSeries:
    """Overload of `filter_data_at_hops` for non-hierarchical input types to indicate non-hierarchical returns.

    Converts raw datasets into a Shnitsel-Style wrapper
    """
    ...


@overload
def filter_data_at_hops(
    active_state_and_data_source: TreeNode[Any, xr.Dataset | DataSeries],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
) -> TreeNode[Any, DataSeries]:
    """Overload of `filter_data_at_hops` for hierarchical input types to indicate hierarchical returns.

    Converts raw datasets into a Shnitsel-Style wrapper
    """


@overload
def filter_data_at_hops(
    active_state_and_data_source: TreeNode[Any, xr.DataArray],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
) -> TreeNode[Any, xr.DataArray]:
    """Overload of `filter_data_at_hops` for hierarchical input types to indicate hierarchical returns.

    Converts raw datasets into a Shnitsel-Style wrapper
    """


def filter_data_at_hops(
    active_state_and_data_source: xr.Dataset
    | xr.DataArray
    | DataSeries
    | TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
) -> (
    DataSeries | xr.DataArray | TreeNode[Any, DataSeries] | TreeNode[Any, xr.DataArray]
):
    """Filter data to only retain data at points where hops of selected transitions occur.

    Needs to be fed either with (hierarchical) trajectory data that has `active_state` (`astate`) information
    or with simple (single) trajectory data with an `active_state` (`astate`) variable.

    If you wish to perform arbitrary filtering, you can employ the `hops_mask_from_active_state()` function
    to just get the boolean mask of hop positions and perform the filtering yourself.

    Parameters
    ----------
    active_state_and_data_source : xr.Dataset |  xr.DataArray | DataSeries | TreeNode[Any, xr.Dataset | DataSeries]| TreeNode[Any, xr.DataArray]
        A source for extracting the active state along a leading dimension and the leading dimension name as well as to filter
        the data from.
    hop_type_selection: StateSelection, optional
        A state selection holding state transitions that should be used in hop filtering.

    Returns
    -------
    DataSeries | xr.DataArray | TreeNode[Any, DataSeries] | TreeNode[Any, xr.DataArray]
        Filtered version of the input data source, where a selection was performed to only retain data at points where
        a hop was happening.
        Each entry along the leading dimension (`time` or `frame`) represents a hop.
        The following coordinates are added along the leading dimension:
            - ``hop_from``: the active state before the hop
            - ``hop_to``: the active state after the hop
    """

    if isinstance(active_state_and_data_source, TreeNode):
        return active_state_and_data_source.map_data(
            lambda x: filter_data_at_hops(x, hop_type_selection=hop_type_selection)
        )
    else:
        if isinstance(active_state_and_data_source, xr.DataArray):
            is_hop_mask = hops_mask_from_active_state(
                active_state_source=active_state_and_data_source,
                hop_type_selection=hop_type_selection,
            )

            hop_dim = list(is_hop_mask.sizes.keys())[0]
            return active_state_and_data_source[{hop_dim: is_hop_mask}]
        else:
            # Frames or Trajectory
            input_dataset = wrap_dataset(active_state_and_data_source, DataSeries)
            is_hop_mask = hops_mask_from_active_state(
                active_state_source=active_state_and_data_source,
                hop_type_selection=hop_type_selection,
            )
            # This introduces the coordinates for is_hop_mask, namely the mask of hopping point flags, the hop_from and hop_to coordinates.
            tmp_dataset = input_dataset.assign_coords(is_hop_mask=is_hop_mask)
            res_dataset = tmp_dataset.sel(is_hop_mask=True)

            # time_step_idxs = np.concat(
            #     [np.arange(traj.sizes['frame']) for _, traj in frames.groupby('trajid')]
            # )
            # hop_tidx = tidxs[is_hop]
            # res = res.assign_coords(
            #     # tidx=('frame', hop_tidx),
            #     hop_from=(frames.astate.shift({'frame': 1}, -1).isel(frame=is_hop)),
            #     hop_to=res.astate,
            # # )
            # if hop_type_selection is not None:
            #     acc = np.full(res.sizes['frame'], False)
            #     for hop_from, hop_to in hop_type_selection:
            #         acc |= (res.hop_from == hop_from) & (res.hop_to == hop_to)
            #     res = res.isel(frame=acc)
            # return res.drop_dims(['trajid_'], errors='ignore')

            # We drop the mask that should be all True values now.
            return wrap_dataset(
                res_dataset.drop("is_hop_mask", errors="ignore"), DataSeries
            )


# TODO: FIXME: Make StateSelection the preferred type for picking hopping types.
def focus_hops(
    frames: xr.Dataset
    | xr.DataArray
    | DataSeries
    | TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
    window: slice | None = None,
):
    """For each hop, create a copy of its trajectory centered on the hop; align these

    Parameters
    ----------
    frames : xr.Dataset | xr.DataArray | DataSeries | TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries]
        An Xarray object (Dataset or DataArray) with a ``frames`` dimension
    hop_type_selection : StateSelection | StateSelectionDescriptor, optional
        Types of hops to include
        See like-named parameter in :py:func:`shnitsel.analyze.hops.hops`
    window : slice, optional
        Clip the range of hop-relative times included to this slice;
        values are interpreted as relative times (not indices), e.g.:
        ``focus_hops(..., window=slice(-1.5, 2.5))``
        clips each trajectorial copy to the region between 1.5 time-units (probably fs)
        before the hop to 2.5 time-units after the hop.

    Returns
    -------
    An object with ``hop`` and ``hop_time`` dimensions.
    Each entry in ``hop`` represents a trajectory;
    there is one trajectory per hop, so possibly more than
    one copy of a given trajectory in the object.
    The following coordinates are added along dimension ``hop_time``, and do not vary
    by hop (i.e. do not contain a ``hop`` dimension):

        - ``hop_time``: the trajectory time coordinate relative to the hop
        - ``hop_tidx``: the trajectory time-step index relative to the hop

    The following coordinates are added along dimensions ``hop_time`` and ``hop``:

        - ``time``: the original time coordinate relative to the start of the trajectory
        - ``tidx``: the trajectory time-step index relative to the start of the trajectory

    The following coordinates are added along dimension ``hop``:

        - ``hop_from``: the active state before the hop
        - ``hop_to``: the active state after the hop
        - ``trajid``: the ID of the trajectory in which the hop occurred
    """
    raise NotImplementedError()
    # TODO: FIXME: Refactor this to new wrapper types
    hop_vals = hops(frames, hop_type_selection=hop_type_selection)
    # If no hops, return empty
    if hop_vals.sizes["frame"] == 0:
        res = frames.isel(frame=[])
        res = res.swap_dims({"frame": "hop_time"})
        res = res.drop_vars(["frame", "trajid", "time"])
        res = res.drop_dims(["trajid_"], errors="ignore")
        res = res.expand_dims("hop").isel(hop=[])
        empty_2d = xr.Variable(("hop", "hop_time"), [[]]).isel(hop=[], hop_time=[])
        res = res.assign_coords(
            {
                "hop_time": ("hop_time", []),
                "hop_tidx": ("hop_time", []),
                "hop_from": ("hop", []),
                "hop_to": ("hop", []),
                "trajid": ("hop", []),
                "time": empty_2d,
                "tidx": empty_2d,
            }
        )
        return res

    to_cat = []
    trajids = []
    for (trajid, time), hop in hop_vals.groupby("frame"):
        traj = sel_trajs(frames, trajid)
        orig_time = traj["time"].data
        hop_time = traj.time - time
        hop_time = hop_time.swap_dims({"frame": "hop_time"})
        hop_time = hop_time.assign_coords(hop_time=hop_time).drop_vars(
            ["frame", "trajid", "time"]
        )
        traj = traj.swap_dims({"frame": "hop_time"})
        traj = traj.assign_coords(hop_time=hop_time).drop_vars(
            ["frame", "trajid", "time"]
        )

        # Add per-hop metadata
        traj = traj.assign_coords(time=(("hop", "hop_time"), orig_time[None, :]))
        tidx = xr.Variable(dims=("hop_time"), data=np.arange(len(orig_time)))
        traj = traj.assign_coords(tidx=tidx.expand_dims("hop"))

        # Add further hop-independent metadata
        traj = traj.assign_coords(hop_tidx=tidx - hop["tidx"].item())

        traj = traj.drop_dims(["trajid_"], errors="ignore")
        if window is not None:
            traj = traj.sel(hop_time=window)

        trajids.append(trajid)
        to_cat.append(traj)

    # FIXME @thevro: xarray 2025.12.0 FutureWarning: data_vars = 'all'->None
    # FIXME @thevro: xarray 2025.12.0 FutureWarning: coords = 'different'->'minimal'
    res = xr.concat(to_cat, "hop", join="outer")
    from_to = (
        hop_vals[["hop_from", "hop_to"]]
        .drop_vars(["frame", "trajid", "time", "tidx"])
        .rename({"frame": "hop"})
    )
    res = res.assign_coords(
        trajid=("hop", trajids), hop_from=from_to["hop_from"], hop_to=from_to["hop_to"]
    )
    return res


@overload
def assign_hop_time(
    frames: xr.Dataset | xr.DataArray | DataSeries,
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
    which: Literal["first", "last"] = "last",
) -> xr.Dataset | xr.DataArray | DataSeries: ...


@overload
def assign_hop_time(
    frames: TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
    which: Literal["first", "last"] = "last",
) -> TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries]: ...


def assign_hop_time(
    frames: xr.Dataset
    | xr.DataArray
    | DataSeries
    | TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries],
    hop_type_selection: StateSelection | StateSelectionDescriptor | None = None,
    which: Literal["first", "last"] = "last",
) -> (
    xr.Dataset
    | xr.DataArray
    | DataSeries
    | TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries]
):
    """Assign a ``hop_time`` coordinate along the ``frames`` axis giving times
    relative to hops

    Parameters
    ----------
    frames : xr.Dataset | xr.DataArray | DataSeries | TreeNode[Any, xr.Dataset | xr.DataArray | DataSeries]
        An Xarray object (Dataset or DataArray) with a ``frames`` dimension
    hop_type_selection : StateSelection | StateSelectionDescriptor, optional
        Types of hops to include
        See like-named parameter in :py:func:`shnitsel.analyze.hops.hops`
    which : {"first", "last"}, default="last"
        Which hop to take, in case multiple hops present within a single trajectory;
        either 'first' or 'last' (the default)

    Returns
    -------
    The same type as the ``frames`` object (potentially as a wrapped dataset) with --

        - the ``hop_time`` coordinate added along the ``frame`` or ``time`` (i.e. leading dim) dimension, containing
          all times relative to one chosen hop in each trajectory,
        - the ``time_at_hop`` coordinate added along the ``trajectory`` dimension,
          containing the time at which each chosen hop occurred

    Both of these coordinates contain ``nan`` for trajectories lacking any hops of the
    types specified
    """
    if isinstance(frames, TreeNode):
        return frames.map_data(
            assign_hop_time,
            hop_type_selection=hop_type_selection,
            which=which,
        )

    leading_dim: DimName

    if 'time' not in frames.coords:
        logging.error(
            "Cannot assign hop time for dataset or array without `time` information."
        )
        raise ValueError("No `time` coordinate provided to `assign_hop_time`.")

    # We cannot add arbitrary coordinates to a data array, so we have to deal with it separately:
    if isinstance(frames, xr.DataArray):
        if 'time' in frames.dims:
            leading_dim = 'time'
        elif 'frame' in frames.dims:
            leading_dim = 'frame'
        else:
            raise ValueError(
                "Could not determine leading dimension of `frames`. Was neither `time` nor `frame`."
            )

        if frames.sizes[leading_dim] == 0:
            return frames.assign_coords(hop_time=(leading_dim, []))
    else:
        wrapped_ds = wrap_dataset(frames, DataSeries)
        leading_dim = wrapped_ds.leading_dim

        # TODO: FIXME: Refactor this to new wrapper types
        if frames.sizes[leading_dim] == 0:
            return frames.assign_coords(hop_time=(leading_dim, []))

    hop_data = filter_data_at_hops(
        frames, hop_type_selection=hop_type_selection
    )  # .reset_index("frame")
    if hop_data.sizes[leading_dim] == 0:
        # TODO: FIXME: This return does not match the description in the `Returns` docstring block. We state there, that we assign nan if no value was found, but we simply skip the trajectory.
        return frames.assign_coords(
            hop_time=(leading_dim, np.full(frames.sizes[leading_dim], np.nan))
        )

    if which == "first":
        fn = min
    elif which == "last":
        fn = max

    time_attrs = dict(frames.coords['time'].attrs)

    time_at_hop_attrs = dict(time_attrs)
    time_at_hop_attrs["long_name"] = "Time at which selected hop happened"
    hop_time_attrs = dict(time_attrs)
    hop_time_attrs["long_name"] = "Time relative to when the selected hop happened"

    # Need to assign per-trajectory data only if we have multiple trajectories
    if 'trajectory' not in frames.dims and 'atrajectory' not in frames.coords:
        # We have a single trajectory:

        # We must have at least one hit or `hop_data` would not have had a hit in the `leading_dim`
        time_at_hop = fn(frames.coords["time"]).item()

        hop_time = (frames.time - time_at_hop).assign_attrs(**hop_time_attrs)

        # Add scalar `time_at_hop` coordinate:
        time_at_hop_da = xr.DataArray(
            time_at_hop, dims=(), name='time_at_hop', attrs=hop_time_attrs
        )
        return frames.assign_coords(time_at_hop=time_at_hop_da, hop_time=hop_time)
    else:
        # If we end up here, we should have a multi-trajectory set.
        if leading_dim == 'frame' and 'atrajectory' in frames.coords:
            # We have a stacked trajectory set
            d_times = {
                trajid: fn(traj.coords["time"]).item()
                for trajid, traj in hop_data.groupby("atrajectory")
            }

            hop_time = (
                frames.time.groupby("atrajectory")
                .map(
                    lambda traj: (
                        traj.time.data
                        - d_times.get(traj['atrajectory'].item(0), np.nan)
                    )
                )
                .assign_attrs(**hop_time_attrs)
            )
        elif 'trajectory' in frames.dims:
            # We have a layered trajectory set
            d_times = {
                trajid: fn(traj.coords["time"]).item()
                for trajid, traj in hop_data.groupby("trajectory")
            }

            hop_time = (
                frames.time.groupby("trajectory")
                .map(
                    lambda traj: (
                        traj.time.data - d_times.get(traj['trajectory'].item(0), np.nan)
                    )
                )
                .assign_attrs(**hop_time_attrs)
            )
        else:
            raise ValueError(
                "Unknown trajectory format: Trajectory could not be identified as either a stacked or layered multi-trajectory format."
            )

        if 'trajectory' in frames.dims:
            do_not_add_coord = False
            if 'trajectory' not in frames.coords:
                if isinstance(frames, xr.DataArray):
                    # Cannot add new dimension to array
                    logging.debug(
                        "Cannot assign new `time_at_hop` coordinate to xr.DataArray with missing `trajectory` dimension."
                    )
                    do_not_add_coord = True
                elif 'atrajectory' in frames.coords:
                    # Can recreate `trajectory` dimension if active trajectory information is present
                    frames = frames.assign_coords(
                        trajectory=np.unique(frames.coords['atrajectory'].data)
                    )
                else:
                    # Cannot recreate `trajectory` dimension if no trajectory data whatsoever is present
                    raise ValueError(
                        "Data is lacking `trajectory` coordinate and there is not enough information to reconstruct the `trajectory` dimension."
                    )

            if not do_not_add_coord:
                # We are not a data array and we have the `trajectory` dimension initialized.
                time_at_hop = [
                    d_times.get(trajid.item(), np.nan)
                    for trajid in frames.coords["trajectory"]
                ]
                frames = frames.assign_coords(
                    time_at_hop=("trajectory", time_at_hop, time_at_hop_attrs)
                )

        return frames.assign_coords(hop_time=hop_time)

    # TODO (thevro): When the data of an input DataArray doesn't have a 'trajectory' dimension,
    # we can't add the `time_at_hop` info as a coordinate. Alternative approaches?
    # It's not exactly hard to do `da.sel(hop_time=0).time` anyway, which gives the same info.
    # if 'trajectory' in frames.dims:
    #     if 'trajectory' not in frames.coords:
    #         frames = frames.assign_coords(
    #             trajectory=np.unique(frames.coords['atrajectory'].data)
    #         )
    #     time_at_hop = [
    #         d_times.get(trajid.item(), np.nan) for trajid in frames.coords["trajectory"]
    #     ]
    #     frames = frames.assign_coords(time_at_hop=("trajectory", time_at_hop))

    # return frames.assign_coords(hop_time=hop_time)
