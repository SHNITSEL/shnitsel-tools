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
        if isinstance(active_state_source, ShnitselDataset):
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
            # if hop_types is not None:
            #     acc = np.full(res.sizes['frame'], False)
            #     for hop_from, hop_to in hop_types:
            #         acc |= (res.hop_from == hop_from) & (res.hop_to == hop_to)
            #     res = res.isel(frame=acc)
            # return res.drop_dims(['trajid_'], errors='ignore')

            # We drop the mask that should be all True values now.
            return wrap_dataset(
                res_dataset.drop("is_hop_mask", errors="ignore"), DataSeries
            )


# TODO: FIXME: Make StateSelection the preferred type for picking hopping types.
def focus_hops(
    frames, hop_types: list[tuple[int, int]] | None = None, window: slice | None = None
):
    """For each hop, create a copy of its trajectory centered on the hop; align these

    Parameters
    ----------
    frames
        An Xarray object (Dataset or DataArray) with a ``frames`` dimension
    hop_types
        Types of hops to include
        See like-named parameter in :py:func:`shnitsel.analyze.hops.hops`
    window
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
    hop_vals = hops(frames, hop_types=hop_types)
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


# TODO: FIXME: Make StateSelection the preferred type for picking hopping types.
def assign_hop_time(
    frames,
    hop_types: list[tuple[int, int]] | None = None,
    which: Literal["first", "last"] = "last",
):
    """Assign a ``hop_time`` coordinate along the ``frames`` axis giving times
    relative to hops

    Parameters
    ----------
    frames
        An Xarray object (Dataset or DataArray) with a ``frames`` dimension
    hop_types
        Types of hops to include
        See like-named parameter in :py:func:`shnitsel.analyze.hops.hops`
    which
        Which hop to take, in case multiple hops present within a single trajectory;
        either 'first' or 'last' (the default)

    Returns
    -------
    The same ``frames`` object with --

        - the ``hop_time`` coordinate added along the ``frame`` dimension, containing
          all times relative to one chosen hop in each trajectory,
        - the ``time_at_hop`` coordinate added along the ``trajectory`` dimension,
          containing the time at which each chosen hop occurred

    Both of these coordinates contain ``nan`` for trajectories lacking any hops of the
    types specified
    """
    if isinstance(frames, (xr.DataTree, TreeNode, ShnitselDB)):
        raise NotImplementedError(
            f"This function is not yet implemented for type {type(frames)}"
        )
    # TODO: FIXME: Refactor this to new wrapper types
    if frames.sizes["frame"] == 0:
        return frames.assign_coords(hop_time=("frame", []))

    hop_vals = hops(frames, hop_types=hop_types).reset_index("frame")
    if hop_vals.sizes["frame"] == 0:
        return frames.assign_coords(
            hop_time=("frame", np.full(frames.sizes["frame"], np.nan))
        )

    if which == "first":
        fn = min
    elif which == "last":
        fn = max
    d_times = {
        trajid: fn(traj.coords["time"]).item()
        for trajid, traj in hop_vals.groupby("atrajectory")
    }

    hop_time = frames.time.groupby("atrajectory").map(
        lambda traj: traj.time.data - d_times.get(traj['atrajectory'].item(0), np.nan)
    )

    # TODO (thevro): When the data of an input DataArray doesn't have a 'trajectory' dimension,
    # we can't add the `time_at_hop` info as a coordinate. Alternative approaches?
    # It's not exactly hard to do `da.sel(hop_time=0).time` anyway, which gives the same info.
    if 'trajectory' in frames.dims:
        if 'trajectory' not in frames.coords:
            frames = frames.assign_coords(
                trajectory=np.unique(frames.coords['atrajectory'].data)
            )
        time_at_hop = [
            d_times.get(trajid.item(), np.nan) for trajid in frames.coords["trajectory"]
        ]
        frames = frames.assign_coords(time_at_hop=("trajectory", time_at_hop))

    return frames.assign_coords(hop_time=hop_time)
