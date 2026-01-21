from typing import Any, Literal, TypeVar, overload
import re

import xarray as xr
import numpy as np

from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.data_series import DataSeries
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.multi_indices import mdiff, sel_trajs
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.state_selection import StateSelection


def _standard_hop_spec(spec):
    # TODO: FIXME: move to state_selection
    search = re.compile("(?P<state_from>.+)(?P<rel>(->)|(<>))(?P<state_to>.+)")
    if not isinstance(spec, str):
        return spec

    subs = re.split(r"\s*,\s*", spec)

    res = []
    for sub in subs:
        found = search.match(sub)
        state_from = int(found.group("state_from"))
        state_to = int(found.group("state_to"))
        rel = found.group("rel")
        if rel == "->":
            res.append((state_from, state_to))
        else:
            res.extend([(state_from, state_to), (state_to, state_from)])
    return res


TrajectoryOrFrames = TypeVar("TrajectoryOrFrames", bound=DataSeries)


# TODO: Finish documentation
@overload
def hops_mask_from_active_state(
    active_state_source: xr.Dataset | DataSeries | xr.DataArray,
    hop_type_selection: StateSelection | None = None,
) -> xr.DataArray:
    """Overload to specify simple return type for simple (flat) input types.

    See `hops_mask_from_active_state()` for details
    """
    ...


@overload
def hops_mask_from_active_state(
    active_state_source: TreeNode[Any, xr.Dataset | DataSeries | xr.DataArray],
    hop_type_selection: StateSelection | None = None,
) -> TreeNode[Any, xr.DataArray]:
    """Overload to specify hierarchical return type for hierarchical input types.

    See other overloads of `hops_mask_from_active_state` for details of how the individual data entries are
    mapped.
    """
    ...


def hops_mask_from_active_state(
    active_state_source: xr.Dataset
    | DataSeries
    | TreeNode[Any, xr.Dataset | DataSeries | xr.DataArray]
    | xr.DataArray,
    hop_type_selection: StateSelection | None = None,
) -> xr.DataArray | TreeNode[Any, xr.DataArray]:
    """Generate boolean masks marking hopping points by identifying changes in the active state of provided
    data source.

    Needs to be fed either with (hierarchical) trajectory data that has `active_state` (`astate`) information
    or directly with the xr.DataArray holding `astate` information.

    Parameters
    ----------
    active_state_source : xr.Dataset | Trajectory | Frames | xr.DataArray | TreeNode[Any, xr.Dataset | Trajectory  |  Frames  |  xr.DataArray]
        A potential source for extracting the active state along a leading dimension and the leading dimension name.
    hop_type_selection: StateSelection, optional
        A state selection holding state transitions that should be used in hop filtering.

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
            active_state_data = active_state_source
            leading_dim = str(active_state_source.dims[0])
        else:
            raise ValueError(
                "Unknown type of provided source for `active_state` data: %s"
                % type(active_state_source)
            )

        is_hop_mask = mdiff(active_state_data, dim=leading_dim) != 0
        # Add prior and current state back as coordinates
        is_hop_mask = is_hop_mask.assign_coords(
            hop_from=(active_state_data.shift({leading_dim: 1}, -1)),
            hop_to=active_state_data,
        )

        if hop_type_selection is not None and leading_dim is not None:
            type_filter = np.full(is_hop_mask.sizes[leading_dim], False)
            for hop_from, hop_to in hop_type_selection.state_combinations:
                type_filter |= (is_hop_mask.hop_from == hop_from) & (
                    is_hop_mask.hop_to == hop_to
                )
            is_hop_mask &= type_filter
        return is_hop_mask


# TODO: FIXME: Make compatible with trees and wrapper datasets
def hops(frames, hop_types: list[tuple[int, int]] | None = None):
    """Select hops

    Parameters
    ----------
    frames
        An Xarray object (Dataset or DataArray) with a ``frames`` dimension
    hop_types
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
    hop_types = _standard_hop_spec(hop_types)
    is_hop = mdiff(frames['astate']) != 0

    res = frames.isel(frame=is_hop)
    tidxs = np.concat(
        [np.arange(traj.sizes['frame']) for _, traj in frames.groupby('trajid')]
    )
    hop_tidx = tidxs[is_hop]
    res = res.assign_coords(
        tidx=('frame', hop_tidx),
        hop_from=(frames['astate'].shift({'frame': 1}, -1).isel(frame=is_hop)),
        hop_to=res['astate'],
    )
    if hop_types is not None:
        acc = np.full(res.sizes['frame'], False)
        for hop_from, hop_to in hop_types:
            acc |= (res.hop_from == hop_from) & (res.hop_to == hop_to)
        res = res.isel(frame=acc)
    if hasattr(res, 'drop_dims'):
        res = res.drop_dims(['trajid_'], errors='ignore')
    return res


@overload
def filter_data_at_hops(
    active_state_and_data_source: xr.Dataset | DataSeries,
    hop_type_selection: StateSelection | None = None,
) -> DataSeries:
    """Overload of `filter_data_at_hops` for non-hierarchical input types to indicate non-hierarchical returns.

    Converts raw datasets into a Shnitsel-Style wrapper
    """
    ...


@overload
def filter_data_at_hops(
    active_state_and_data_source: TreeNode[Any, xr.Dataset | DataSeries],
    hop_type_selection: StateSelection | None = None,
) -> TreeNode[Any, DataSeries]:
    """Overload of `filter_data_at_hops` for hierarchical input types to indicate hierarchical returns.

    Converts raw datasets into a Shnitsel-Style wrapper
    """


def filter_data_at_hops(
    active_state_and_data_source: xr.Dataset
    | DataSeries
    | TreeNode[Any, xr.Dataset | DataSeries],
    hop_type_selection: StateSelection | None = None,
) -> DataSeries | TreeNode[Any, DataSeries]:
    """Filter data to only retain data at points where hops of selected transitions occur.

    Needs to be fed either with (hierarchical) trajectory data that has `active_state` (`astate`) information
    or with simple (single) trajectory data with an `active_state` (`astate`) variable.

    If you wish to perform arbitrary filtering, you can employ the `hops_mask_from_active_state()` function
    to just get the boolean mask of hop positions and perform the filtering yourself.

    Parameters
    ----------
    active_state_and_data_source : xr.Dataset | DataSeries | TreeNode[Any, xr.Dataset | DataSeries]
        A source for extracting the active state along a leading dimension and the leading dimension name as well as to filter
        the data from.
    hop_type_selection: StateSelection, optional
        A state selection holding state transitions that should be used in hop filtering.

    Returns
    -------
    Trajectory | Frames | TreeNode[Any, Trajectory | Frames]
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
        # Frames or Trajectory
        input_dataset = wrap_dataset(active_state_and_data_source, DataSeries)
        is_hop_mask = hops_mask_from_active_state(
            active_state_source=active_state_and_data_source
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
        - the ``time_at_hop`` coordinate added along the ``trajid_`` dimension,
          containing the time at which each chosen hop occurred

    Both of these coordinates contain ``nan`` for trajectories lacking any hops of the
    types specified
    """
    raise NotImplementedError()
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
        for trajid, traj in hop_vals.groupby("trajid")
    }

    hop_time = frames.time.groupby("trajid").map(
        lambda traj: traj.time.data - d_times.get(traj.trajid.item(0), np.nan)
    )

    time_at_hop = [
        d_times.get(trajid.item(), np.nan) for trajid in frames.coords["trajid_"]
    ]

    return frames.assign_coords(hop_time=hop_time, time_at_hop=("trajid_", time_at_hop))
