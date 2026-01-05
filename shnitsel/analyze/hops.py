from typing import Literal, TypeVar, overload

import xarray as xr
import numpy as np

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.multi_indices import mdiff
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.state_selection import StateSelection

TrajectoryOrFrames = TypeVar("TrajectoryOrFrames", bound=Trajectory | Frames)


# TODO: Finish documentation
@overload
def hops_mask_from_active_state(
    active_state_source: Trajectory | Frames | xr.DataArray,
    hop_type_selection: StateSelection | None = None,
) -> xr.DataArray:
    """Overload to specify simple return type for simple (flat) input types.

    See `hops_mask_from_active_state()` for details

    Parameters
    ----------
    active_state_source : Trajectory | Frames | xr.DataArray
        The simple object to extract active state information and calculate hopping points from.
        From Trajectory or Frames objects, extracts the `active_state` (`astate`) variable.
        plain DataArray objects are assumed to hold integer information on the active state.
    hop_type_selection: StateSelection, optional
        A state selection holding state transitions that should be used in hop filtering.

    Returns
    -------
    xr.DataArray
        The boolean mask of points along the leading dimension where hops occur.
        Has the following coordinates attached:
            - ``hop_from``: The state from which the hop at this point would occur
            - ``hop_to``: The state to which the hop would be.
        The flag is set in the frame along the leading dimension after the frame boundary where
        the active state changes, meaning hops may occur in the last frame of a trajectory but
        never in the first frame of a trajectory.

    Raises
    ------
    ValueError
        Unsupported input type
    """
    ...


@overload
def hops_mask_from_active_state(
    active_state_source: ShnitselDB[Trajectory | Frames | xr.DataArray],
    hop_type_selection: StateSelection | None = None,
) -> ShnitselDB[xr.DataArray]:
    """Overload to specify hierarchical return type for hierarchical input types.

    See other overloads of `hops_mask_from_active_state` for details of how the individual data entries are
    mapped

    Parameters
    ----------
    active_state_source : ShnitselDB[Trajectory | Frames | xr.DataArray]
        A hierarchical set of sources of active state data that should be mapped using `hops_mask_from_active_state`
    hop_type_selection: StateSelection, optional
        A state selection holding state transitions that should be used in hop filtering.

    Returns
    -------
    ShnitselDB[xr.DataArray]
        The hierarchically structured result of running the mapping of `hops_mask_from_active_state()` over entries in the
        input tree.

    Raises
    ------
    ValueError
        Invalid input type in one of the data leaves of the tree.
    """
    ...


def hops_mask_from_active_state(
    active_state_source: Trajectory
    | Frames
    | ShnitselDB[Trajectory | Frames | xr.DataArray]
    | xr.DataArray,
    hop_type_selection: StateSelection | None = None,
) -> xr.DataArray | ShnitselDB[xr.DataArray]:
    """Generate boolean masks marking hopping points by identifying changes in the active state of provided
    data source.

    Needs to be fed either with (hierarchical) trajectory data that has `active_state` (`astate`) information
    or directly with the xr.DataArray holding `astate` information.

    Parameters
    ----------
    active_state_source : Trajectory | Frames | ShnitselDB[Trajectory  |  Frames  |  xr.DataArray] | xr.DataArray
        A potential source for extracting the active state along a leading dimension and the leading dimension name.
    hop_type_selection: StateSelection, optional
        A state selection holding state transitions that should be used in hop filtering.

    Returns
    -------
    xr.DataArray | ShnitselDB[xr.DataArray]
        Either the flat boolean mask of leading dimension instances where a hop happens or a hierarchical structure holding
        such a flat mask for every original data entry in the hierarchical input data

    Raises
    ------
    ValueError
        If an unsupported input type was provided
    """
    if isinstance(active_state_source, ShnitselDB):
        return active_state_source.map_data(hops_mask_from_active_state)

    else:
        active_state_data: xr.DataArray
        leading_dim: str | None
        if isinstance(active_state_source, Trajectory) or isinstance(
            active_state_source, Frames
        ):
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


@overload
def filter_data_at_hops(
    active_state_and_data_source: TrajectoryOrFrames,
    hop_type_selection: StateSelection | None = None,
) -> TrajectoryOrFrames:
    """Overload of `filter_data_at_hops` for non-hierarchical input types to indicate non-hierarchical returns.

    _extended_summary_

    Parameters
    ----------
    active_state_and_data_source : Trajectory | Frames
        _description_
    hop_type_selection : StateSelection | None, optional
        _description_, by default None

    Returns
    -------
    Trajectory | Frames
        _description_
    """
    ...


@overload
def filter_data_at_hops(
    active_state_and_data_source: ShnitselDB[TrajectoryOrFrames],
    hop_type_selection: StateSelection | None = None,
) -> ShnitselDB[TrajectoryOrFrames]:
    """Overload of `filter_data_at_hops` for hierarchical input types to indicate hierarchical returns.

    _extended_summary_

    Parameters
    ----------
    active_state_and_data_source : ShnitselDB[TrajectoryOrFrames]
        _description_
    hop_type_selection : StateSelection | None, optional
        _description_, by default None

    Returns
    -------
    ShnitselDB[TrajectoryOrFrames]
        _description_
    """


def filter_data_at_hops(
    active_state_and_data_source: Trajectory | Frames | ShnitselDB[Trajectory | Frames],
    hop_type_selection: StateSelection | None = None,
) -> Trajectory | Frames | ShnitselDB[Trajectory | Frames]:
    """Filter data to only retain data at points where hops of selected transitions occur.

    Needs to be fed either with (hierarchical) trajectory data that has `active_state` (`astate`) information
    or with simple (single) trajectory data with an `active_state` (`astate`) variable.

    If you wish to perform arbitrary filtering, you can employ the `hops_mask_from_active_state()` function
    to just get the boolean mask of hop positions and perform the filtering yourself.

    Parameters
    ----------
    active_state_and_data_source : Trajectory | Frames | ShnitselDB[Trajectory  |  Frames  |  xr.DataArray] | xr.DataArray
        A source for extracting the active state along a leading dimension and the leading dimension name as well as to filter
        the data from.
    hop_type_selection: StateSelection, optional
        A state selection holding state transitions that should be used in hop filtering.

    Returns
    -------
    Trajectory | Frames | ShnitselDB[Trajectory | Frames]
        Filtered version of the input data source, where a selection was performed to only retain data at points where
        a hop was happening.
        Each entry along the leading dimension (`time` or `frame`) represents a hop.
        The following coordinates are added along the leading dimension:
            - ``hop_from``: the active state before the hop
            - ``hop_to``: the active state after the hop
    """

    if isinstance(active_state_and_data_source, ShnitselDB):
        return active_state_and_data_source.map_data(
            lambda x: filter_data_at_hops(x, hop_type_selection=hop_type_selection)
        )
    else:
        # Frames or Trajectory
        input_dataset = active_state_and_data_source.dataset
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
        return type(active_state_and_data_source)(
            res_dataset.drop('is_hop_mask', errors='ignore')
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
    The following coordinates are added along dimension ``hop_time``:

        - ``hop_time``: the trajectory time coordinate relative to the hop
        - ``time``: the original time coordinate relative to the start of the trajectory
        - ``hop_tidx``: the trajectory time-step index relative to the hop
        - ``tidx``: the trajectory time-step index relative to the start of the trajectory

    The following coordinates are added along dimension ``hop``:

        - ``hop_from``: the active state before the hop
        - ``hop_to``: the active state after the hop
    """
    raise NotImplementedError()
    # TODO: FIXME: Refactor this to new wrapper types
    hop_vals = hops(frames, hop_types=hop_types)
    to_cat = []
    trajids = []
    for (trajid, time), hop in hop_vals.groupby('frame'):
        traj = frames.st.sel_trajs(trajid)
        orig_time = traj.time
        hop_time = traj.time - time
        hop_time = hop_time.swap_dims({'frame': 'hop_time'})
        hop_time = hop_time.assign_coords(hop_time=hop_time).drop_vars(
            ['frame', 'trajid', 'time']
        )
        traj = traj.swap_dims({'frame': 'hop_time'})
        traj = traj.assign_coords(hop_time=hop_time).drop_vars(
            ['frame', 'trajid', 'time']
        )

        # Add useful metadata
        traj = traj.assign_coords(time=('hop_time', orig_time.data))
        traj = traj.assign_coords(tidx=('hop_time', np.arange(traj.sizes['hop_time'])))
        traj = traj.assign_coords(hop_tidx=traj['tidx'] - hop['tidx'].item())

        trajids.append(traj.coords['trajid_'].item())
        traj = traj.drop_dims('trajid_')
        if window is not None:
            traj = traj.sel(hop_time=window)

        to_cat.append(traj)

    res = xr.concat(to_cat, 'hop', join='outer')
    from_to = (
        hop_vals[['hop_from', 'hop_to']]
        .drop_vars(['frame', 'trajid', 'time'])
        .rename({'frame': 'hop'})
    )
    res = res.assign_coords(
        trajid=('hop', trajids), hop_from=from_to['hop_from'], hop_to=from_to['hop_to']
    )
    return res


# TODO: FIXME: Make StateSelection the preferred type for picking hopping types.
def assign_hop_time(
    frames,
    hop_types: list[tuple[int, int]] | None = None,
    which: Literal['first', 'last'] = 'last',
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
    The same ``frames`` object with the ``hop_time`` coordinate added along the
    ``frame`` dimension, containing times relative to one chosen hop in each
    trajectory, and containing ``nan`` in trajectories lacking any hops of the
    types specified.
    """
    raise NotImplementedError()
    # TODO: FIXME: Refactor this to new wrapper types
    hop_vals = hops(frames, hop_types=hop_types).reset_index('frame')

    if which == 'first':
        fn = min
    elif which == 'last':
        fn = max
    hop_time = {
        trajid: fn(traj.coords['time']).item()
        for trajid, traj in hop_vals.groupby('trajid')
    }

    hop_time = frames.time.groupby('trajid').map(
        lambda traj: traj.time.data - hop_time.get(traj.trajid.item(0), np.nan)
    )

    return frames.assign_coords(hop_time=hop_time)
