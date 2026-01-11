from typing import Literal

import xarray as xr
import numpy as np

from shnitsel.data.multi_indices import mdiff, sel_trajs


# TODO: type-hinting of appropriate generality for first argument
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

    Returns
    -------
    An indexed version of ``frames``, where each entry in the
    ``frames`` dimension represents a hop.
    The following coordinates are added along ``frames``:

        - ``tidx``: the time-step index of the hop in its trajectory
        - ``hop_from``: the active state before the hop
        - ``hop_to``: the active state after the hop
    """
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
    return res.drop_dims(['trajid_'], errors='ignore')


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
    hop_vals = hops(frames, hop_types=hop_types)
    # If no hops, return empty
    if hop_vals.sizes['frame'] == 0:
        res = frames.isel(frame=[])
        res = res.swap_dims({'frame': 'hop_time'})
        res = res.drop_vars(['frame', 'trajid', 'time'])
        res = res.drop_dims(['trajid_'], errors='ignore')
        res = res.expand_dims('hop').isel(hop=[])
        empty_2d = xr.Variable(('hop', 'hop_time'), [[]]).isel(hop=[], hop_time=[])
        res = res.assign_coords(
            {
                'hop_time': ('hop_time', []),
                'hop_tidx': ('hop_time', []),
                'hop_from': ('hop', []),
                'hop_to': ('hop', []),
                'trajid': ('hop', []),
                'time': empty_2d,
                'tidx': empty_2d,
            }
        )
        return res

    to_cat = []
    trajids = []
    for (trajid, time), hop in hop_vals.groupby('frame'):
        traj = sel_trajs(frames, trajid)
        orig_time = traj['time'].data
        hop_time = traj.time - time
        hop_time = hop_time.swap_dims({'frame': 'hop_time'})
        hop_time = hop_time.assign_coords(hop_time=hop_time).drop_vars(
            ['frame', 'trajid', 'time']
        )
        traj = traj.swap_dims({'frame': 'hop_time'})
        traj = traj.assign_coords(hop_time=hop_time).drop_vars(
            ['frame', 'trajid', 'time']
        )

        # Add per-hop metadata
        traj = traj.assign_coords(time=(('hop', 'hop_time'), orig_time[None, :]))
        tidx = xr.Variable(dims=('hop_time'), data=np.arange(len(orig_time)))
        traj = traj.assign_coords(tidx=tidx.expand_dims('hop'))

        # Add further hop-independent metadata
        traj = traj.assign_coords(hop_tidx=tidx - hop['tidx'].item())

        traj = traj.drop_dims(['trajid_'], errors='ignore')
        if window is not None:
            traj = traj.sel(hop_time=window)

        trajids.append(trajid)
        to_cat.append(traj)

    # FIXME @thevro: xarray 2025.12.0 FutureWarning: data_vars = 'all'->None
    # FIXME @thevro: xarray 2025.12.0 FutureWarning: coords = 'different'->'minimal'
    res = xr.concat(to_cat, 'hop', join='outer')
    from_to = (
        hop_vals[['hop_from', 'hop_to']]
        .drop_vars(['frame', 'trajid', 'time', 'tidx'])
        .rename({'frame': 'hop'})
    )
    res = res.assign_coords(
        trajid=('hop', trajids), hop_from=from_to['hop_from'], hop_to=from_to['hop_to']
    )
    return res


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
    The same ``frames`` object with --

        - the ``hop_time`` coordinate added along the ``frame`` dimension, containing
          all times relative to one chosen hop in each trajectory,
        - the ``time_at_hop`` coordinate added along the ``trajid_`` dimension,
          containing the time at which each chosen hop occurred

    Both of these coordinates contain ``nan`` for trajectories lacking any hops of the
    types specified
    """
    if frames.sizes['frame'] == 0:
        return frames.assign_coords(hop_time=('frame', []))

    hop_vals = hops(frames, hop_types=hop_types).reset_index('frame')
    if hop_vals.sizes['frame'] == 0:
        return frames.assign_coords(
            hop_time=('frame', np.full(frames.sizes['frame'], np.nan))
        )

    if which == 'first':
        fn = min
    elif which == 'last':
        fn = max
    d_times = {
        trajid: fn(traj.coords['time']).item()
        for trajid, traj in hop_vals.groupby('trajid')
    }

    hop_time = frames.time.groupby('trajid').map(
        lambda traj: traj.time.data - d_times.get(traj.trajid.item(0), np.nan)
    )

    time_at_hop = [
        d_times.get(trajid.item(), np.nan) for trajid in frames.coords['trajid_']
    ]

    return frames.assign_coords(hop_time=hop_time, time_at_hop=('trajid_', time_at_hop))