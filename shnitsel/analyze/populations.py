import logging
from typing import Iterable, overload

import numpy as np
import xarray as xr

from shnitsel.core._api_info import API, internal
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.units.definitions import time, unit_dimensions

from .._contracts import needs


@overload
def calc_classical_populations(
    data: ShnitselDB[Trajectory | Frames],
) -> ShnitselDB[xr.DataArray]:
    """Specialized version of the population calculation where a tree hierarchy of Trajectory data
    is mapped to a tree hierarchy of population statistics.

    The hierarchy will first be grouped by metadata, then a population statistics `xr.DataArray` will
    be calculated for each flat group.

    Parameters
    ----------
    data : ShnitselDB[Trajectory  |  Frames]
        The tree-structured trajectory data

    Returns
    -------
    ShnitselDB[xr.DataArray]
        The tree structure holding population statistics for grouped data out of the input.
        Results contain inidividual states' absolute population numbers for every time step.
    """
    ...


@overload
def calc_classical_populations(
    data: Trajectory | Frames | xr.Dataset,
) -> xr.DataArray:
    """Specialized version of the population calculation where a single Trajectory or Frameset instance
    is mapped to population statistics for their different states.

    Parameters
    ----------
    data : Trajectory | Frames | xr.Dataset
        The input dataset to calculate population statistics along the `time` dimension for.

    Returns
    -------
    xr.DataArray
        The multi-dimensional array with coordinates and annotations that holds the absolute population data for states.
    """


@API()
@needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
def calc_classical_populations(
    frames: ShnitselDB[Trajectory | Frames] | Trajectory | Frames | xr.Dataset,
) -> xr.DataArray | ShnitselDB[xr.DataArray]:
    """Function to calculate classical state populations from the active state information in `astate` of the dataset `frames.

    Does not use the partial QM coefficients of the states.

    Args:
        frames (Frames): The dataset holding the active state information in a variable `astate`.

    Returns:
        xr.DataArray: The array holding the ratio of trajectories in each respective state.
        ShnitselDB[xr.DataArray]: The tree holding the hierarchical structure of data to be processed. Will calculate populations across grouped subtrees.
    """
    if isinstance(frames, ShnitselDB):
        # TODO: convert the shnitsel db subtrees to frames and perform some aggregation

        # Calculation for trajcectory:
        # num_ts = frames.sizes['time']
        # num_state = frames.sizes['state']

        # populations = np.zeros((num_ts, num_state), dtype=np.float32)
        # for a in range(frames.sizes['time']):
        #     if frames.astate[a] > 0:
        #         populations[a][frames.astate[a] - 1] = 1.0

        # pops = xr.DataArray(
        #     populations, dims=['time', 'state'], coords={'time': data.coords['time']}
        # )
        db: ShnitselDB[Trajectory | Frames] = frames.group_data_by_metadata()

        def _map_prepare(frames: Trajectory | Frames) -> xr.DataArray:
            if isinstance(frames, Trajectory) and not isinstance(frames, Frames):
                frames = frames.as_frames

            return _calc_classical_populations(frames)

        mapped_pop_data = db.map_data(_map_prepare, keep_empty_branches=False)

        def _combine_func(pop_data: Iterable[xr.DataArray]) -> xr.DataArray:
            res: xr.DataArray | None = None
            res_timelen: int = -1
            for data in pop_data:
                if data is None:
                    continue

                if res is None:
                    # Copy
                    res = data
                    res_timelen = res.sizes['time']
                else:
                    new_timelen = data.sizes['time']

                    time_values: xr.DataArray | None = None

                    if new_timelen < res_timelen:
                        data = data.pad(
                            {'time': (0, res_timelen - new_timelen)},
                            constant_values=0.0,
                            keep_attrs=True,
                        )
                        time_values = res.time
                    else:
                        res = res.pad(
                            {'time': (0, new_timelen - res_timelen)},
                            constant_values=0.0,
                            keep_attrs=True,
                        )
                        time_values = data.time()

                    res = res + data
                    if time_values is not None:
                        res = res.assign_coords(
                            {'time': ('time', time_values, time_values.attrs)}
                        )
            if res is None:
                # No population data, empty array
                return xr.DataArray()
            # Yield the absolute population for this group
            return res

        reduced_pop_data = mapped_pop_data.map_flat_group_data(_combine_func)
        return reduced_pop_data
    else:
        eventual_frames: Frames
        if isinstance(frames, xr.Dataset):
            # Convert to Trajectory or Frames if provided a Dataset
            if 'time' in frames.dims:
                eventual_frames = Trajectory(frames).as_frames
            elif 'frame' in frames.dims:
                eventual_frames = Frames(frames)
            else:
                raise ValueError(
                    "Dataset had neither `time` nor `frame` coordinate. Cannot calculate population statistics"
                )
        elif isinstance(frames, Trajectory):
            eventual_frames = frames.as_frames

        population_results = _calc_classical_populations(eventual_frames)
        return population_results


@internal()
def _calc_classical_populations(
    frames: Frames,
) -> xr.DataArray:
    """Function to calculate classical state populations from the active state information in `astate` of the dataset `frames.

    Does not use the partial QM coefficients of the states.

    Args:
        frames (Frames): The dataset holding the active state information in a variable `astate`.

    Returns:
        xr.DataArray: The array holding the absolute number of trajectories in each respective state.
    """
    # TODO: FIXME: Make this able to deal with ShnitselDB/tree data directly. This should not be too much of an issue?
    data = frames.active_state
    if -1 in data:
        logging.warning(
            "`active_state` data contains the placeholder value `-1`, "
            "indicating missing state information.  "
            "The frames in question will be excluded from the "
            "population count altogether."
        )
        data = data.sel(frame=(data != -1))

    nstates = frames.sizes['state']
    # zero_or_one = int(frames.coords['state'].min())
    lowest_state_id = 1  # TODO: For now, assume lowest state is 1
    assert lowest_state_id in {0, 1}

    # Add dummy time coordinate if missing
    if 'time' not in data.coords:
        data = data.assign_coords(
            {
                'time': (
                    'frame',
                    xr.zeros_like(data.coords['frame']),
                    {'units': time.femto_seconds, 'unit_dim': unit_dimensions.time},
                )
            }
        )
        logging.warning(
            "No `time` coordinate found. Added dummy time coordinate of 0fs to all frames."
        )

    input_core_dims = [['frame']]
    pops = data.groupby('time').map(
        lambda group: xr.apply_ufunc(
            lambda values: np.bincount(values, minlength=nstates + lowest_state_id)[
                lowest_state_id:
            ],
            group,
            input_core_dims=input_core_dims,
            output_core_dims=[['state']],
        )
    )
    pops.name = "abs_population"
    # TODO: FIXME: Do we want absolute populations as well?
    return pops.assign_coords(  # (pops / pops.sum('state'))
        state=frames.state_ids
    ).assign_attrs({'long_name': "Absolute population of different states over times"})


# Alternative name of the function to calculate population statistics
calc_pops = calc_classical_populations
