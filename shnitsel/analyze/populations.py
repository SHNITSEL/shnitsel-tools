import logging

import numpy as np
import xarray as xr

from shnitsel.core._api_info import API, internal
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.units.definitions import time, unit_dimensions

from .._contracts import needs


@API()
@needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
def calc_classical_populations(
    frames: ShnitselDB[Trajectory | Frames] | Trajectory | Frames | xr.Dataset,
) -> xr.DataArray:
    """Function to calculate classical state populations from the active state information in `astate` of the dataset `frames.

    Does not use the partial QM coefficients of the states.

    Args:
        frames (Frames): The dataset holding the active state information in a variable `astate`.

    Returns:
        xr.DataArray: The array holding the ratio of trajectories in each respective state.
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
        pass
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
        xr.DataArray: The array holding the ratio of trajectories in each respective state.
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
    # TODO: FIXME: Do we want absolute populations as well?
    return (
        (pops / pops.sum('state'))
        .assign_coords(state=frames.state_ids)
        .assign_attrs({'long_name': "Population of different states over times"})
    )


# Alternative name of the function to calculate population statistics
calc_pops = calc_classical_populations
