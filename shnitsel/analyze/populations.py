from logging import warning

import numpy as np
import xarray as xr

from shnitsel.core.typedefs import Frames

from .._contracts import needs


@needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
def calc_classical_populations(frames: Frames) -> xr.DataArray:
    """Function to calculate classical state populations from the active state information in `astate` of the dataset `frames.

    Does not use the partial QM coefficients of the states.

    Args:
        frames (Frames): The dataset holding the active state information in a variable `astate`.

    Returns:
        xr.DataArray: The array holding the ratio of trajectories in each respective state.
    """
    # TODO: FIXME: Make this able to deal with ShnitselDB/tree data directly. This should not be too much of an issue?
    data = frames['astate']
    if -1 in frames['astate']:
        warning(
            "`frames['astate']` contains the placeholder value `-1`, "
            "indicating missing state information.  "
            "The frames in question will be excluded from the "
            "population count altogether."
        )
        data = data.sel(frame=(data != -1))
    nstates = frames.sizes['state']
    # zero_or_one = int(frames.coords['state'].min())
    lowest_state_id = 1  # TODO: For now, assume lowest state is 1
    assert lowest_state_id in {0, 1}
    if 'frame' in data.dims:
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
    elif 'time' in data.dims:
        num_ts = frames.sizes['time']
        num_state = frames.sizes['state']

        populations = np.zeros((num_ts, num_state), dtype=np.float32)
        for a in range(frames.sizes['time']):
            if frames.astate[a] > 0:
                populations[a][frames.astate[a] - 1] = 1.0

        pops = xr.DataArray(
            populations, dims=['time', 'state'], coords={'time': data.coords['time']}
        )
    else:
        raise ValueError("Neither time nor frame dimension in dataset.")
    return (
        (pops / pops.sum('state'))
        .assign_coords(state=frames['state'])
        .assign_attrs({'long_name': "Population of different states over times"})
    )


# Alternative name of the function to ca
calc_pops = calc_classical_populations
