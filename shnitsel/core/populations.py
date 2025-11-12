from logging import warning
from typing import TypeAlias

import numpy as np
import xarray as xr

from .._contracts import needs

Frames: TypeAlias = xr.Dataset


@needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
def classical(frames: Frames) -> xr.DataArray:
    """Fast way to calculate populations
    Requires states ids to be small integers
    """
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
    zero_or_one = 1  # TODO: For now, assume lowest state is 1
    assert zero_or_one in {0,1}
    pops = data.groupby('time').map(
        lambda group: xr.apply_ufunc(
            lambda values: np.bincount(values, minlength=nstates + zero_or_one)[
                zero_or_one:
            ],
            group,
            input_core_dims=[['frame']],
            output_core_dims=[['state']],
        )
    )
    return (pops / pops.sum('state')).assign_coords(state=frames['state'])

calc_pops = classical