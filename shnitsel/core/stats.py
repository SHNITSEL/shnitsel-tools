from typing import Hashable, TypeAlias

import numpy as np
import numpy.typing as npt
import scipy.stats as st
import xarray as xr

from .._contacts import needs

DimName: TypeAlias = Hashable


#####################################################
# For calculating confidence intervals, the following
# functions offer varying levels of abstraction
# TODO make naming consistent


def calc_ci(a: npt.NDArray, confidence: float = 0.95) -> npt.NDArray:
    if np.array(a).ndim != 1:
        raise ValueError("This function accepts 1D input only")
    return np.stack(
        st.t.interval(confidence, len(a) - 1, loc=np.mean(a), scale=st.sem(a))
    )


def ci_agg_last_dim(a, confidence=0.95):
    outer_shape = tuple(a.shape[:-1])
    res = np.full(outer_shape + (3,), np.nan)
    for idxs in np.ndindex(outer_shape):
        res[idxs, :2] = calc_ci(a[idxs], confidence=confidence)
        res[idxs, 2] = np.mean(a[idxs])
    return res


def xr_calc_ci(a: xr.DataArray, dim: DimName, confidence: float = 0.95) -> xr.Dataset:
    res_da: xr.DataArray = xr.apply_ufunc(
        ci_agg_last_dim,
        a,
        kwargs={'confidence': confidence},
        output_core_dims=[['bound']],
        input_core_dims=[[dim]],
    )
    return res_da.assign_coords(  #
        dict(bound=['lower', 'upper', 'mean'])
    ).to_dataset('bound')


@needs(groupable={'time'}, dims={'frame'})
def time_grouped_ci(x: xr.DataArray, confidence: float = 0.9) -> xr.Dataset:
    return x.groupby('time').map(
        lambda x: xr_calc_ci(x, dim='frame', confidence=confidence)
    )