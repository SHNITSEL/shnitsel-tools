from logging import warning
import logging
from typing import Hashable

import numpy as np
import numpy.typing as npt
import scipy.stats as st
import xarray as xr

from shnitsel.io.shared.trajectory_setup import get_statecomb_coordinate

from .generic import keep_norming, subtract_combinations as subtract_combinations
from .spectra import assign_fosc

from .._contracts import needs

from .typedefs import DimName, Frames, InterState, PerState


#####################################################
# For calculating confidence intervals, the following
# functions offer varying levels of abstraction
# TODO make naming consistent


def calc_confidence_interval(
    data_array: npt.NDArray, confidence: float = 0.95
) -> npt.NDArray:
    """Function to calculate the confidence interval for a variable array `a`.

    The result is a numpy array with stacked entries with the lower and upper limits of the confidence interval.

    Args:
        a (npt.NDArray): The Numpy array to calculate the confidence interval for.
        confidence (float, optional): The confidence level to get the confidence interval for. Defaults to 0.95.

    Raises:
        ValueError: Raised if the provided `data_array` is not one-dimensional

    Returns:
        npt.NDArray: Numpy array with lower and upper bounds of the confidence interval
    """
    if np.array(data_array).ndim != 1:
        raise ValueError("This function accepts 1D input only")
    return np.stack(
        st.t.interval(
            confidence,
            len(data_array) - 1,
            loc=np.mean(data_array),
            scale=st.sem(data_array),
        )
    )


def confidence_interval_aggregate_last_dim(data_array: npt.NDArray, confidence=0.95) -> npt.NDArray:
    """Calculate the confidence interval from statistics aggregated across the last dimension. 

    For our purposes, this should amount to the trajectory being averaged over.

    Args:
        data_array (npt.NDArray): The numpy data array to calculate the confidence interval for. 
        confidence (float, optional): The confidence level to use for calculations. Defaults to 0.95.

    Returns:
        npt.NDArray: A numpy array with (lower_bound,upper_bound,mean) of the confidence interval in the last dimension. Otherwise same shape as data_array.
    """
    outer_shape = tuple(data_array.shape[:-1])
    res = np.full(outer_shape + (3,), np.nan)
    for idxs in np.ndindex(outer_shape):
        res[idxs, :2] = calc_confidence_interval(data_array[idxs], confidence=confidence)
        res[idxs, 2] = np.mean(data_array[idxs])
    return res


def xr_calc_confidence_interval(
    data_array: xr.DataArray, dim: DimName, confidence: float = 0.95
) -> xr.Dataset:
    """Function to calculate confidence interval data for the input data_array. 
    Results are then repackaged back into an xr.Dataset.

    The dimension denoted by `dim` will be aggregated across.

    Args:
        data_array (xr.DataArray): Input data to have confidence intervals calculated for.
        dim (DimName): Dimension to calculate the confidence interval data from.
        confidence (float, optional): Confidence level for Confidence interval calculation. Defaults to 0.95.

    Returns:
        xr.Dataset: Dataset wrapping a Dataarray with coordinate `bound` with values 'lower', 'upper', and 'mean, which refer to the lower and the upper bound of the confidence interval of this and 
    """
    res_da: xr.DataArray = xr.apply_ufunc(
        confidence_interval_aggregate_last_dim,
        data_array,
        kwargs={'confidence': confidence},
        output_core_dims=[['bound']],
        input_core_dims=[[dim]],
    )
    return res_da.assign_coords(  #
        dict(bound=['lower', 'upper', 'mean'])
    ).to_dataset('bound')


@needs(groupable={'time'}, dims={'frame'})
def time_grouped_confidence_interval(
    x: xr.DataArray, confidence: float = 0.9
) -> xr.Dataset:
    # TODO: document
    """Function to calculate the 

    Args:
        x (xr.DataArray): _description_
        confidence (float, optional): _description_. Defaults to 0.9.

    Returns:
        xr.Dataset: _description_
    """
    return x.groupby('time').map(
        lambda x: xr_calc_confidence_interval(x, dim='frame', confidence=confidence)
    )


@needs(dims={'state'})
def get_per_state(frames: Frames) -> PerState:
    props_per = {'energy', 'forces', 'dip_perm'}.intersection(frames.keys())
    per_state = frames[props_per].map(keep_norming, keep_attrs=False)
    per_state['forces'] = per_state['forces'].where(per_state['forces'] != 0)

    per_state['energy'].attrs['long_name'] = r'$E$'
    per_state['forces'].attrs['long_name'] = r'$\mathbf{F}$'
    if 'dip_perm' in per_state:
        per_state['dip_perm'].attrs['long_name'] = r'$\mathbf{\mu}_i$'
    return per_state


@needs(dims={'state'}, coords={'state'})
def get_inter_state(frames: Frames) -> InterState:
    """Calculate inter-state properties of a dataset for certain observables.

    Currently calculates inter-state levels of energy differences.
    Will calculate Differences between the values of these observables indexed by state.
    If no `statecomb` dimension exists, will create one.

    Args:
        frames (Frames): The basis Dataset to calculate the interstate properties for

    Returns:
        InterState: A Dataset containing interstate properties
    """
    prop: Hashable
    inter_base_props = ['energy']
    available_inter_base_props = []
    # TODO: Check that energy is actually the only inter-state property. We already have statecomb for nacs and dip_trans and astate is not state-dependent.
    for prop in inter_base_props:
        if prop in frames:
            available_inter_base_props.append(prop)
        else:
            warning(f"Dataset does not contain variable '{prop}'")

    if 'statecomb' not in frames.sizes:
        logging.info(
            "Creating a new `statecomb` dimension because it was not yet set when calculating inter-state properties."
        )
        frames = frames.assign_coords(get_statecomb_coordinate(frames.state))
        frames['statecomb'].attrs['long_name'] = "State combinations"

    inter_state = frames.copy()
    for prop in available_inter_base_props:
        if 'state' in frames[prop].dims and 'statecomb' not in frames[prop].dims:
            inter_state_res = subtract_combinations(
                frames[prop], dim='state', labels=False
            )
            inter_state[str(prop) + "_interstate"] = inter_state_res
            inter_state[str(prop) + "_interstate"].attrs["long_name"] = (
                f"Derived inter-state differences of variable `{prop}`"
            )

    # TODO: FIXME: We can't just redefine the statecomb dimension. If it is there, we need to keep it.
    # def state_renamer(lo, hi):
    #     if isinstance(lo, int):
    #         lower_str = f"S_{lo-1}"
    #     else:
    #         lower_str = lo
    #     if isinstance(hi, int):
    #         higher_str = f"S_{hi-1}"
    #     else:
    #         higher_str = hi
    #     f'${higher_str} - {lower_str}$'

    # if 'statecomb' in frames:
    #
    #     warning(
    #         "'statecomb' already exists as an index, variable or coordinate"
    #         " in the dataset, hence it will be removed before recomputation"
    #     )
    # inter_state = flatten_midx(inter_state, 'statecomb', state_renamer)
    # inter_state['statecomb'].attrs['long_name'] = "State combinations"

    if {'energy_interstate', 'dip_trans'}.issubset(frames.variables.keys()):
        inter_state = assign_fosc(inter_state)

    return inter_state
