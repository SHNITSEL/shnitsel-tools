from logging import warning
import logging
from typing import Hashable, TypeAlias

import numpy as np
import numpy.typing as npt
import scipy.stats as st
import xarray as xr

from shnitsel.io.shared.trajectory_setup import get_statecomb_coordinate

from .midx import flatten_midx
from .generic import keep_norming, subtract_combinations as subtract_combinations
from .spectra import assign_fosc

from .._contracts import needs

from .typedefs import DimName, Frames, InterState, PerState


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

    iprops = []
    # TODO: Check that energy is actually the only inter-state property. We already have statecomb for nacs and dip_trans and astate is not state-dependent.
    for prop in ['energy']:  # ['energy', 'nacs', 'astate', 'dip_trans']:
        if prop in frames:
            iprops.append(prop)
        else:
            warning(f"Dataset does not contain variable '{prop}'")

    if 'statecomb' not in frames.sizes:
        logging.info(
            "Creating a new `statecomb` dimension because it was not yet set when calculating inter-state properties."
        )
        frames = frames.assign_coords(get_statecomb_coordinate(frames.state))
        frames['statecomb'].attrs['long_name'] = "State combinations"

    inter_state = frames[iprops]
    for prop in inter_state:
        if 'state' in inter_state[prop].dims:
            inter_state[prop] = subtract_combinations(
                inter_state[prop], dim='state', labels=True
            )
    inter_state = inter_state.map(keep_norming)

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

    if {'energy', 'dip_trans'}.issubset(iprops):
        inter_state = assign_fosc(inter_state)

    return inter_state
