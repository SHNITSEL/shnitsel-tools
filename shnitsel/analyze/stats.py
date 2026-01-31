from logging import warning
import logging
from typing import Hashable

import numpy as np
import numpy.typing as npt
import scipy.stats as st
import xarray as xr

from shnitsel.core._api_info import internal
from shnitsel.io.shared.trajectory_setup import get_statecomb_coordinate
from shnitsel.units.definitions import energy
from shnitsel.units.conversion import convert_energy

from .generic import keep_norming, subtract_combinations as subtract_combinations
from .spectra import calculate_fosc

from .._contracts import needs

from ..core.typedefs import DimName

from ..data.dataset_containers import (
    Frames,
    Trajectory,
    MultiSeriesDataset,
    wrap_dataset,
)


#####################################################
# For calculating confidence intervals, the following
# functions offer varying levels of abstraction
# TODO make naming consistent


def calc_confidence_interval(
    data_array: npt.NDArray, confidence: float = 0.95
) -> npt.NDArray:
    """Function to calculate the confidence interval for a variable array `a`.

    The result is a numpy array with stacked entries with the lower and upper limits of the confidence interval.

    Parameters
    ----------
    a : npt.NDArray
        The Numpy array to calculate the confidence interval for.
    confidence : float, optional
        The confidence level to get the confidence interval for. Defaults to 0.95.

    Raises
    ------
    ValueError
        Raised if the provided `data_array` is not one-dimensional

    Returns
    -------
    npt.NDArray
        Numpy array with lower and upper bounds of the confidence interval
    """
    if np.array(data_array).ndim != 1:
        raise ValueError("This function accepts 1D input only")

    data_mean = np.mean(data_array)
    if len(data_array) < 3:
        lower = np.min(data_array)
        upper = np.max(data_array)
        return np.array([lower, upper])
    else:
        std_error = st.sem(data_array)
        res_interval = st.t.interval(
            confidence,
            len(data_array) - 1,
            loc=data_mean,
            scale=std_error,
        )
        # print(res_interval)
        # print(np.stack(res_interval))
        # raise ValueError("stop")
        return np.stack(res_interval)

    return np.stack(
        st.t.interval(
            confidence,
            len(data_array) - 1,
            loc=np.mean(data_array),
            scale=st.sem(data_array),
        )
    )


def confidence_interval_aggregate_last_dim(
    data_array: npt.NDArray, confidence=0.95
) -> npt.NDArray:
    """Calculate the confidence interval from statistics aggregated across the last dimension.

    For our purposes, this should amount to the trajectory being averaged over.

    Parameters
    ----------
    data_array : npt.NDArray
        The numpy data array to calculate the confidence interval for.
    confidence : float, optional
        The confidence level to use for calculations. Defaults to 0.95.

    Returns
    -------
    npt.NDArray
        A numpy array with (lower_bound,upper_bound,mean) of the confidence interval in the last dimension. Otherwise same shape as data_array.
    """
    outer_shape = tuple(data_array.shape[:-1])
    res = np.full(outer_shape + (3,), np.nan)
    for idxs in np.ndindex(outer_shape):
        res[idxs, :2] = calc_confidence_interval(
            data_array[idxs], confidence=confidence
        )
        res[idxs, 2] = np.mean(data_array[idxs])
    return res


@internal()
def calc_confidence_interval_in_array_dimensions(
    data_array: xr.DataArray, dim: DimName, confidence: float = 0.95
) -> xr.DataArray:
    """Function to calculate confidence interval data for the input data_array.
    Results are then repackaged back into an xr.DataArray, where the dimension `bound` allows to choose between confidence interval limits and the mean of the distribution.

    The dimension denoted by `dim` will be aggregated across.

    Parameters
    ----------
    data_array : xr.DataArray
        Input data to have confidence intervals calculated for.
    dim : DimName
        Dimension to calculate the confidence interval data from.
    confidence : float, optional
        Confidence level for Confidence interval calculation. Defaults to 0.95.

    Returns
    -------
    xr.DataArray
        DataArray with coordinate `bound` with values 'lower', 'upper', and 'mean', which refer to the lower and the upper bound of the confidence interval of this and
    """
    res_da: xr.DataArray = xr.apply_ufunc(
        confidence_interval_aggregate_last_dim,
        data_array,
        kwargs={"confidence": confidence},
        output_core_dims=[["bound"]],
        input_core_dims=[[dim]],
    )
    return res_da.assign_coords(  #
        dict(bound=["lower", "upper", "mean"])
    )


@needs(groupable={"time"}, dims={"frame"})
def time_grouped_confidence_interval(
    data_array: xr.DataArray, confidence: float = 0.9
) -> xr.Dataset:
    """Function to calculate the per-time confidence interval of a DataArray that is groupable by the `time` coordinate.

    Parameters
    ----------
    data_array : xr.DataArray
        Data Array for whose data the confidence intervals should be calculated
    confidence : float, optional
        The confidence level to calculate the interval bounds for. Defaults to 0.9.

    Returns
    -------
    xr.Dataset
        A new Dataset, where variables 'lower', 'upper' and 'mean' contain the lower and upper bounds of the confidence interval in each time step and mean is the mean at each point in time.
    """
    if "frame" in data_array.dims:
        return (
            data_array.groupby("time")
            .map(
                lambda x: calc_confidence_interval_in_array_dimensions(
                    x, dim="frame", confidence=confidence
                )
            )
            .to_dataset("bound")
        )
    elif "time" in data_array.dims:
        return (
            data_array.groupby("time")
            .map(lambda x: xr.DataArray(np.array([x, x, x]), dims=["bound"]))
            .assign_coords(dict(bound=["lower", "upper", "mean"]))
            .to_dataset("bound")
        )
    else:
        raise ValueError("Data contained neither time nor frame dimension.")


@needs(dims={"state"})
def get_per_state(
    frames: Frames | Trajectory | MultiSeriesDataset | xr.Dataset,
) -> xr.Dataset:
    """Isolate the standard per-state properties (energy, forces, permanent dipoles)
    from an xr.Dataset, and take their norm over all array dimensions other than 'state'
    so that the resulting variables can be easily plotted against another.

    Parameters
    ----------
    frames : Frames | Trajectory | MultiSeriesDataset | xr.Dataset
        An xr.Dataset object containing at least 'energy', 'forces' and 'dip_perm' variables

    Returns
    -------
    xr.Dataset
        An xr.Dataset object containing only 'energy', 'forces' and 'dip_perm' variables
    """
    # TODO: FIXME: Attributes need to be kept.
    # And why create a new dataset instead of amending the original one?
    props_per = {"energy", "forces", "dip_perm"}.intersection(frames.keys())
    per_base_props_to_norm = {"dip_perm", "forces"}.intersection(frames.keys())

    per_state = frames[props_per]

    for prop in per_base_props_to_norm:
        if prop in frames:
            per_state[str(prop) + "_norm"] = keep_norming(frames[prop])

    if "forces" in per_state:
        per_state["forces"] = per_state["forces"].where(per_state["forces"] != 0)
        per_state["forces"].attrs["long_name"] = r"$\mathbf{F}$"

    if "energy" in per_state:
        per_state["energy"].attrs["long_name"] = r"$E$"

    if "dip_perm" in per_state:
        per_state["dip_perm"].attrs["long_name"] = r"$\mathbf{\mu}_i$"

    per_state.attrs.update(frames.attrs)
    return per_state


@needs(dims={"state"}, coords={"state"})
def get_inter_state(
    frames: Frames | Trajectory | MultiSeriesDataset | xr.Dataset,
) -> xr.Dataset:
    """Calculate inter-state properties of a dataset for certain observables.

    Currently calculates inter-state levels of energy differences.
    Will calculate Differences between the values of these observables indexed by state.
    If no `statecomb` dimension exists, will create one.

    Parameters
    ----------
    frames : Frames | Trajectory | MultiSeriesDataset | xr.Dataset
        The basis Dataset to calculate the interstate properties for

    Returns
    -------
    xr.Dataset
        A Dataset containing interstate properties
    """
    prop: Hashable
    inter_base_props = ["energy"]
    available_inter_base_props = []
    inter_base_props_to_norm = ["dip_trans", "nacs", "socs"]
    # TODO: Check that energy is actually the only inter-state property. We already have statecomb for nacs and dip_trans and astate is not state-dependent.
    for prop in inter_base_props:
        if prop in frames:
            available_inter_base_props.append(prop)
        else:
            warning(f"Dataset does not contain variable '{prop}'")

    if "statecomb" not in frames.sizes:
        logging.info(
            "Creating a new `statecomb` dimension because it was not yet set when calculating inter-state properties."
        )
        statecomb_coords = get_statecomb_coordinate(frames.state)
        frames = frames.assign_coords(statecomb_coords)
        frames["statecomb"].attrs["long_name"] = "State combinations"

    inter_state = frames.copy()
    for prop in available_inter_base_props:
        if "state" in frames[prop].dims and "statecomb" not in frames[prop].dims:
            inter_state_res = subtract_combinations(
                frames[prop], dim="state", add_labels=False
            )
            inter_state[str(prop) + "_interstate"] = inter_state_res
            inter_state[str(prop) + "_interstate"].attrs["long_name"] = (
                f"Derived inter-state differences of variable `{prop}`"
            )

    for prop in inter_base_props_to_norm:
        if prop in frames:
            inter_state[str(prop) + "_norm"] = keep_norming(frames[prop])

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

    if {"energy_interstate", "dip_trans"}.issubset(inter_state.variables.keys()):
        fosc_data = calculate_fosc(
            inter_state.energy_interstate, inter_state.dip_trans_norm
        )
        inter_state = inter_state.assign(fosc=fosc_data)

    if "energy_interstate" in inter_state:
        inter_state["energy_interstate"] = convert_energy(
            inter_state["energy_interstate"], to=energy.eV
        )
        inter_state["energy_interstate"].attrs["long_name"] = (
            "Energy delta between the energy levels of various states derived from `energy`"
        )

    # TODO: FIXME: Consider whether the result should contain only statecomb
    # inter_state = inter_state.drop_dims(['state'], errors='ignore')
    return inter_state
