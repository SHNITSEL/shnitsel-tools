from itertools import product
import logging
from typing import Any, Iterable, Literal, Sequence, overload

import numpy as np
import xarray as xr

from shnitsel.core._api_info import internal
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.inter_state import InterState
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.state_selection import StateSelection
from shnitsel.units.definitions import energy, dipole

from shnitsel.core.typedefs import DimName, SpectraDictType, StateCombination

from .generic import keep_norming, subtract_combinations
from .._contracts import needs
from ..units import convert_energy, convert_dipole

import ase.units as si


@internal()
def _calculate_fosc(
    energy_interstate: xr.DataArray, dip_trans_norm: xr.DataArray
) -> xr.DataArray:
    """Internal function to actually calculate the oscillator frequency for energies and transition dipoles.

    Parameters
    ----------
    energy_interstate : DataArray
        The Array of Energies in the system.
    dip_trans_norm : DataArray
        The array of associated norm of transition dipoles in the system.


    Returns
    -------
    DataArray
        The resulting oscillation frequency (f_osc) array.

    Notes
    -----
        We use the following unitless form of f_osc:
        $$f_{osc} = \\frac{2}{3} \\frac{m_e}{\\hbar^2} \\cdot \\Delta E \\cdot \\frac{\\mu^2}{e^2}$$
    """

    return (
        2
        / 3
        * si._me
        / (si._hbar**2)
        * (si.Bohr / si.m) ** 2
        * convert_energy(energy_interstate, to=energy.J)
        * convert_dipole(dip_trans_norm, to=dipole.au)
        ** 2  # Until here, we have unit of Bohr^2/m^2 without the  (si.Bohr/si.m)**2 conversion
    )


def calculate_fosc(
    energy_interstate: xr.DataArray, dip_trans_norm: xr.DataArray
) -> xr.DataArray:
    """Function to obtain a dataarray containing the oscillator strength as a dataarray.

    Parameters
    ----------
    energy_interstate : DataArray
        The array of inter-state energies (delta_E) in the system.
    dip_trans_norm : DataArray
        The array of associated transition dipoles in the system with their norm calculated across the direction dimension.

    Returns
    -------
    DataArray
        The resulting datarray of oscillator strength f_osc
    """
    assert energy_interstate.values.shape == dip_trans_norm.values.shape, (
        f"Energy and dip_trans do not have the same shapes: {energy_interstate.values.shape} <-> {dip_trans_norm.values.shape}"
    )

    da = _calculate_fosc(energy_interstate, dip_trans_norm)
    da.name = 'fosc'
    da.attrs.update(
        {
            "long_name": r"$f_{\mathrm{osc}}$",
            "description": "derived from 'energy_interstate' and 'dip_trans' variables",
        }
    )
    return da


@overload
def get_fosc(data: Trajectory | Frames | InterState) -> xr.DataArray: ...


@overload
def get_fosc(
    data: ShnitselDB[Trajectory | Frames | InterState],
) -> ShnitselDB[xr.DataArray]: ...


def get_fosc(
    data: Trajectory
    | Frames
    | InterState
    | TreeNode[Any, Trajectory | Frames | InterState],
) -> xr.DataArray | TreeNode[Any, xr.DataArray]:
    """Function to calculate the strength of the oscillator for state-to-state transitions.

    If provided a simple data type like a Trajectory of Frames, will extract InterState information from it.
    If provided a hierarchical tree structure, will map calculation over data entries.
    Uses energy delta and `dip_trans_norm`, i.e. the norm of the transition dipole for calculation and will only yield results if those are available.
    Otherwise, the function call will fail.

    Parameters
    ----------
    data :  Trajectory
            | Frames
            | InterState
            |TreeNode[Any, Trajectory | Frames | InterState]
        The data from which interstate information can be deduced (transition dipoles and energy deltas).
        Alternatively, the interstate information can be provided directly.
        The data can also be supplied in a hierarchical tree format. In that case, the operation will be applied to data entries
        in the tree structure.

    Returns
    -------
    xr.DataArray | TreeNode[Any, Trajectory | Frames | InterState]
        Either the array of fosc values for simple input structures or the tree structure with the fosc results for individual data entries after mapping
        over the tree.
    """
    if isinstance(data, TreeNode):
        return data.map_data(get_fosc)
    else:
        # Trajectory, Frames or Interstate
        interstate_data: InterState
        if isinstance(data, InterState):
            interstate_data = data
        elif isinstance(data, Frames):
            interstate_data = InterState(data)
        elif isinstance(data, Trajectory):
            interstate_data = InterState(data)
        else:
            raise ValueError(
                "Invalid input type provided for fosc calculation: %s" % type(data)
            )
        return calculate_fosc(
            interstate_data.energy_interstate, interstate_data.dipole_transition_norm
        )


# TODO: deprecate (made redundant by DerivedProperties)
# @needs(data_vars={'energy_interstate', 'dip_trans'}, coords={'statecomb'})
# def assign_fosc(ds: xr.Dataset | Trajectory | Frames) -> xr.Dataset:
#     """Function to calculate oscillator strength fosc and create a new dataset with this variable assigned.

#     Args:
#         ds (xr.Dataset): Dataset from which to calculate fosc

#     Returns:
#         xr.Dataset: Dataset with the member variable fosc set
#     """
#     if 'dip_trans_norm' not in ds:
#         ds['dip_trans_norm'] = keep_norming(ds['dip_trans'])

#     da = _fosc(ds['energy_interstate'], ds['dip_trans_norm'])
#     res = ds.assign(fosc=da)
#     return res


def apply_gauss_broadening(
    delta_E: xr.DataArray,
    fosc: xr.DataArray,
    agg_dim: DimName | None = None,
    *,
    width_in_eV: float = 0.5,  # in eV
    nsamples: int = 1000,
    min_energy_range: float = 0,
    max_energy_range: float | None = None,
) -> xr.DataArray:
    r"""
    Applies a gaussian smoothing kernel to the fosc data.

    Aggregation is performed along the `agg_dim` dimension.

    Parameters
    ----------
    delta_E : xr.DataArray
        Values used for the x-axis, presumably $E_i$
    fosc : xr.DataArray
        Values used for the y-axis, presumably $f_\mathrm{osc}$
    agg_dim : DimName, optional
        Dimension along which to aggregate the many Gaussian distributions,
        by default None, which means, no automatic mean calculation is applied.
    width_in_eV : float, optional
        the width of the Gaussian distributions used, by default 0.5 eV
    nsamples : int, optional
        number of evenly spaced x-values over which to sample the distribution,
        by default 1000
    max_energy_range : float
        The minimum x-value used for the energy-spectrum, by default 0.0
    max_energy_range : float, optional
        the maximum x-value, by default 3 standard deviations
        beyond the pre-broadened maximum
    """

    # TODO: FIXME: Should this remain as-is?
    stdev_in_eV = width_in_eV / 2

    delta_E_eV = convert_energy(delta_E, to=energy.eV)

    def gaussian_filter(x):
        nonlocal stdev_in_eV
        return (
            1
            / (np.sqrt(2 * np.pi) * stdev_in_eV)
            * np.exp(-(x**2) / (2 * stdev_in_eV**2))
        )

    xname = getattr(delta_E_eV, 'name', 'energy_interstate') or 'xdim'
    yname = getattr(fosc, 'name', 'fosc') or 'ydim'

    if max_energy_range is None:
        # TODO: FIXME: The calculation does not fit the statement of the comment above it. Is stdev relative?
        # broadening could visibly overshoot the former maximum by 3 standard deviations
        max_energy_range = delta_E_eV.max().item() + 3 * stdev_in_eV

        assert max_energy_range is not None, (
            "Could not calculate maximum of the provided energy"
        )

    energy_range = np.linspace(min_energy_range, max_energy_range, num=nsamples)
    Espace = xr.DataArray(energy_range, dims=[xname], attrs=delta_E_eV.attrs)
    res: xr.DataArray = gaussian_filter(Espace - delta_E_eV) * fosc

    if agg_dim is not None:
        assert agg_dim in delta_E.sizes, (
            f"E does not have required dimension {agg_dim} for aggregation"
        )
        res = res.mean(dim=agg_dim)

    # print(res)
    res.name = yname
    res.attrs = fosc.attrs
    for cname, coord in res.coords.items():
        if cname in fosc.coords:
            coord.attrs = fosc.coords[cname].attrs
    return res.assign_coords({xname: Espace})


@overload
def get_fosc_gauss_broadened(
    interstate_data: InterState | Trajectory | Frames | xr.Dataset,
    width_in_eV: float = 0.5,
    nsamples: int = 1000,
    max_energy_range: float | None = None,
) -> xr.DataArray: ...


@overload
def get_fosc_gauss_broadened(
    interstate_data: ShnitselDB[InterState | Trajectory | Frames | xr.Dataset],
    width_in_eV: float = 0.5,
    nsamples: int = 1000,
    max_energy_range: float | None = None,
) -> ShnitselDB[xr.DataArray]: ...


@needs(data_vars={'energy_interstate', 'fosc'})
def get_fosc_gauss_broadened(
    interstate_data_source: xr.Dataset
    | InterState
    | Trajectory
    | Frames
    | TreeNode[Any, InterState | Trajectory | Frames | xr.Dataset],
    width_in_eV: float = 0.5,
    nsamples: int = 1000,
    max_energy_range: float | None = None,
) -> xr.DataArray | TreeNode[Any, xr.DataArray]:
    """Function to get the broadened spectrum of the interstate energy and oscillator strength data to plot
    a nice and smooth spectrum.

    Width of the smoothing kernel is given in eV and the energy is assumed to be in eV or will be converted to eV.

    Parameters
    ----------
    interstate_data_source : InterState | Trajectory | Frames | xr.Dataset | TreeNode[Any, InterState  |  Trajectory  |  Frames  |  xr.Dataset]
        Interstate dataset or source for such data with `energy_interstate` and `fosc` information.
            If provided as Frames or Trajectory, must provide `energy` and `dip_trans` data.
            If provided as tree, operation will be mapped over data.
    width_in_eV : float, optional
        Width of the gaussian smoothing kernel in eV, by default 0.5
    nsamples : int, optional
        Number of samples/steps in the range of the energy spectrum, by default 1000
    max_energy_range : float | None, optional
        Maximum of the energy range to consider for the spectrum, by default None

    Returns
    -------
    xr.DataArray | TreeNode[Any, xr.DataArray]
        Resulting broadened spectrum statistics either of the dataset input or mapped over the entire input tree
    """
    if isinstance(interstate_data_source, TreeNode):
        return interstate_data_source.map_data(get_fosc_gauss_broadened)
    else:
        interstate_dataset: xr.Dataset
        if isinstance(interstate_data_source, Trajectory) or isinstance(
            interstate_data_source, Frames
        ):
            interstate_dataset = interstate_data_source.inter_state.dataset
        elif isinstance(interstate_data_source, InterState):
            interstate_dataset = interstate_data_source.dataset
        elif isinstance(interstate_data_source, xr.Dataset):
            interstate_dataset = interstate_data_source
        else:
            raise ValueError(
                "Unsupported type for fosc generation and gauss broadening: %s"
                % type(interstate_data_source)
            )

        return apply_gauss_broadening(
            interstate_dataset.delta_energy,
            interstate_dataset.fosc,
            width_in_eV=width_in_eV,
            nsamples=nsamples,
            max_energy_range=max_energy_range,
        )


@overload
def get_spectrum_at_time(
    interstate_data: ShnitselDB[InterState | xr.Dataset],
    t: float,
    sc: StateCombination,
    rel_cutoff: float = 0.01,
) -> ShnitselDB[xr.Dataset] | None: ...
@overload
def get_spectrum_at_time(
    interstate_data: InterState | xr.Dataset | Sequence[InterState | xr.Dataset],
    t: float,
    sc: StateCombination,
    rel_cutoff: float = 0.01,
) -> xr.DataArray | None: ...


@needs(data_vars={'energy_interstate', 'fosc'}, coords={"statecomb", "time"})
def get_spectrum_at_time(
    interstate_data: InterState
    | xr.Dataset
    | Sequence[InterState | xr.Dataset]
    | TreeNode[Any, InterState | xr.Dataset],
    t: float,
    sc: StateCombination,
    rel_cutoff: float = 0.01,
) -> xr.DataArray | TreeNode[Any, xr.DataArray] | None:
    """Function to calculate a gaussian-smoothed spectrum of an interstate dataset at one specific point in time and for one specific state transition

    _extended_summary_

    Parameters
    ----------
    interstate_data : InterState | xr.Dataset | Sequence[InterState | xr.Dataset ] | TreeNode[Any, InterState | xr.Dataset]
        An InterState dataset with fosc and energy data.
        Alternatively, a sequence of Interstate data or a hierarchical set of interstate data.
    t : float
        The time at which to evaluate the spectrum
    sc : StateCombination
        State combination identifier. Provided as a tuple ``(from, to)`` of state indices.
    rel_cutoff : float, optional
        Relative cutoff threshold. Values below the max of the resulting spectrum times this scale will be ignored, by default 0.01

    Returns
    -------
    xr.DataArray
        The Gauss-broadened spectrum of the provided `data` system.
        If broadening across trajectories could not be performed, just returns the fosc array.
        If a sequence of interstate data was provided, aggregation will be performed across
        the spectra of the different datasets.
    TreeNode[Any, xr.DataArray]
        Hierarchically mapped data, where spectrum calculation was performed across flat groups.
        If a hierarchical tree set was provided, the data will be grouped and aggregation performed across the spectra of data within the same group
    None
        If the spectrum could not be calculated for whatever reason. Most likely due to missing data for the spectrum calculation.

    Raises
    ------
    ValueError
        Unsupported type provided.
    """
    if isinstance(interstate_data, TreeNode):
        interstate_data_grouped = interstate_data.group_data_by_metadata()
        # mapped_spectrum_data = interstate_data_grouped.map_data(
        #     lambda x: get_spectrum_at_time(x, t, sc, 0.0)
        # )

        spectrum_tree = interstate_data_grouped.map_flat_group_data(
            lambda x: get_spectrum_at_time(list(x), t, sc, rel_cutoff=rel_cutoff)
        )
        return spectrum_tree

    elif isinstance(interstate_data, Sequence):
        spectra_results: list[xr.DataArray] = []
        for traj_data in interstate_data:
            if isinstance(traj_data, InterState):
                traj_data = traj_data.dataset

            if t > traj_data.coords['time'].max():
                # Don't produce a spectrum after the end of the trajectory
                continue

            if t not in traj_data.coords['time']:
                times = np.unique(traj_data.coords['time'])
                diffs = np.abs(times - t)
                curr_t = times[np.argmin(diffs)]
            else:
                curr_t = t

            tmp_res = get_spectrum_at_time(traj_data, curr_t, sc, 0.00)
            if tmp_res is not None:
                spectra_results.append(tmp_res)

        num_spectra = len(spectra_results)
        if num_spectra == 0:
            return None
        else:
            gauss_broadened_point_spectrum = spectra_results[0]
            for i in range(1, num_spectra):
                gauss_broadened_point_spectrum += spectra_results[i]
            gauss_broadened_point_spectrum /= float(num_spectra)
            gauss_broadened_point_spectrum.name = spectra_results[0].name
            gauss_broadened_point_spectrum.attrs.update(spectra_results[0].attrs)

            # Perform cutoff filtering:
            max_ = gauss_broadened_point_spectrum.max().item()
            non_negligible = gauss_broadened_point_spectrum.where(
                gauss_broadened_point_spectrum > rel_cutoff * max_, drop=True
            ).energy_interstate
            if len(non_negligible) == 0:
                return gauss_broadened_point_spectrum.sel(
                    energy_interstate=non_negligible
                )
            return gauss_broadened_point_spectrum.sel(
                energy_interstate=slice(non_negligible.min(), non_negligible.max())
            )
    else:
        interstate_ds: xr.Dataset
        if isinstance(interstate_data, InterState):
            interstate_ds = interstate_data.dataset
        elif isinstance(interstate_data, xr.Dataset):
            interstate_ds = interstate_data
        else:
            raise ValueError(
                "Unsupported type provided to spectrum calculation: %s"
                % type(interstate_data)
            )

        # TODO: FIXME: Deal with hierarchical data being processed, i.e. take a set of Interstate data and do manual gauss broadening
        # following required because `method='nearest'` doesn't work for MultiIndex
        if t not in interstate_ds.coords['time']:
            times = np.unique(interstate_ds.coords['time'])
            diffs = np.abs(times - t)
            t = times[np.argmin(diffs)]

        # Only take one timestep and one state combination
        interstate_ds = interstate_ds.sel(time=t, statecomb=sc)

        # Figure out how the trajectory is indexed across multiple trajectories or whether a single trajectory is provided
        trajid_dim = None
        if "active_trajectory" in interstate_ds.energy_interstate.sizes:
            trajid_dim = "active_trajectory"
        elif "trajectory" in interstate_ds.energy_interstate.sizes:
            trajid_dim = "trajectory"

        if trajid_dim is None:
            logging.info("Single Trajectory provided, no aggregation performed")

        gauss_broadened_point_spectrum: xr.DataArray = apply_gauss_broadening(
            interstate_ds.energy_interstate, interstate_ds.fosc, agg_dim=trajid_dim
        )
        max_ = gauss_broadened_point_spectrum.max().item()
        non_negligible = gauss_broadened_point_spectrum.where(
            gauss_broadened_point_spectrum > rel_cutoff * max_, drop=True
        ).energy_interstate
        if len(non_negligible) == 0:
            return gauss_broadened_point_spectrum.sel(energy_interstate=non_negligible)
        return gauss_broadened_point_spectrum.sel(
            energy_interstate=slice(non_negligible.min(), non_negligible.max())
        )


@overload
def get_spectra(
    interstate_data: InterState | xr.Dataset | Sequence[InterState | xr.Dataset],
    state_selection: StateSelection | None = None,
    times: Iterable[float] | Literal['all'] | None = None,
    rel_cutoff: float = 0.01,
) -> SpectraDictType: ...


@overload
def get_spectra(
    interstate_data: ShnitselDB[InterState | xr.Dataset],
    state_selection: StateSelection | None = None,
    times: Iterable[float] | Literal['all'] | None = None,
    rel_cutoff: float = 0.01,
) -> ShnitselDB[SpectraDictType]: ...


@needs(data_vars={'energy', 'fosc'}, coords={"statecomb", "time"})
def get_spectra(
    interstate_data: InterState
    | xr.Dataset
    | Sequence[InterState | xr.Dataset]
    | TreeNode[Any, InterState | xr.Dataset],
    state_selection: StateSelection | None = None,
    times: Iterable[float] | Literal['all'] | None = None,
    rel_cutoff: float = 0.01,
) -> SpectraDictType | TreeNode[Any, SpectraDictType]:
    """Function to calculate (gauss-broadened) spectra at multiple (or all) points in time

    Uses strength of oscillator (`fosc`) and energy deltas (`energy_interstate`) to calculate
    a smoothened spectrum across the energy phase space.
    The times at which the spectrum is calculated and the list of state transitions to consider
    can be controlled via the parameters.

    Parameters
    ----------
    interstate : InterState | xr.Dataset | Sequence[InterState | xr.Dataset] | TreeNode[Any, InterState | xr.Dataset]
        The data source for interstate data that needs to provide interstate energy differenses and fosc data,
        which can be derived from `energy_interstate` and `dip_trans`.
        If the necessary data is not provided, an error will be raised.
        If a sequence is provided, state selection defaults will be generated from the first data set in the selection.
        If a hierarchical tree is provided, data will be grouped and the spectra calculation will be mapped over
        the data within each group.
    state_selection : StateSelection, optional
        State combination selection provided either as a `StateSelection` instance or, by default will consider all state combinations.
    times : Iterable[float] | Literal['all'] | None, optional
        Specific times at which the spectrum should be calculated or `all` if the spectrum should be extracted at all times, by default None, which means that
        a set of times will be chosen automatically. (Note: will currently be initialized as [0,10,20,30] in arbitrary time units)
    rel_cutoff : float, optional
        Factor for the cutoff of broadened/smoothened spectrum relative to maximum to be considered, by default 0.01

    Raises
    ------
    ValueError
        Unsupported type provided.

    Returns
    -------
    SpectraDictType
        The resulting spectrum as a mapping from `(t,sc)` pairs to the resulting, broadened spectrum.
        The first item in the key is the time at which the spectrum was extracted.
        The second item is the state combination descriptor for which it was calculated.
    TreeNode[Any, SpectraDictType]
        If a hierarchical input is provided, the spectra results will be provided in hierarchical structure with one result for every grouped set of data.
    """
    if times is None:
        # TODO: FIXME: Make choice of times more sophisticated
        times = [0, 10, 20, 30]

    if isinstance(interstate_data, TreeNode):
        interstate_data_grouped = interstate_data.group_data_by_metadata()

        spectrum_tree = interstate_data_grouped.map_flat_group_data(
            lambda x: get_spectra(
                list(x),
                state_selection=state_selection,
                times=times,
                rel_cutoff=rel_cutoff,
            )
        )
        return spectrum_tree
    elif isinstance(interstate_data, Sequence):
        if times == 'all':
            times = set()
            for data in interstate_data:
                times.update(data.coords['time'])
            times = list(times)
            times.sort()

        if state_selection is None:
            # Use all combinations if no selection provided
            if isinstance(interstate_data[0], InterState):
                state_selection = StateSelection.init_from_dataset(
                    interstate_data[0].dataset
                )
            else:
                state_selection = StateSelection.init_from_dataset(interstate_data[0])

        sc_values: Iterable[tuple[int, int]] = state_selection.state_combinations

        # TODO: FIXME: We should consider making this a data array with nice dimensions instead
        res: SpectraDictType = {
            (t, sc): res_spec
            for t, sc in product(times, sc_values)
            if (
                res_spec := get_spectrum_at_time(
                    interstate_data, t=t, sc=sc, rel_cutoff=rel_cutoff
                )
            )
            is not None
        }
        return res
    else:
        # TODO: FIXME: Allow tree as input and apply to collection in grouped data manner
        interstate_dataset: xr.Dataset

        if isinstance(interstate_data, InterState):
            interstate_dataset = interstate_data.dataset
        elif isinstance(interstate_data, xr.Dataset):
            interstate_dataset = interstate_data
        else:
            raise ValueError(
                "Unsupported type provided for interstate data: %s"
                % type(interstate_data)
            )

        if times == 'all':
            times = interstate_dataset.coords['time'].values

        if state_selection is None:
            # Use all combinations if no selection provided
            state_selection = StateSelection.init_from_dataset(interstate_dataset)

        sc_values: Iterable[tuple[int, int]] = state_selection.state_combinations

        # TODO: FIXME: We should consider making this a data array with nice dimensions instead
        res: SpectraDictType = {
            (t, sc): res_spec
            for t, sc in product(times, sc_values)
            if (
                res_spec := get_spectrum_at_time(
                    interstate_data, t=t, sc=sc, rel_cutoff=rel_cutoff
                )
            )
            is not None
        }
        return res


def get_spectra_groups(
    spectra: SpectraDictType,
    state_selection: StateSelection | None = None,
) -> tuple[SpectraDictType, SpectraDictType]:
    """Group spectra results into spectra involving the ground state or only excited states.

    _extended_summary_

    Parameters
    ----------
    spectra : SpectraDictType
        The Spectral calculation results, e.g. from `calc_spectra()`. Indexed by (timestep, state_combination) and yielding the associated spectrum.
    state_selection : StateSelection | None, optional
        The selection of states to consider as ground and active states, by default None.
        If not provided, all state transitions with state ids above 1 will be considered `excited` all others `ground` state transitions.
        If provided, the excited state combinations will be extracted using `state_selection.excited_state_transitions()` with default parameters.

    Returns
    -------
    tuple[ SpectraDictType, SpectraDictType]
        First the spectra involving the ground state
        Second the spectra involving only excited states.
    """
    ground, excited = {}, {}

    if state_selection is None:
        for (t, (sc_from, sc_to)), v in spectra.items():
            if sc_from > 1 and sc_to > 1:
                excited[t, (sc_from, sc_to)] = v
            else:
                ground[t, (sc_from, sc_to)] = v
    else:
        excited_transitions = (
            state_selection.excited_state_transitions().state_combinations
        )
        for (t, sc), v in spectra.items():
            if sc in excited_transitions:
                excited[t, sc] = v
            else:
                ground[t, sc] = v

    sgroups = (ground, excited)
    return sgroups


# TODO: FIXME: This looks like it should be covered by get_spectra?
# @needs(data_vars={'energy', 'fosc'}, coords={'frame', 'trajid_'})
# def spectra_all_times(inter_state: InterState | xr.Dataset) -> xr.DataArray:
#     """Function to calculate the spectra at all times.

#     Does not return a dict with only the relevant (t,sc) combinations as above but instead a full
#     xr.DataArray with a time dimension that has spectrum data for all times within the dataset averaged across trajectories.

#     Args:
#         inter_state (xr.Dataset): The InterState transformed Dataset.

#     Raises:
#         ValueError: If required variables or dimensions are missing

#     Returns:
#         xr.DataArray: The resulting spectra across all times.
#     """
#     assert isinstance(inter_state, xr.Dataset)
#     if 'energy' not in inter_state.data_vars:
#         raise ValueError("Missing required variable 'energy'")
#     if 'fosc' not in inter_state.data_vars:
#         raise ValueError("Missing required variable 'fosc'")
#     assert 'frame' in inter_state and 'active_trajectory' in inter_state, (
#         "Missing required dimensions"
#     )
#     # TODO: FIXME: This probably should not have to unstack here? We should just accept a tree and use each trajectory individually and then aggregate over trajectories?
#     data = inter_state.unstack('frame')
#     return apply_gauss_broadening(data.energy, data.fosc, agg_dim='active_trajectory')
