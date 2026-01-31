from functools import cached_property
import logging
from typing import Any, Iterable, overload

import numpy as np
import xarray as xr

from shnitsel.core._api_info import API, internal
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.units.definitions import time, unit_dimensions

from .._contracts import needs


class PopulationStatistics:
    """Class to hold population statistics of a homogeneous/compatible set
    of trajectories.

    Holds absolute and relative statistics, with the latter being derived from the prior.
    """

    _absolute_statistics: xr.DataArray

    def __init__(self, base_traj_or_data: Frames | xr.DataArray):
        if isinstance(base_traj_or_data, xr.DataArray):
            self._absolute_statistics = base_traj_or_data
        else:
            self._absolute_statistics = (
                PopulationStatistics._calc_classical_populations(base_traj_or_data)
            )

    @property
    def absolute(self) -> xr.DataArray:
        """Returns the absolute population statistics across `time`.

        Values should be integer and of the shape (time, state)
        where the value in [t,s] is the number of trajectories with active state `s`
        at time `t`

        Returns
        -------
        xr.DataArray
            The array with absolute population statistics over time
        """
        return self._absolute_statistics

    @cached_property
    def relative(self) -> xr.DataArray:
        """Returns the relative population statistics across `time`, i.e. the ratio of trajectories in every state.

        Values should be float in the range `[0,1]` and of the shape (time, state)
        where the value at position [t,s] is the ratio of trajectories still running up to this point with active state `s`
        at time `t`.

        Note for interpretation that only having one trajectory left at time `t` with state `s` would make the value at `[t,s]` a full 1.

        Returns
        -------
        xr.DataArray
            The array with relative population statistics over time
        """
        abs_data = self.absolute
        res = (
            (abs_data / abs_data.sum('state'))
            .assign_coords(abs_data.coords)
            .assign_attrs(
                {'long_name': "Relative population of different states over times"}
            )
        )
        res.name = "rel_populations"
        return res

    @internal()
    @staticmethod
    def _calc_classical_populations(
        frames: Frames,
    ) -> xr.DataArray:
        """Function to calculate classical state populations from the active state information in `astate` of the dataset `frames.

        Does not use the partial QM coefficients of the states.

        Parameters
        ----------
        frames : Frames
            The dataset holding the active state information in a variable `astate`.

        Returns
        -------
        xr.DataArray
            The array holding the absolute number of trajectories in each respective state.
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
        ).assign_attrs(
            {'long_name': "Absolute population of different states over times"}
        )


@overload
def calc_classical_populations(
    data: ShnitselDB[Trajectory | Frames],
) -> ShnitselDB[PopulationStatistics]:
    """Specialized version of the population calculation where a tree hierarchy of Trajectory data
    is mapped to a tree hierarchy of population statistics.

    The hierarchy will first be grouped by metadata, then a population statistics `xr.DataArray` will
    be calculated for each flat group.>

    Parameters
    ----------
    data : ShnitselDB[Trajectory | Frames]
        The tree-structured trajectory data

    Returns
    -------
    ShnitselDB[PopulationStatistics]
        The tree structure holding population statistics for grouped data out of the input.
        Results contain inidividual states' absolute population numbers for every time step.
    """
    ...


@overload
def calc_classical_populations(
    data: Trajectory | Frames | xr.Dataset,
) -> PopulationStatistics:
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
    frames: TreeNode[Any, Trajectory | Frames | MultiSeriesDataset | xr.Dataset]
    | Trajectory
    | Frames
    | MultiSeriesDataset
    | xr.Dataset,
) -> PopulationStatistics | TreeNode[Any, PopulationStatistics]:
    """Function to calculate classical state populations from the active state information in `astate` of the dataset `frames.

    Does not use the partial QM coefficients of the states.

    Parameters
    ----------
    frames : Frames
        The dataset holding the active state information in a variable `astate`.

    Returns
    -------
    PopulationStatistics
        The object holding population statistics (absolute+relative) in each respective state.
    ShnitselDB[PopulationStatistics]
        The tree holding the hierarchical population statistics for each flat group in the tree with compatible metadata.
    """
    if isinstance(frames, TreeNode):
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
        db: TreeNode[Any, MultiSeriesDataset | Trajectory | Frames | xr.Dataset] = (
            frames.group_data_by_metadata()
        )

        def _map_prepare(
            frames: MultiSeriesDataset | Trajectory | Frames | xr.Dataset,
        ) -> PopulationStatistics:
            wrapped_ds = wrap_dataset(frames, MultiSeriesDataset | Trajectory | Frames)
            pop_frames: Frames
            if isinstance(wrapped_ds, Trajectory) and not isinstance(
                wrapped_ds, Frames
            ):
                pop_frames = frames.as_frames
            elif isinstance(wrapped_ds, MultiSeriesDataset):
                pop_frames = wrapped_ds.as_stacked
            else:
                pop_frames = wrapped_ds

            return PopulationStatistics(pop_frames)

        mapped_pop_data = db.map_data(_map_prepare, keep_empty_branches=False)

        def _combine_func(
            pop_data: Iterable[PopulationStatistics],
        ) -> PopulationStatistics:
            res_abs: xr.DataArray | None = None
            res_timelen: int = -1
            for data in pop_data:
                data_abs = data.absolute
                if data is None:
                    continue

                if res_abs is None:
                    # Copy into result for the first match
                    res_abs = data_abs
                    res_timelen = res_abs.sizes['time']
                else:
                    new_timelen = data_abs.sizes['time']

                    time_values: xr.DataArray | None = None

                    if new_timelen < res_timelen:
                        data_abs = data_abs.pad(
                            {'time': (0, res_timelen - new_timelen)},
                            constant_values=0.0,
                            keep_attrs=True,
                        )
                        time_values = res_abs.time
                    else:
                        res_abs = res_abs.pad(
                            {'time': (0, new_timelen - res_timelen)},
                            constant_values=0.0,
                            keep_attrs=True,
                        )
                        time_values = data_abs.time()

                    res_abs = res_abs + data_abs
                    if time_values is not None:
                        res_abs = res_abs.assign_coords(
                            {'time': ('time', time_values, time_values.attrs)}
                        )
            if res_abs is None:
                # No population data, empty array
                return PopulationStatistics(xr.DataArray())
            # Yield the absolute population for this group
            return PopulationStatistics(res_abs)

        reduced_pop_data = mapped_pop_data.map_flat_group_data(_combine_func)
        return reduced_pop_data
    else:
        eventual_frames: Frames

        wrapped_dataset = wrap_dataset(frames, Trajectory | Frames | MultiSeriesDataset)
        if isinstance(wrapped_dataset, MultiSeriesDataset):
            eventual_frames = wrapped_dataset.as_stacked
        else:
            eventual_frames = wrapped_dataset.as_frames

        population_results = PopulationStatistics(eventual_frames)
        return population_results


# Alternative name of the function to calculate population statistics
calc_pops = calc_classical_populations
