import logging
from typing import TypeVar
import xarray as xr

from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree import TreeNode

from shnitsel.data.dataset_containers import Trajectory, Frames, wrap_dataset
from shnitsel.io.shared.helpers import LoadingParameters
from shnitsel.data.state_helpers import (
    default_state_name_assigner,
    default_state_type_assigner,
)
from shnitsel.io.shared.variable_flagging import (
    clean_unassigned_variables,
    is_variable_assigned,
)
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults

NodeType = TypeVar("NodeType", bound=TreeNode)
DataType = TypeVar("DataType", bound=xr.Dataset | Trajectory | Frames)


def finalize_loaded_trajectory(
    dataset: DataType | NodeType | None,
    loading_parameters: LoadingParameters | None,
) -> DataType | Trajectory | Frames | NodeType | None:
    """Function to apply some final postprocessing common to all input routines that allow reading of single trajectories from input formats.

    Parameters
    ----------
    dataset : xr.Dataset | Trajectory | Frames | NodeType | None
        The dataset to perform finalization on. Only updates Dataset, Trajectory and Frames data.
        All other data will be returned unchanged
    loading_parameters : LoadingParameters | None
        Parameters to set some defaults.

    Returns
    -------
    xr.Dataset | Trajectory | Frames | NodeType | None
        The same type as the original `dataset` parameter but potentially with some
        default values and conversions applied.
    """
    if dataset is not None:
        # logging.debug(f"Finalizing: {repr(dataset)}")
        if isinstance(dataset, (xr.Dataset, Trajectory, Frames)):
            rebuild_type = None
            if isinstance(dataset, (Trajectory, Frames)):
                rebuild_type = type(dataset)
                dataset = dataset.dataset

            dataset = set_state_defaults(dataset, loading_parameters)
            # Clean up variables if the variables are not assigned yet.
            dataset = clean_unassigned_variables(dataset)
            dataset = convert_all_units_to_shnitsel_defaults(dataset)

            if rebuild_type:
                return rebuild_type(dataset)
            else:
                return wrap_dataset(dataset, (Trajectory | Frames))
        else:
            # Should we post-process all individual trajectories just in case?
            if isinstance(dataset, TreeNode):
                return dataset.map_data(
                    lambda x: finalize_loaded_trajectory(
                        x, loading_parameters=loading_parameters
                    ),
                    keep_empty_branches=True,
                )
            return dataset

    return None


def set_state_defaults(
    dataset: xr.Dataset, loading_parameters: LoadingParameters | None
) -> xr.Dataset:
    """Helper function to apply default settings to dataset variables
    for state names and state types if they have not been assigned at some point
    earlier during the configuration process.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset to set state name and state type information on.
    loading_parameters : LoadingParameters | None
        Currently unused settings to be applied to trajectory import.

    Returns
    -------
    xr.Dataset
        The dataset but with default values for state types and state names
    """
    # TODO: FIXME: apply configured names from loading_parameters

    if is_variable_assigned(dataset.state_types) and is_variable_assigned(
        dataset.state_names
    ):
        logging.debug(
            "Types and names of state already set for dataset in finalization."
        )

        # logging.debug(f"Types: {dataset.state_types}")
        # logging.debug(f"Names: {dataset.state_names}")
        return dataset

    logging.debug("Assigning default state names and/or types.")

    if not is_variable_assigned(dataset.state_types):
        dataset = default_state_type_assigner(dataset)
    if not is_variable_assigned(dataset.state_names):
        dataset = default_state_name_assigner(dataset)
    return dataset


def normalize_dataset(ds: xr.Dataset) -> xr.Dataset:
    """Helper method to perform some standardized renaming operations as well as some
    restructuring, e.g. if a multi-index is missing.

    May also convert some legacy attributes to promoted dimensionless variables.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to normalize according to current Shnitsel standards

    Returns
    -------
    xr.Dataset
        The renamed dataset.
    """
    # Rename legacy uses of `trajid` and `trajid_`
    if 'trajid' in ds.dims:
        ds = ds.swap_dims(trajid='trajectory')
    if 'trajid_' in ds.dims:
        ds = ds.swap_dims(trajid='trajectory')

    if 'trajid' in ds.coords:
        if 'trajectory' not in ds.coords['trajid'].dims:
            ds = ds.rename(trajid="atrajectory")
        else:
            ds = ds.rename(trajid="trajectory")

    if 'trajid_' in ds.coords:
        if 'trajectory' not in ds.coords['trajid_'].dims:
            ds = ds.rename(trajid_="atrajectory")
        else:
            ds = ds.rename(trajid_="trajectory")

    # Check if frameset has a multi-index
    if 'frame' in ds.dims and 'time' in ds.coords and 'frame' not in ds.xindexes:
        if 'frame' in ds.coords['time'].dims:
            # Make frame a multi-index
            ds = ds.set_xindex('frame')

    # Turn delta_t and t_max and max_ts into a variable
    if 'time' in ds.coords:
        time_unit = ds.time.attrs.get('units', 'fs')
    else:
        time_unit = 'fs'

    for var, has_unit in [('delta_t', True), ('max_ts', False), ('t_max', True)]:
        if var not in ds.coords:
            if var in ds.data_vars:
                ds = ds.set_coords(var)
            elif var in ds.attrs:
                dt_arr = xr.DataArray(float(ds.attrs.get(var, -1)), dims=(), name=var)
                ds = ds.assign_coords(delta_t=dt_arr)

            if has_unit and var in ds and 'units' not in ds.delta_t.attrs:
                ds.delta_t.attrs['unitdim'] = 'time'
                ds.delta_t.attrs['units'] = time_unit

    return ds
