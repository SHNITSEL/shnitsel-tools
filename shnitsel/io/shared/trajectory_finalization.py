import logging
from typing import TypeVar
import xarray as xr

from shnitsel.data.tree import TreeNode

from shnitsel.data.dataset_containers import Trajectory, Frames
from shnitsel.io.shared.helpers import LoadingParameters
from shnitsel.data.state_helpers import (
    default_state_name_assigner,
    default_state_type_assigner,
)
from shnitsel.io.shared.variable_flagging import (
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

            # TODO: FIXME: use loading_parameters to configure state names
            dataset = set_state_defaults(dataset, loading_parameters)

            # TODO: FIXME: Configure cleaning to only run on initial construction trajectories, not when loaded from Shnitsel files.
            unset_vars = []
            for var in dataset.variables:
                if is_variable_assigned(dataset[var]):
                    # Remove tags
                    del dataset[var].attrs["__assigned"]
                else:
                    unset_vars.append(var)

            logging.debug(f"Dropping unset variables: {unset_vars}")
            dataset = dataset.drop_vars(unset_vars)
            convert_all_units_to_shnitsel_defaults(dataset)

            if rebuild_type:
                return rebuild_type(dataset)
            else:
                try:
                    return Trajectory(dataset)
                except:
                    try:
                        return Frames(dataset)
                    except:
                        return dataset
        else:
            # TODO: FIXME: Should we post-process all individual trajectories just in case?
            return dataset

    return None


def set_state_defaults(
    dataset: xr.Dataset, loading_parameters: LoadingParameters | None
) -> xr.Dataset:
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

    logging.debug("Assigning default state names and/or.")

    if not is_variable_assigned(dataset.state_types):
        dataset = default_state_type_assigner(dataset)
    if not is_variable_assigned(dataset.state_names):
        dataset = default_state_name_assigner(dataset)
    return dataset
