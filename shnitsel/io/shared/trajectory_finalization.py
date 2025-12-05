import logging
import xarray as xr

from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.data.trajectory_format import Trajectory
from shnitsel.data.trajectory_format_wrapper import wrap_trajectory
from shnitsel.io.shared.helpers import LoadingParameters
from shnitsel.data.state_helpers import (
    default_state_name_assigner,
    default_state_type_assigner,
)
from shnitsel.io.shared.variable_flagging import (
    is_variable_assigned,
)
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults


def finalize_loaded_trajectory(
    dataset: xr.Dataset | ShnitselDB | None,
    loading_parameters: LoadingParameters | None,
) -> Trajectory | ShnitselDB | None:
    """Function to apply some final postprocessing common to all input routines.

    Args:
        dataset (xr.Dataset | ShnitselDB | None): The dataset, ShnitselDB object or None, if reading failed. Only individual trajectories will be post-processed here.
        loading_parameters (LoadingParameters | None): Parameters to set some defaults.

    Returns:
        Trajectory | None: _description_
    """
    if dataset is not None:
        # logging.debug(f"Finalizing: {repr(dataset)}")
        if isinstance(dataset, xr.Dataset):
            # TODO: FIXME: use loading_parameters to configure state names
            dataset = set_state_defaults(dataset, loading_parameters)

            unset_vars = []
            for var in dataset.variables:
                if is_variable_assigned(dataset[var]):
                    # Remove tags
                    del dataset[var].attrs["__assigned"]
                else:
                    unset_vars.append(var)

            logging.debug(f"Dropping unset variables: {unset_vars}")
            dataset = dataset.drop_vars(unset_vars)

            return wrap_trajectory(convert_all_units_to_shnitsel_defaults(dataset))
        elif isinstance(dataset, ShnitselDB):
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
