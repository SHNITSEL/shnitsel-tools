import logging
import xarray as xr


def mark_variable_assigned(var: xr.DataArray) -> None:
    """
    Function to set a flag on the variable in the dataset to mark it as actually available and not just filled with default values.

    Should only be called on variables that had a non-default value assigned.
    Variables that have not been flagged, may be dropped upon finalization of the loading routine.

    Parameters
    ----------
    var : xr.DataArray
        The variable to set the flag on to mark it as assigned to.
    """
    var.attrs["__assigned"] = True


def mark_dataset_for_cleanup(ds: xr.Dataset) -> None:
    """
    Function to set a flag on the dataset to remove all unassigned variables

    Parameters
    ----------
    ds : xr.Dataset
        The variable to set the flag on to mark it as assigned to.
    """
    ds.attrs["__shnitsel_setup_for_cleanup"] = True


def is_marked_for_cleanup(ds: xr.Dataset) -> bool:
    """
    Function to set a flag on the dataset to remove all unassigned variables

    Parameters
    ----------
    ds : xr.Dataset
        The variable to check the cleanup flag on

    Returns
    -------
    Whether a cleanup flag is set on a dataset.
    """
    return ds.attrs.get("__shnitsel_setup_for_cleanup", False)


def is_variable_assigned(var: xr.DataArray) -> bool:
    """
    Function to check a flag on a variable in a dataset whether it has been assigned with actual values.


    Parameters
    ----------
    var : xr.DataArray
        The variable to check for a set "__assigned" flag.

    Returns
    -------
    bool
        Returns True if the variable has been marked as having been assigned a value previously.
    """
    return "__assigned" in var.attrs


def clean_unassigned_variables(dataset: xr.Dataset) -> xr.Dataset:
    """Helper function to clean all unassigned data variables from a dataset.

    To make sure we only ever clean data after import, a flag is checked in the `_shnitsel_setup_for_cleanup`
    attribute so that we do not accidentally clean out data from previously serialized data formats,
    where such a cleanup is not neede.
    The flag is normally set in `create_initial_dataset()` and not necessary to set manually.

    For all unflagged inputs, the return value will be the unchanged input parameter.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset to check for cleaning and potentially drop unassigned variables from.

    Returns
    -------
    xr.Dataset
        The dataset either unchanged if the flag for cleaning was not set or with unassigned variables
        left over from the dummy setup process being removed.
    """

    unset_vars = []
    for var in dataset.variables:
        if is_variable_assigned(dataset[var]):
            # Remove tags
            del dataset[var].attrs["__assigned"]
        else:
            unset_vars.append(var)

    if is_marked_for_cleanup(dataset):
        logging.debug(f"Dropping unset variables: {unset_vars}")
        dataset = dataset.drop_vars(unset_vars)
        del dataset.attrs["__shnitsel_setup_for_cleanup"]

    return dataset
