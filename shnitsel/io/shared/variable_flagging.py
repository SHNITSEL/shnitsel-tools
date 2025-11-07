import xarray as xr


def mark_variable_assigned(var: xr.DataArray) -> None:
    """Function to set a flag on the variable in the dataset to mark it as actually available and not just filled with default values.

    Should only be called on variables that had a non-default value assigned.
    Variables that have not been flagged, may be dropped upon finalization of the loading routine.

    Args:
        var (xr.DataArray): The variable to set the flag on to mark it as assigned to.
    """
    var.attrs["__assigned"] = True


def is_variable_assigned(var: xr.DataArray) -> bool:
    """Function to check a flag on a variable in a dataset whether it has been assigned with actual values.

    Args:
        var (xr.DataArray): The variable to check for a set "__assigned" flag.
    """
    return "__assigned" in var.attrs
