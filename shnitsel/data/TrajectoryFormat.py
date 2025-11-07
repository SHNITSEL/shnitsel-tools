from shnitsel.data.Proxy import Proxy
import xarray as xr


class Trajectory(Proxy):
    """Class to wrap trajectory information in a shnitsel-conform format.

    Used to keep track of original data while the trajectory data is modified
    """

    def __init__(self, initial_ds: xr.Dataset) -> None:
        super().__init__(initial_ds)
        self.attrs["__original_dataset"] = initial_ds.copy(deep=True)

    def get_current_raw(self) -> xr.Dataset:
        return object.__getattribute__(self, "_obj")

    def get_original_raw(self) -> xr.Dataset:
        # TODO: Make Proxy wrapper functions return wrappers as well
        return self.attrs["__original_dataset"]


def wrap_trajectory(ds: xr.Dataset | Trajectory) -> Trajectory:
    """Function to wrap a dataset in a trajectory format.

    Used to hide the actual type logic of the Trajectory wrapper.

    Args:
        ds (xr.Dataset | Trajectory): The dataset or trajectory to (potentially) wrap in a Trajectory instance.

    Returns:
        Trajectory: The dataset wrapped in a Trajectory object or the original Trajectory instance.
    """
    if isinstance(ds, Trajectory):
        return ds
    else:
        return Trajectory(ds)
