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
