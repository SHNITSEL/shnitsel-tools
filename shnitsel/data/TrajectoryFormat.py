import xarray as xr


class Trajectory:
    """Class to wrap read trajectory information in a shnitsel-conform format.

    Used to keep track of original data while the trajectory data is modified
    """

    def __init__(self, initial_ds: xr.Dataset) -> None:
        self.original_dataset = initial_ds
        self.current_derived_dataset = initial_ds.copy(deep=True)

    def get_current_raw(self) -> xr.Dataset:
        return self.current_derived_dataset

    def get_original_raw(self) -> xr.Dataset:
        return self.original_dataset
