from typing import TypeAlias
from shnitsel.data.proxy_class import Proxy
import xarray as xr


Trajectory: TypeAlias = xr.Dataset


class _Trajectory(Proxy):
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


def wrap_trajectory(ds: Trajectory | xr.DataTree) -> Trajectory | xr.DataTree:
    """Function to wrap a dataset or tree of datasets in a trajectory format.

    Used to hide the actual type logic of the Trajectory wrapper.

    Args:
        ds ( Trajectory|xr.DataTree): The dataset to wrap in a Trajectory instance. If provided a tree, all datasets in the tree will be wrapped.

    Returns:
        Trajectory|xr.DataTree: The dataset wrapped in a Trajectory object or the original Trajectory instance. Alternatively, the DataTree with each individual dataset wrapped.
    """

    if isinstance(ds, xr.DataTree):
        return ds.map_over_datasets(func=_wrap_single_trajectory)
    else:
        return _wrap_single_trajectory(ds)


def _wrap_single_trajectory(ds: Trajectory) -> Trajectory:
    """Function to wrap a single dataset in a trajectory format.

    Used to hide the actual type logic of the Trajectory wrapper.

    Args:
        ds (xr.Dataset | Trajectory): The dataset or trajectory to (potentially) wrap in a Trajectory instance.

    Returns:
        Trajectory: The dataset wrapped in a Trajectory object or the original Trajectory instance.
    """

    if "__original_dataset" not in ds.attrs:
        ds.attrs["__original_dataset"] = ds.copy(deep=True)
    return ds
