from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.data.trajectory_format import Trajectory


def wrap_trajectory(ds: Trajectory | ShnitselDB) -> Trajectory | ShnitselDB:
    """Function to wrap a dataset or tree of datasets in a trajectory format.

    Used to hide the actual type logic of the Trajectory wrapper.

    Args:
        ds ( Trajectory|ShnitselDB): The dataset to wrap in a Trajectory instance. If provided a ShnitselDB tree, all datasets in the tree will be wrapped.

    Returns:
        Trajectory|ShnitselDB: The dataset wrapped in a Trajectory object or the original Trajectory instance. Alternatively, the DataTree with each individual dataset wrapped.
    """

    if isinstance(ds, ShnitselDB):
        # TODO: FIXME: For the time being, I have no simple solution to getting a ShnitselDB back if we
        # apply a map to the tree because all attributes will be lost.
        return ds  # .map_over_datasets(func=_wrap_single_trajectory)
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
