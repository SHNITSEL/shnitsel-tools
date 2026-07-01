from typing import Any

from shnitsel.bridges import to_xyz, traj_to_xyz
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.xr_io_compatibility import SupportsToXrConversion
from shnitsel.io.shared.helpers import PathOptionsType
import xarray as xr
from shnitsel.core._api_info import API


@API()
def write_xyz(
    dataset: xr.Dataset
    | xr.DataArray
    | SupportsToXrConversion
    | TreeNode[Any, xr.Dataset | xr.DataArray | SupportsToXrConversion],
    savepath: PathOptionsType,
):
    """Function to write structures to an xyz file

    Parameters
    ----------
    dataset : xr.Dataset | xr.DataArray | SupportsToXrConversion | TreeNode[Any, xr.Dataset  |  xr.DataArray  |  SupportsToXrConversion]
        The data source from which to extract the geometries to be written
    savepath : PathOptionsType
        The path to write the output format to.

    Returns
    -------
    bool
        True upon success, False otherwise
    """
    with open(savepath, "w") as fout:

        def map_ds_write(
            ds: xr.DataArray | xr.Dataset | SupportsToXrConversion,
        ) -> None:
            if isinstance(ds, SupportsToXrConversion):
                _, ds, _ = ds.as_xr_dataset()
            if isinstance(ds, xr.DataArray):
                fout.write(to_xyz(ds))
            elif isinstance(ds, xr.Dataset):
                fout.write(traj_to_xyz(ds.atXYZ))
            else:
                logging.error("Unsupported dataset type: %s", type(ds))
            return None

        if isinstance(dataset, TreeNode):
            dataset.map_data(lambda ds: map_ds_write(ds), keep_empty_branches=False)
        else:
            map_ds_write(dataset)

    return True
