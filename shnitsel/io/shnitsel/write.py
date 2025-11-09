import logging
import os
import numpy as np

from shnitsel.data.trajectory_format import Trajectory
import xarray as xr
import json

from shnitsel.io.helpers import PathOptionsType


def write_shnitsel_file(
    dataset: xr.Dataset | Trajectory, savepath: PathOptionsType, complevel: int = 9
):
    """Function to write a trajectory in Shnitsel format (xr.) to a ntcdf hdf5 file format.

    Strips all internal attributes first to avoid errors during writing.
    When writing directly with to_netcdf, errors might occur due to internally set attributes with problematic types.

    Args:
        dataset (xr.Dataset | Trajectory): The dataset or trajectory to write (omit if using accessor)
        savepath (PathOptionsType): The path at which to save the trajectory file

    Returns:
        Unknown: Returns the result of the final call to xr.Dataset.to_netcdf()
    """
    """Save a ``Dataset``, presumably (but not necessarily) consisting of frames of trajectories, to a file at ``path``.

    Parameters
    ----------
    frames (omit if using accessor)
        The ``Dataset`` to save
    path
        The path at which to save it
    complevel, optional
        The level of ``gzip`` compression which will be applied to all variables in the ``Dataset``, by default 9

    Notes
    -----
    This function/accessor method wraps :py:meth:`xarray.Dataset.to_netcdf` but not :py:func:`numpy.any`.
    """
    cleaned_ds = dataset.copy()  # Shallow copy to avoid adding attrs etc. to original
    encoding = {
        var: {"compression": "gzip", "compression_opts": complevel}
        for var in cleaned_ds
    }

    # NetCDF does not support booleans
    for data_var in cleaned_ds.data_vars:
        if np.issubdtype(cleaned_ds.data_vars[data_var].dtype, np.bool_):
            cleaned_ds = cleaned_ds.assign(
                {data_var: cleaned_ds.data_vars[data_var].astype('i1')}
            )
    for coord in cleaned_ds.coords:
        if np.issubdtype(cleaned_ds.coords[coord].dtype, np.bool_):
            cleaned_ds = cleaned_ds.assign_coords(
                {coord: cleaned_ds.coords[coord].astype('i1')}
            )

    # NetCDF does not support MultiIndex
    # Keep a record of the level names in the attrs
    midx_names = []
    for name, index in cleaned_ds.indexes.items():
        if index.name == name and len(index.names) > 1:
            midx_names.append(name)
            midx_levels = list(index.names)
            cleaned_ds.attrs[f'_MultiIndex_levels_for_{name}'] = midx_levels
    cleaned_ds.attrs['_MultiIndex_levels_from_attrs'] = 1

    def ndarray_to_json_ser(value):
        return {"__ndarray:": {"entries": value.tolist(), "dtype": value.dtype.descr}}

    remove_attrs = []

    for attr in cleaned_ds.attrs:
        # Strip internal attributes
        if str(attr).startswith("__"):
            # logging.debug(f"Mark for removing {attr}")
            remove_attrs.append(attr)
        else:
            value = cleaned_ds.attrs[attr]
            if isinstance(value, np.ndarray):
                value = ndarray_to_json_ser(value)
            cleaned_ds.attrs[attr] = json.dumps(value)

    for attr in remove_attrs:
        del cleaned_ds.attrs[attr]
        logging.debug(f"Stripping attribute {attr}")

    for data_var in cleaned_ds.variables:
        # If we delete while iterating, an error will occur.
        remove_attrs = []
        for attr in cleaned_ds[data_var].attrs:
            # Strip internal attributes
            if str(attr).startswith("__"):
                # logging.debug(f"Mark for removing {data_var}.{attr}")
                remove_attrs.append(attr)
            else:
                value = cleaned_ds[data_var].attrs[attr]
                if isinstance(value, np.ndarray):
                    value = ndarray_to_json_ser(value)
                cleaned_ds[data_var].attrs[attr] = json.dumps(value)
        for attr in remove_attrs:
            logging.debug(f"Stripping attribute {data_var}.{attr}")
            del cleaned_ds[data_var].attrs[attr]

        # if np.issubdtype(np.asarray(cleaned_ds.attrs[attr]).dtype, np.bool_):
        #    cleaned_ds.attrs[attr] = int(cleaned_ds.attrs[attr])

    cleaned_ds.attrs["__attrs_json_encoded"] = 1
    cleaned_ds.attrs["__shnitsel_format_version"] = "v1.0"

    return cleaned_ds.reset_index(midx_names).to_netcdf(
        savepath, engine='h5netcdf', encoding=encoding
    )
