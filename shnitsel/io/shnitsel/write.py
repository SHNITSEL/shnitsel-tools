
import os
import numpy as np


def write_shnitsel_file(frames, path: os.PathLike, complevel=9):
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
    frames = frames.copy()  # Shallow copy to avoid adding attrs etc. to original
    encoding = {
        var: {"compression": "gzip", "compression_opts": complevel} for var in frames
    }

    # NetCDF does not support booleans
    for data_var in frames.data_vars:
        if np.issubdtype(frames.data_vars[data_var].dtype, np.bool_):
            frames = frames.assign(
                {data_var: frames.data_vars[data_var].astype('i1')})
    for coord in frames.coords:
        if np.issubdtype(frames.coords[coord].dtype, np.bool_):
            frames = frames.assign_coords(
                {coord: frames.coords[coord].astype('i1')})
    for attr in frames.attrs:
        if np.issubdtype(np.asarray(frames.attrs[attr]).dtype, np.bool_):
            frames.attrs[attr] = int(frames.attrs[attr])

    # NetCDF does not support MultiIndex
    # Keep a record of the level names in the attrs
    midx_names = []
    for name, index in frames.indexes.items():
        if index.name == name and len(index.names) > 1:
            midx_names.append(name)
            midx_levels = list(index.names)
            frames.attrs[f'_MultiIndex_levels_for_{name}'] = midx_levels
    frames.attrs['_MultiIndex_levels_from_attrs'] = 1
    frames.reset_index(midx_names).to_netcdf(
        path, engine='h5netcdf', encoding=encoding)
