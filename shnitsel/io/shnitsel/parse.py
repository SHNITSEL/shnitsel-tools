import json
import logging
import os
import pathlib
import xarray as xr

from shnitsel.io.helpers import LoadingParameters, PathOptionsType

# TODO: We probably need a fallback version for old and new shnitsel file formats

# def open_frames(path):


def read_shnitsel_file(
    path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Opens a NetCDF4 file saved by shnitsel-tools, specially interpreting certain attributes.

    Parameters
    ----------
    path (PathOptionsType):
        The path of the file to open.
    loading_parameters (LoadingParameters,optional): Parameter settings for e.g. standard units or state names.

    Returns
    -------
        An :py:class:`xarray.Dataset` with any MultiIndex restored.

    Raises
    ------
    FileNotFoundError
        If there is is nothing at ``path``, or ``path`` is not a file.
    ValueError (or other exception)
        Raised by the underlying `h5netcdf <https://h5netcdf.org/>`_ engine if the file is corrupted.
    """
    # TODO: FIXME: use loading_parameters to configure units and state names
    # The error raised for a missing file can be misleading
    try:
        frames = xr.open_dataset(path)
    except ValueError as err:
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        else:
            raise err

    # Decode json encoded attributes if json encoding is recorded
    if "__attrs_json_encoded" in frames.attrs:
        del frames.attrs["__attrs_json_encoded"]

        for attr in frames.attrs:
            frames.attrs[attr] = json.loads(frames.attrs[attr])

        for data_var in frames.variables:
            for attr in frames[data_var].attrs:
                frames[data_var].attrs[attr] = json.loads(
                    frames.attrs[attr])

    # Restore MultiIndexes
    indicator = "_MultiIndex_levels_from_attrs"
    if frames.attrs.get(indicator, False):
        # New way: get level names from attrs
        del frames.attrs[indicator]
        for k, v in frames.attrs.items():
            if k.startswith("_MultiIndex_levels_for_"):
                frames = frames.set_xindex(v)
                del frames.attrs[k]
    else:
        # TODO: FIXME: rename time trajectory to same name everywhere
        # Old way: hardcoded level names
        tcoord = None
        if "time" in frames.coords:
            tcoord = "time"
        elif "ts" in frames.coords:
            logging.info(
                "Renaming 'ts' dimension to 'time' to make trajectory conform to standard shnitsel format."
            )
            frames.rename({"ts": "time"})
            tcoord = "time"

        if tcoord is not None:
            frames = frames.set_xindex(["trajid", tcoord])

        if "from" in frames.coords and "to" in frames.coords:
            frames = frames.set_xindex(["from", "to"])

    return frames
