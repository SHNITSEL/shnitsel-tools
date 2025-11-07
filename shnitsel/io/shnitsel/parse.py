import json
import logging
import os
import pathlib
from typing import Any, Callable, Dict
import numpy as np
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

    if "__shnitsel_format_version" in frames.attrs:
        shnitsel_format_version = frames.attrs["__shnitsel_format_version"]
        del frames.attrs["__shnitsel_format_version"]
    else:
        shnitsel_format_version = "v1.0"

    if shnitsel_format_version in __SHNITSEL_READERS:
        return __SHNITSEL_READERS[shnitsel_format_version](frames, loading_parameters)
    else:
        logging.error(
            f"Attempting to load a shnitsel file with unknown format {shnitsel_format_version}. \n"
            "This file might have been created with a later version of the `shnitsel-tools` package. \n"
            "Please update the `shnitsel-tools` package and attempt to read the file again."
        )


def _parse_shnitsel_file_v1_0(
    frames: xr.Dataset, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Internal function to do a best-effort attempt to load the original shnitsel file format.

    Will print a warning that you should be using shnitsel v1.1 files to have full type information.

    Args:
        frames (xr.Dataset): The loaded Dataset from the netcdf file that needs to be post-processed.
        loading_parameters (LoadingParameters | None, optional): Optional loading parameters for setting units. Defaults to None.

    Returns:
        xr.Dataset: The post-processed shnitsel trajectory
    """
    logging.warning(
        f"You are opening a Shnitsel file of format v1.0. This format did not contain full unit information for all observables. \n"
        f"You should either regenerate the shnitsel file from the input data with a later version of the shnitsel-tools package or attempt to retrieve a later version of the file."
    )

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


def _parse_shnitsel_file_v1_1(
    frames: xr.Dataset, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset:
    """Internal function to parse the revised shnitsel format v1.1 with better attribute encoding and more extensive unit declarations.

    Args:
        frames (xr.Dataset): The loaded Dataset from the netcdf file that needs to be post-processed.
        loading_parameters (LoadingParameters | None, optional): Optional loading parameters for setting units. Defaults to None.

    Returns:
        xr.Dataset: The post-processed shnitsel trajectory
    """
    # Decode json encoded attributes if json encoding is recorded
    if "__attrs_json_encoded" in frames.attrs:
        del frames.attrs["__attrs_json_encoded"]

        def json_deserialize_ndarray(value: str) -> Any:
            value_d = json.loads(value)
            if isinstance(value_d, dict):
                if "__ndarray" in value_d:
                    config = value_d["__ndarray"]

                    entries = config["entries"]
                    dtype_descr = np.dtype([tuple(i) for i in config["dtype"]])

                    value_d = np.array(entries, dtype=dtype_descr)
            return value_d

        for attr in frames.attrs:
            frames.attrs[attr] = json_deserialize_ndarray(frames.attrs[attr])

        for data_var in frames.variables:
            for attr in frames[data_var].attrs:
                frames[data_var].attrs[attr] = json_deserialize_ndarray(
                    frames.attrs[attr]
                )

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


__SHNITSEL_READERS: Dict[
    str, Callable[[xr.Dataset, LoadingParameters | None], xr.Dataset]
] = {
    "v1.0": _parse_shnitsel_file_v1_0,
    "v1.1": _parse_shnitsel_file_v1_1,
}
