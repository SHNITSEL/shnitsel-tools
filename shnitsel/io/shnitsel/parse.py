import json
import logging
import os
import pathlib
from typing import Any, Callable, Dict, TypeVar
import numpy as np
import xarray as xr
import sys

from shnitsel.data.shnitsel_db_format import (
    build_shnitsel_db,
)
from shnitsel.io.helpers import LoadingParameters, PathOptionsType

# def open_frames(path):

T = TypeVar('T')


def read_shnitsel_file(
    path: PathOptionsType, loading_parameters: LoadingParameters | None = None
) -> xr.Dataset | xr.DataTree:
    """Opens a NetCDF4 file saved by shnitsel-tools, specially interpreting certain attributes.

    Parameters
    ----------
    path (PathOptionsType):
        The path of the file to open.
    loading_parameters (LoadingParameters,optional): Parameter settings for e.g. standard units or state names.

    Returns
    -------
        An :py:class:`xarray.Dataset` with any MultiIndex restored.
        A :py:class:`ShnitselDB` with any MultiIndex restored and attributes decoded.

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
    except ValueError as ds_err:
        dataset_info = sys.exc_info()
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        else:
            try:
                frames = xr.open_datatree(path)
            except ValueError as dt_err:
                datatree_info = sys.exc_info()
                message = f"Failed to load file as either Dataset or DataTree: {ds_err} \n{dataset_info}\n {dt_err}\n {datatree_info}"
                logging.error(message)
                raise ValueError(message)

    if "__shnitsel_format_version" in frames.attrs:
        shnitsel_format_version = frames.attrs["__shnitsel_format_version"]
        del frames.attrs["__shnitsel_format_version"]
    else:
        shnitsel_format_version = "v1.0"

    if shnitsel_format_version in __SHNITSEL_READERS:
        return __SHNITSEL_READERS[shnitsel_format_version](frames, loading_parameters)
    else:
        message = (
            f"Attempting to load a shnitsel file with unknown format {shnitsel_format_version}. \n"
            "This file might have been created with a later version of the `shnitsel-tools` package. \n"
            "Please update the `shnitsel-tools` package and attempt to read the file again."
        )
        logging.error(message)
        raise ValueError(message)


def _parse_shnitsel_file_v1_0(
    frames: T, loading_parameters: LoadingParameters | None = None
) -> T:
    """Internal function to do a best-effort attempt to load the original shnitsel file format.

    Will print a warning that you should be using shnitsel v1.1 files to have full type information.

    Args:
        frames (xr.Dataset): The loaded Dataset from the netcdf file that needs to be post-processed.
        loading_parameters (LoadingParameters | None, optional): Optional loading parameters for setting units. Defaults to None.

    Returns:
        xr.Dataset: The post-processed shnitsel trajectory
    """
    if not isinstance(frames, xr.Dataset):
        raise ValueError(
            "A version 1.0 shnitsel file can only contain xr.Dataset entries."
        )

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
            frames = frames.rename({"ts": "time"})
            tcoord = "time"

        if tcoord is not None:
            frames = frames.set_xindex(["trajid", tcoord])

        if "from" in frames.coords and "to" in frames.coords:
            frames = frames.set_xindex(["from", "to"])

    return frames


def _decode_shnitsel_v1_1_dataset(dataset: xr.Dataset) -> xr.Dataset:
    """Function to decode encoded attributes and MultiIndices of a dataset

    Args:
        dataset (xr.Dataset): The dataset to process

    Returns
        xr.Dataset: Copy of the dataset with attributes and MultiIndex instances decoded
    """
    if "__attrs_json_encoded" in dataset.attrs:
        del dataset.attrs["__attrs_json_encoded"]

        def json_deserialize_ndarray(value: str) -> Any:
            value_d = json.loads(value)
            if isinstance(value_d, dict):
                if "__ndarray" in value_d:
                    config = value_d["__ndarray"]

                    entries = config["entries"]
                    dtype_descr = np.dtype([tuple(i) for i in config["dtype"]])

                    value_d = np.array(entries, dtype=dtype_descr)
            return value_d

        for attr in dataset.attrs:
            dataset.attrs[attr] = json_deserialize_ndarray(dataset.attrs[attr])

        for data_var in dataset.variables:
            for attr in dataset[data_var].attrs:
                dataset[data_var].attrs[attr] = json_deserialize_ndarray(
                    dataset.attrs[attr]
                )

    # Restore MultiIndexes
    indicator = "_MultiIndex_levels_from_attrs"
    if dataset.attrs.get(indicator, False):
        # New way: get level names from attrs
        del dataset.attrs[indicator]
        for k, v in dataset.attrs.items():
            if k.startswith("_MultiIndex_levels_for_"):
                dataset = dataset.set_xindex(v)
                del dataset.attrs[k]
    else:
        # TODO: FIXME: rename time trajectory to same name everywhere
        # Old way: hardcoded level names
        tcoord = None
        if "time" in dataset.coords:
            tcoord = "time"
        elif "ts" in dataset.coords:
            logging.info(
                "Renaming 'ts' dimension to 'time' to make trajectory conform to standard shnitsel format."
            )
            dataset = dataset.rename({"ts": "time"})
            tcoord = "time"

        if tcoord is not None:
            dataset = dataset.set_xindex(["trajid", tcoord])

        if "from" in dataset.coords and "to" in dataset.coords:
            dataset = dataset.set_xindex(["from", "to"])

    return dataset


def _parse_shnitsel_file_v1_1(
    frames: T, loading_parameters: LoadingParameters | None = None
) -> T:
    """Internal function to parse the revised shnitsel format v1.1 with better attribute encoding and more extensive unit declarations.

    Args:
        frames (xr.Dataset|xr.DataTree): The loaded Dataset from the netcdf file that needs to be post-processed.
        loading_parameters (LoadingParameters | None, optional): Optional loading parameters for setting units. Defaults to None.

    Returns:
        xr.Dataset: The post-processed shnitsel trajectory
    """
    if not isinstance(frames, xr.Dataset) and not isinstance(frames, xr.DataTree):
        raise ValueError(
            "A version 1.1 shnitsel file can only contain xr.Dataset or xr.DataTree entries."
        )
    if isinstance(frames, xr.Dataset):
        # Decode json encoded attributes if json encoding is recorded
        return _decode_shnitsel_v1_1_dataset(frames)
    elif isinstance(frames, xr.DataTree):
        import pprint

        # Decode json encoded attributes if json encoding is recorded
        decoded_frames = frames.map_over_datasets(_decode_shnitsel_v1_1_dataset)

        pprint.pprint(decoded_frames)

        return build_shnitsel_db(frames)


__SHNITSEL_READERS: Dict[
    str,
    Callable[
        [xr.Dataset | xr.DataTree, LoadingParameters | None], xr.Dataset | xr.DataTree
    ],
] = {
    "v1.0": _parse_shnitsel_file_v1_0,
    "v1.1": _parse_shnitsel_file_v1_1,
}
