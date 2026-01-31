import json
import logging
import os
import traceback
from typing import Any, Callable, Dict, TypeVar
import numpy as np
import xarray as xr
import sys

from shnitsel.core._api_info import internal
from shnitsel.io.shared.trajectory_finalization import normalize_dataset
from shnitsel.io.shared.trajectory_setup import fill_missing_dataset_variables
from ...data.tree.datatree_level import _datatree_level_attribute_key

from shnitsel.io.shared.variable_flagging import mark_variable_assigned

from shnitsel.io.shared.helpers import LoadingParameters, PathOptionsType

# def open_frames(path):

T = TypeVar("T")


@internal()
def read_shnitsel_file(
    path: PathOptionsType,
    loading_parameters: LoadingParameters | None = None,
) -> xr.Dataset | xr.DataTree | None:
    """Opens a NetCDF4 file saved by shnitsel-tools, specially interpreting certain attributes.

    Parameters
    ----------
    path : PathOptionsType
        The path of the file to open.
    loading_parameters : LoadingParameters, optional
        Parameter settings for e.g. standard units or state names.

    Returns
    -------
    xr.Dataset | xr.DataTree | None
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
        frames = xr.open_datatree(path)

        # Unpack the dataset if the file did not contain a tree
        if (
            _datatree_level_attribute_key not in frames.attrs
            and (
                "_shnitsel_tree_indicator" not in frames.attrs
                or frames.attrs["_shnitsel_tree_indicator"] != "TREE"
            )
            and len(frames.children) == 0
        ):
            # We have a simple data type, unwrap the tree.
            # logging.debug(f"{frames.attrs=}")
            frames = frames.dataset
    except ValueError as ds_err:
        dataset_info = sys.exc_info()
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        else:
            try:
                frames = xr.open_dataset(path)
            except ValueError as dt_err:
                datatree_info = sys.exc_info()
                message = "Failed to load file as either Dataset or DataTree: %(ds_err)s \n %(ds_info)s \n %(dt_err)s \n %(dt_info)s"
                params = {
                    "ds_err": ds_err,
                    "ds_info": dataset_info,
                    "dt_err": dt_err,
                    "dt_info": datatree_info,
                }
                logging.error(message, params)
                raise ValueError(message % params)

    if "__shnitsel_format_version" in frames.attrs:
        shnitsel_format_version = frames.attrs["__shnitsel_format_version"]
        del frames.attrs["__shnitsel_format_version"]
    else:
        shnitsel_format_version = "v1.0"

    if shnitsel_format_version in _SHNITSEL_READERS:
        return _SHNITSEL_READERS[shnitsel_format_version](frames, loading_parameters)
    else:
        message = (
            "Attempting to load a shnitsel file with unknown format %(format_version)s. \n"
            "This file might have been created with a later version of the `shnitsel-tools` package. \n"
            "Please update the `shnitsel-tools` package and attempt to read the file again."
        )
        params = {"format_version": {shnitsel_format_version}}
        logging.error(message, params)
        raise ValueError(message % params)


def _parse_shnitsel_file_v1_0(
    frames: T, loading_parameters: LoadingParameters | None = None
) -> T:
    """Internal function to do a best-effort attempt to load the original shnitsel file format.

    Will print a warning that you should be using shnitsel v1.1 files to have full type information.

    Parameters
    ----------
    frames : xr.Dataset
        The loaded Dataset from the netcdf file that needs to be post-processed.
    loading_parameters : LoadingParameters | None, optional
        Optional loading parameters for setting units. Defaults to None.

    Returns
    -------
    xr.Dataset
        The post-processed shnitsel trajectory
    """
    if not isinstance(frames, xr.Dataset):
        raise ValueError(
            "A version 1.0 shnitsel file can only contain xr.Dataset entries."
        )

    logging.warning(
        "You are opening a Shnitsel file of format v1.0. This format did not contain full unit information for all observables. \n"
        "You should either regenerate the shnitsel file from the input data with a later version of the shnitsel-tools package or attempt to retrieve a later version of the file."
    )

    # Rename time coordinate to same name everywhere
    tcoord = None
    if "time" in frames.coords:
        tcoord = "time"
        if "units" not in frames.time.attrs:
            logging.warning("Guessing time dimension to be unitless. (units=`1`)")
            frames.time.attrs.update({"units": "1", "unitdim": "time"})
    elif "ts" in frames.coords:
        logging.info(
            "Renaming 'ts' dimension to 'time' to make trajectory conform to standard shnitsel format."
        )
        frames = frames.rename({"ts": "time"})
        if "delta_t" in frames.attrs:
            frames.time.values *= float(frames.attrs["delta_t"])
            logging.warning("Guessing time unit to be `fs`")
            frames.time.attrs.update({"units": "fs", "unitdim": "time"})
        else:
            logging.warning("Guessing time dimension to be unitless. (units=`1`)")
            frames.time.attrs.update({"units": "1", "unitdim": "time"})
        tcoord = "time"

    # Restore MultiIndexes
    indicator = "_MultiIndex_levels_from_attrs"
    level_prefix = "_MultiIndex_levels_for_"
    if frames.attrs.get(indicator, False):
        # New way: get level names from attrs
        del frames.attrs[indicator]
        for k, v in frames.attrs.items():
            if k.startswith(level_prefix):
                index_name = k[len(level_prefix) :]
                # print(f"Index {index_name=} : levels={v}")
                if len(set([frames.coords[level_name].size for level_name in v])) == 1:
                    # all levels have the same length:
                    frames = frames.set_xindex(v)
                    mark_variable_assigned(frames[index_name])

                del frames.attrs[k]
    else:
        # Old way: hardcoded level names

        if tcoord is not None and tcoord in frames.coords and "trajid" in frames:
            # Clear existing indexes before setting the multi-index
            if tcoord in frames.indexes:
                frames = frames.reset_index(tcoord)
            if "trajid" in frames.indexes:
                frames = frames.reset_index("trajid")

            frames = frames.set_xindex(["trajid", tcoord])

        if "from" in frames.coords and "to" in frames.coords:
            if "from" in frames.indexes:
                frames = frames.reset_index("from")
            if "to" in frames.indexes:
                frames = frames.reset_index("to")
            frames = frames.set_xindex(["from", "to"])

    return normalize_dataset(fill_missing_dataset_variables(frames))


def _decode_shnitsel_v1_1_dataset(dataset: xr.Dataset) -> xr.Dataset:
    """Function to decode encoded attributes and MultiIndices of a dataset

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset to process

    Returns
    -------
    xr.Dataset
        Copy of the dataset with attributes and MultiIndex instances decoded
    """
    if "__attrs_json_encoded" in dataset.attrs:
        del dataset.attrs["__attrs_json_encoded"]

        def json_deserialize_ndarray(key: str, value: str) -> Any:
            if key.startswith("__") or not isinstance(value, str):
                return value
            try:
                value_d = json.loads(value)
                if isinstance(value_d, dict):
                    if "__ndarray" in value_d:
                        config = value_d["__ndarray"]

                        entries = config["entries"]
                        dtype_descr = np.dtype([tuple(i) for i in config["dtype"]])

                        value_d = np.array(entries, dtype=dtype_descr)
                return value_d
            except TypeError as e:
                params = {"e": e, "traceback": traceback.format_exc()}
                logging.debug(
                    "Encountered error during json decode: %(e)s \n %(traceback)s",
                    params,
                )
                return value

        for attr in dataset.attrs:
            dataset.attrs[attr] = json_deserialize_ndarray(
                str(attr), dataset.attrs[attr]
            )

        for data_var in dataset.variables:
            # Mark variables as assigned so they are not stripped in finalization
            mark_variable_assigned(dataset[data_var])
            for attr in dataset[data_var].attrs:
                dataset[data_var].attrs[attr] = json_deserialize_ndarray(
                    str(attr), dataset[data_var].attrs[attr]
                )

    # Rename time coordinate to same name everywhere
    tcoord = None
    if "time" in dataset.coords:
        tcoord = "time"
        if "units" not in dataset.time.attrs:
            logging.warning("Guessing time dimension to be unitless. (units=`1`)")
            dataset.time.attrs.update({"units": "1", "unitdim": "time"})
    elif "ts" in dataset.coords:
        logging.info(
            "Renaming 'ts' dimension to 'time' to make trajectory conform to standard shnitsel format."
        )
        frames = dataset.rename({"ts": "time"})
        if "delta_t" in frames.attrs:
            frames.time.values *= float(frames.attrs["delta_t"])
            logging.warning("Guessing time unit to be `fs`")
            frames.time.attrs.update({"units": "fs", "unitdim": "time"})
        else:
            logging.warning("Guessing time dimension to be unitless. (units=`1`)")
            frames.time.attrs.update({"units": "1", "unitdim": "time"})
        tcoord = "time"

    # Restore MultiIndexes
    indicator = "_MultiIndex_levels_from_attrs"
    level_prefix = "_MultiIndex_levels_for_"
    if dataset.attrs.get(indicator, False):
        # New way: get level names from attrs
        del dataset.attrs[indicator]
        for k, v in dataset.attrs.items():
            if k.startswith(level_prefix):
                index_name = k[len(level_prefix) :]
                # print(f"Index {index_name=} : levels={v}")
                if len(set([dataset.coords[level_name].size for level_name in v])) == 1:
                    # all levels have the same length:
                    dataset = dataset.set_xindex(v)
                    mark_variable_assigned(dataset[index_name])

                del dataset.attrs[k]
    else:
        # Old way: hardcoded level names
        if tcoord is not None and tcoord in dataset.coords and "trajid" in dataset:
            # Clear existing indexes before setting the multi-index
            if tcoord in dataset.indexes:
                dataset = dataset.reset_index(tcoord)
            if "trajid" in dataset.indexes:
                dataset = dataset.reset_index("trajid")
            dataset = dataset.set_xindex(["trajid", tcoord])

        if "from" in dataset.coords and "to" in dataset.coords:
            if "from" in dataset.indexes:
                dataset = dataset.reset_index("from")
            if "to" in dataset.indexes:
                dataset = dataset.reset_index("to")
            dataset = dataset.set_xindex(["from", "to"])

    return normalize_dataset(fill_missing_dataset_variables(dataset))


XRAttrType = TypeVar("XRAttrType", bound=xr.Dataset | xr.DataArray | xr.DataTree)


def decode_attrs(obj: XRAttrType) -> XRAttrType:
    """Helper function to decode attributes of an object
    that may have been serialized during the writing process.

    Parameters
    ----------
    obj : xr.Dataset | xr.DataArray | xr.DataTree
        The object whose attributes should be decoded

    Returns
    -------
    xr.Dataset | xr.DataArray | xr.DataTree
        The same type as the input `obj` but with deserialized attributes
    """
    if "__attrs_json_encoded" in obj.attrs:
        del obj.attrs["__attrs_json_encoded"]

        def json_deserialize_ndarray(key: str, value: str | Any) -> Any:
            if key.startswith("__") or not isinstance(value, str):
                return value
            try:
                value_d = json.loads(value)
                if isinstance(value_d, dict):
                    if "__ndarray" in value_d:
                        config = value_d["__ndarray"]

                        entries = config["entries"]
                        dtype_descr = np.dtype([tuple(i) for i in config["dtype"]])

                        value_d = np.array(entries, dtype=dtype_descr)
                return value_d
            except TypeError as e:
                params = {"e": e, "traceback": traceback.format_exc()}
                logging.debug(
                    "Encountered error during json decode: %(e)s \n %(traceback)s",
                    params,
                )
                return value

        for attr in obj.attrs:
            obj.attrs[attr] = json_deserialize_ndarray(str(attr), obj.attrs[attr])

    return obj


def _decode_shnitsel_v1_1_datatree(datatree: xr.DataTree) -> xr.DataTree:
    """Decoder for v1.1 versions of shnitsel datatree formats.

    Parameters
    ----------
    datatree : xr.DataTree
        The Datatree as read from the file.

    Returns
    -------
    xr.DataTree
        The datatree after initial processing and decoding of attributes
    """
    res = datatree.copy()
    if res.has_data:
        res.dataset = _decode_shnitsel_v1_1_dataset(res.dataset)

    return res.assign(
        {k: _decode_shnitsel_v1_1_datatree(v) for k, v in res.children.items()}
    )


def _decode_shnitsel_v1_2_datatree(datatree: xr.DataTree) -> xr.DataTree:
    """Decoder for v1.2 versions of shnitsel datatree formats.

    Parameters
    ----------
    datatree : xr.DataTree
        The Datatree as read from the file.

    Returns
    -------
    xr.DataTree
        The datatree after initial processing and decoding of attributes
    """
    res = datatree.copy()
    if res.has_data:
        res.dataset = _decode_shnitsel_v1_1_dataset(res.dataset)

    res = decode_attrs(res)

    return res.assign(
        {k: _decode_shnitsel_v1_2_datatree(v) for k, v in res.children.items()}
    )


def _parse_shnitsel_file_v1_1(
    frames: T, loading_parameters: LoadingParameters | None = None
) -> T:
    """Internal function to parse the revised shnitsel format v1.1 with better attribute encoding and more extensive unit declarations.

    Parameters
    ----------
    frames : xr.Dataset|xr.DataTree
        The loaded Dataset or tree from the netcdf file that needs to be post-processed.
    loading_parameters : LoadingParameters | None, optional
        Optional loading parameters for setting units. Defaults to None.

    Returns
    -------
    xr.Dataset
        The post-processed shnitsel trajectory
    """
    from shnitsel.data.tree import (
        complete_shnitsel_tree,
    )

    if not isinstance(frames, xr.Dataset) and not isinstance(frames, xr.DataTree):
        raise ValueError(
            "A version 1.1 shnitsel file can only contain xr.Dataset or xr.DataTree entries."
        )
    if isinstance(frames, xr.DataTree):
        # import pprint
        decoded_tree = _decode_shnitsel_v1_1_datatree(frames)
        # shnitsel_db = build_shnitsel_db(frames)
        # pprint.pprint(shnitsel_db)

        # Decode json encoded attributes if json encoding is recorded
        # return shnitsel_db.map_over_trajectories(_decode_shnitsel_v1_1_dataset)  # type: ignore
        return decoded_tree
    if isinstance(frames, xr.Dataset):
        # Decode json encoded attributes if json encoding is recorded
        return _decode_shnitsel_v1_1_dataset(frames)


def _parse_shnitsel_file_v1_2(
    frames: T, loading_parameters: LoadingParameters | None = None
) -> T:
    """Internal function to parse the revised shnitsel format v1.2 with better attribute encoding and more extensive unit declarations using the DataTree format.

    Parameters
    ----------
    frames : xr.DataTree
        The loaded Dataset from the netcdf file that needs to be post-processed.
    loading_parameters : LoadingParameters | None, optional
        Optional loading parameters for setting units. Defaults to None.

    Returns
    -------
    xr.DataTree
        The post-processed shnitsel trajectory
    """
    from shnitsel.data.tree import (
        complete_shnitsel_tree,
    )

    if not isinstance(frames, xr.Dataset) and not isinstance(frames, xr.DataTree):
        raise ValueError(
            "A version 1.2 shnitsel file can only contain xr.Dataset or xr.DataTree entries."
        )
    if isinstance(frames, xr.DataTree):
        # import pprint

        # pprint.pprint(frames)
        decoded_tree = _decode_shnitsel_v1_2_datatree(frames)
        # pprint.pprint(shnitsel_db)
        # decoded_tree = build_shnitsel_db(decoded_tree)
        # pprint.pprint(shnitsel_db)
        return decoded_tree
    if isinstance(frames, xr.Dataset):
        # Decode json encoded attributes if json encoding is recorded
        return _decode_shnitsel_v1_1_dataset(frames)


_SHNITSEL_READERS: Dict[
    str,
    Callable[
        [xr.Dataset | xr.DataTree, LoadingParameters | None], xr.Dataset | xr.DataTree
    ],
] = {
    "v1.0": _parse_shnitsel_file_v1_0,
    "v1.1": _parse_shnitsel_file_v1_1,
    "v1.2": _parse_shnitsel_file_v1_2,
    "v1.3": _parse_shnitsel_file_v1_2,  # it's the 1.2 format but with more extensive in-tree metdatata
}
