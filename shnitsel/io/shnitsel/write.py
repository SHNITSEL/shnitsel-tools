import logging
import pathlib
from typing import Any, Dict, Hashable
import numpy as np

import xarray as xr
import json

from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.xr_conversion import (
    data_to_xarray_dataset,
    tree_to_xarray_datatree,
)
from shnitsel.io.shared.helpers import PathOptionsType, make_uniform_path
from shnitsel.data.xr_io_compatibility import SupportsToXrConversion


class NumpyDataEncoder(json.JSONEncoder):
    """Custom encoder for numpy data types"""

    def default(self, o):
        if isinstance(o, np.integer):  # (np.int_, np.intc, np.intp, np.int8,
            # np.int16, np.int32, np.int64, np.uint8,
            # np.uint16, np.uint32, np.uint64)):

            return int(o)

        elif isinstance(
            o, np.floating
        ):  # (np.float_, np.float16, np.float32, np.float64)):
            return float(o)

        elif isinstance(o, (np.complexfloating, np.complex128)):
            return {'real': o.real, 'imag': o.imag}

        elif isinstance(o, (np.ndarray,)):
            return o.tolist()

        elif isinstance(o, (np.bool_)):
            return bool(o)

        elif isinstance(o, (np.void)):
            return None

        return json.JSONEncoder.default(self, o)


def ndarray_to_json_ser(value):
    return {"__ndarray": {"entries": value.tolist(), "dtype": value.dtype.descr}}


def _prepare_datatree(node: xr.DataTree) -> xr.DataTree:
    cleaned_node = node.copy()

    if cleaned_node.has_data:
        cleaned_node.dataset = _prepare_dataset(cleaned_node.dataset)

    cleaned_node = cleaned_node.assign(
        {k: _prepare_datatree(v) for k, v in cleaned_node.children.items()}
    )

    if "__attrs_json_encoded" not in cleaned_node.attrs:
        encode_attrs(cleaned_node)

    cleaned_node.attrs["__attrs_json_encoded"] = 1
    return cleaned_node


def _prepare_dataset(dataset: xr.Dataset) -> xr.Dataset:
    """Function to prepare a dataset for encoding by re-encoding most of its attributes to account for types not supported by NetCDF.

    Also removed internal settings and re-encodes multi-indices

    Parameters
    ----------
    dataset : xr.Dataset
        Dataset to process

    Returns
    -------
    xr.Dataset
        A copy of the Dataset with internal attributes removed, attributes appropriately encoded and multi-indices re-encoded.
    """
    cleaned_ds = dataset.copy()  # Shallow copy to avoid adding attrs etc. to original

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
        if str(coord).startswith("__"):
            # Try and get rid of runtime-only data.
            cleaned_ds = cleaned_ds.drop_vars(coord)
            # TODO: FIXME: Convert `__mol` object into SMILES to restore later?

    # NetCDF does not support MultiIndex
    # Keep a record of the level names in the attrs
    midx_names = []
    for name, index in cleaned_ds.indexes.items():
        if index.name == name and len(index.names) > 1:
            midx_names.append(name)
            midx_levels = list(index.names)
            cleaned_ds.attrs[f'_MultiIndex_levels_for_{name}'] = midx_levels
    cleaned_ds.attrs['_MultiIndex_levels_from_attrs'] = 1

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
            try:
                cleaned_ds.attrs[attr] = json.dumps(value, cls=NumpyDataEncoder)
            except ValueError as e:
                print(f"ds.attrs['{attr}']={cleaned_ds.attrs[attr]} -> {e}")

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

                try:
                    cleaned_ds[data_var].attrs[attr] = json.dumps(
                        value, cls=NumpyDataEncoder
                    )
                except ValueError as e:
                    print(f"ds['{data_var}'].attrs['{attr}']={cleaned_ds[data_var].attrs[attr]} -> {e}")
        for attr in remove_attrs:
            logging.debug(f"Stripping attribute {data_var}.{attr}")
            del cleaned_ds[data_var].attrs[attr]

        # if np.issubdtype(np.asarray(cleaned_ds.attrs[attr]).dtype, np.bool_):
        #    cleaned_ds.attrs[attr] = int(cleaned_ds.attrs[attr])

    cleaned_ds.attrs["__attrs_json_encoded"] = 1
    return cleaned_ds.reset_index(midx_names)


def encode_attrs(obj):
    remove_attrs = []

    for attr in obj.attrs:
        # Strip internal attributes
        if str(attr).startswith("__"):
            # logging.debug(f"Mark for removing {attr}")
            remove_attrs.append(attr)
        else:
            value = obj.attrs[attr]
            if isinstance(value, np.ndarray):
                value = ndarray_to_json_ser(value)
            obj.attrs[attr] = json.dumps(value, cls=NumpyDataEncoder)

    for attr in remove_attrs:
        del obj.attrs[attr]
        logging.debug(f"Stripping attribute {attr}")


def _dataset_to_encoding(dataset: xr.Dataset, complevel: int) -> dict[Hashable, Any]:
    """Generate encoding information for NetCDF4 encoding from a dataset

    Parameters
    ----------
    dataset : xr.Dataset
        Dataset to generate encoding information for
    complevel : int
        The compression level to apply to arrays

    Returns
    -------
    dict[Hashable, Any]
        Resulting encoding settings
    """
    encoding = {
        var: {"compression": "gzip", "compression_opts": complevel} for var in dataset
    }
    return encoding


def write_shnitsel_file(
    dataset: xr.Dataset
    | xr.DataArray
    | SupportsToXrConversion
    | TreeNode[Any, xr.Dataset | xr.DataArray | SupportsToXrConversion],
    savepath: PathOptionsType,
    complevel: int = 9,
):
    """Function to write a trajectory in Shnitsel format (xr.) to a ntcdf hdf5 file format.

    Strips all internal attributes first to avoid errors during writing.
    When writing directly with to_netcdf, errors might occur due to internally set attributes with problematic types.

    Parameters
    ----------
    dataset : xr.Dataset | xr.DataArray | ShnitselDataset | SupportsToXrConversion | TreeNode[Any, xr.Dataset | xr.DataArray | SupportsToXrConversion]
        The dataset or trajectory to write (omit if using accessor).
    savepath : PathOptionsType
        The path at which to save the trajectory file.
    complevel : int, optional
        The compression level to apply during saving, by default 9

    Returns
    -------
    Unknown
        Returns the result of the final call to xr.Dataset.to_netcdf() or xr.DataTree.to_netcdf()
    """
    savepath_obj: pathlib.Path = make_uniform_path(savepath)  # type: ignore

    # Make sure the extension is appropriately set.
    if not savepath_obj.name.endswith(".nc"):
        savepath_obj = savepath_obj.parent / (savepath_obj.name + ".nc")

    if isinstance(dataset, TreeNode):
        tmp_res = tree_to_xarray_datatree(dataset)
        if tmp_res is None:
            raise ValueError("Tree could not be converted to netcdf conforming format.")

        cleaned_tree = _prepare_datatree(tmp_res)

        encoding = {}
        for leaf in cleaned_tree.leaves:
            if leaf.has_data:
                encoding[leaf.path] = _dataset_to_encoding(leaf.dataset, complevel)

        cleaned_tree.attrs["__shnitsel_format_version"] = "v1.3"
        # import pprint

        # pprint.pprint(cleaned_tree)
        return cleaned_tree.to_netcdf(savepath, engine='h5netcdf', encoding=encoding)
    elif isinstance(dataset, xr.Dataset):
        cleaned_ds = _prepare_dataset(dataset)
        encoding = _dataset_to_encoding(cleaned_ds, complevel)

        cleaned_ds.attrs["__shnitsel_format_version"] = "v1.3"

        return cleaned_ds.to_netcdf(savepath, engine='h5netcdf', encoding=encoding)
    else:
        ds, metadata = data_to_xarray_dataset(dataset, dict())
        if ds is None:
            raise ValueError("Data not be converted to netcdf conforming format.")
        ds.attrs["_shnitsel_io_meta"] = metadata

        return write_shnitsel_file(ds, savepath=savepath)
