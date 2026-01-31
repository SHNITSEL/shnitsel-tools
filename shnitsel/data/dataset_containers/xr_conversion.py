import logging
from typing import Any

from shnitsel.data.dataset_containers.shared import ShnitselDataset
from ..xr_io_compatibility import (
    SupportsToXrConversion,
    SupportsFromXrConversion,
)
from ..xr_io_registry import get_registered_input_handler
from . import Trajectory, Frames, InterState, PerState, wrap_dataset
import xarray as xr


def data_to_xarray_dataset(
    raw_data: (
        xr.DataArray | xr.Dataset | ShnitselDataset | None | SupportsToXrConversion
    ),
    metadata: dict[str, Any],
) -> tuple[xr.Dataset | None, dict[str, Any]]:
    """Support function to convert simple, non-hierarchical data structures (that could be in trees)
    to xarray datasets.

    Parameters
    ----------
    raw_data : xr.DataArray | xr.Dataset | Trajectory | Frames | InterState | PerState | None | SupportsToXrConversion
        Any type of raw data that should be converted. See type hints for supported types.
        All explicitly listed types will also be supported by deserialization.
        If the type is only supported through the `SupportsToXrConversion` interface, the type needs to be
        registered using the `xr_io_compatibility.register_custom_xr_input_type()` function before the conversion is attempted.
        Types not explicitly listed and not supporting the io conversion protocol will trigger an error.
    metadata : dict[str, Any]
        Metadata attributes to update according to the performed conversion.

    Returns
    -------
    xr.Dataset | None
        The converted dataset result or none if no conversion was possible
    dict[str, Any]
        The updated metadata dictionary

    Raises
    ------
    ValueError
        If an unsupported type is provided for conversion.
    """
    if raw_data is None:
        tree_data = None
    elif isinstance(raw_data, xr.DataArray):
        metadata["_shnitsel_data_type"] = "xarray::DataAray"
        metadata["_shnitsel_data_var"] = "data"
        tree_data = raw_data.to_dataset(name="data")
    elif isinstance(raw_data, xr.Dataset):
        metadata["_shnitsel_data_type"] = "xarray::Dataset"
        tree_data = raw_data
    # elif isinstance(raw_data, Frames):
    #     metadata["_shnitsel_data_type"] = "shnitsel::Frames"
    #     tree_data = raw_data.dataset
    # elif isinstance(raw_data, Trajectory):
    #     metadata["_shnitsel_data_type"] = "shnitsel::Trajectory"
    #     tree_data = raw_data.dataset
    # elif isinstance(raw_data, InterState):
    #     metadata["_shnitsel_data_type"] = "shnitsel::InterState"
    #     tree_data = raw_data.dataset
    # elif isinstance(raw_data, PerState):
    #     metadata["_shnitsel_data_type"] = "shnitsel::PerState"
    #     tree_data = raw_data.dataset
    elif isinstance(raw_data, SupportsToXrConversion):
        io_type_key, serialized_data, conversion_metadata = raw_data.as_xr_dataset()
        metadata["_shnitsel_data_type"] = (
            io_type_key if io_type_key is not None else "xarray::Dataset"
        )

        # Just a precaution
        if "_shnitsel_data_type_meta" in metadata:
            del metadata["_shnitsel_data_type_meta"]

        metadata["_shnitsel_data_type_meta"] = conversion_metadata

        tree_data = serialized_data

    else:
        logging.error(
            "Currently unsupported type %s found in data to be converted to xarray dataset.",
            type(raw_data),
        )
        raise ValueError(
            "Currently unsupported type %s found in data to be converted to xarray dataset."
            % type(raw_data)
        )
    return (tree_data, metadata)


def xr_dataset_to_shnitsel_format(
    raw_data: xr.Dataset,
    metadata: dict[str, Any],
) -> ShnitselDataset | xr.Dataset | xr.DataArray | SupportsFromXrConversion:
    """Function to support abstract deserialization of various types that could appear
    in a shnitsel style hierarchical tree structure.

    Will attempt to read type information from the `metadata` dict

    Parameters
    ----------
    raw_data : xr.Dataset
        The raw dataset input that should be deserialized.
    metadtata: dict[str, Any]
        Metadata dictionary for this entry. Should have a field `_shnitsel_data_type` with
        a type tag identifying the type that has been serialized in the provided `raw_data`
        parameter.
        For custom types, should also have a field `_shnitsel_data_type_meta` if metadata is
        required for deserialization but will be initialized to an empty dict if not provided.

    Returns
    -------
    Trajectory | Frames | InterState | PerState | xr.Dataset | xr.DataArray | SupportsFromXrConversion
        A deserialized object of a core shnitsel tools type, core xarray type or a type registered with
        `xr_io_compatibility.register_custom_xr_input_type()`.
        If deserialization failed, will be returned as dataset
    """

    if "_shnitsel_data_type" in metadata:
        shnitsel_type_hint = metadata["_shnitsel_data_type"]
        del metadata['_shnitsel_data_type']

        if shnitsel_type_hint == "xarray::Dataset":
            return raw_data
        elif shnitsel_type_hint == "xarray::DataAray":
            array_var_name = metadata.get("_shnitsel_data_var", "data")
            return raw_data[array_var_name]
        # elif shnitsel_type_hint == "shnitsel::Frames":
        #     return Frames(raw_data)
        # elif shnitsel_type_hint == "shnitsel::Trajectory":
        #     return Trajectory(raw_data)
        # elif shnitsel_type_hint == "shnitsel::InterState":
        #     return InterState(direct_interstate_data=raw_data)
        # elif shnitsel_type_hint == "shnitsel::PerState":
        #     return PerState(direct_perstate_data=raw_data)
        else:
            io_type_tag = shnitsel_type_hint
            io_metadata = metadata.get("_shnitsel_data_type_meta", {})

            if '_shnitsel_data_type_meta' in metadata:
                del metadata['_shnitsel_data_type_meta']

            input_handler = get_registered_input_handler(io_type_tag)
            if input_handler is None:
                logging.error(
                    "Unknown shnitsel deserialization type: %s. Returned as Dataset instead",
                    shnitsel_type_hint,
                )
                return raw_data
            else:
                try:
                    res = input_handler.from_xr_dataset(raw_data, io_metadata)
                except TypeError as e:
                    logging.error(
                        "Deserialization with input handler of type %s failed for data marked with tag %s. Will return as dataset. Error: %s.",
                        input_handler,
                        io_type_tag,
                        e,
                    )
                    return raw_data
                return res
    else:
        # Make best effort guess for the type:
        if (
            "energy_interstate" in raw_data
            or "dip_trans_norm" in raw_data
            or "fosc" in raw_data
            or "nacs_norm" in raw_data
            or "socs_norm" in raw_data
        ):
            try:
                return InterState(direct_interstate_data=raw_data)
            except AssertionError:
                pass

        if "forces_norm" in raw_data or "dip_perm_norm" in raw_data:
            try:
                return PerState(direct_perstate_data=raw_data)
            except AssertionError:
                pass
        # Try our best to wrap it in any type
        return wrap_dataset(raw_data)
