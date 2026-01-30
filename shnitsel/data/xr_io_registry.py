import logging
from typing import TypeVar

from .xr_io_compatibility import SupportsFromXrConversion


Readable = TypeVar("Readable", bound=SupportsFromXrConversion)

INPUT_TYPE_REGISTRY: dict[str, type[SupportsFromXrConversion]] = {}


def register_custom_xr_input_type(
    cls: type[Readable], io_type_tag: str | None = None
) -> bool:
    """Function to register a custom type that can be parsed from an xarray Dataset
    using the `SupportsFromXrConversion` protocol definition to support input from
    xarray style netcdf files.

    Parameters
    ----------
    cls : type[Readable]
        A class supporting the `SupportsFromXrConversion` protocol to be invoked to deserialize an object from
        a `xr.Dataset` instance when reading a netcdf file.
    io_type_tag : str, optional
        The string type tag to be used to mark this class as the executing instance for
        deserialization.
        If not set, will use `get_type_marker()` on the registered class.

    Returns
    -------
    bool
        True if registration succeeded. False if there was a clash with the `io_type_tag` of an existing type.
    """
    if io_type_tag is None:
        io_type_tag = cls.get_type_marker()

    if io_type_tag in INPUT_TYPE_REGISTRY:
        logging.error("IO type tag already in use: %s", io_type_tag)
        return False
    else:
        INPUT_TYPE_REGISTRY[io_type_tag] = cls
        return True


def get_registered_input_handler(
    io_type_tag: str,
) -> type[SupportsFromXrConversion] | None:
    """Function to look up a potentially registered input handler for a previously serialized data type.

    If no input handler is registered, will return `None`

    Parameters
    ----------
    io_type_tag : str
        The type tag under which the type was registered in the system with `register_custom_xr_input_type()`.

    Returns
    -------
    type[SupportsFromXrConversion]
        Either a class object that supports the protocol of an input handler or `None` if no handler was found.
    """
    if io_type_tag in INPUT_TYPE_REGISTRY:
        return INPUT_TYPE_REGISTRY[io_type_tag]
    else:
        return None


def setup_defaults():
    """Helper function to register all default shnitsel data types with the xr io compatibility registry."""
    from .dataset_containers import (
        ShnitselDataset,
        DataSeries,
        Trajectory,
        Frames,
        PerState,
        InterState,
        MultiSeriesDataset,
        MultiSeriesLayered,
        MultiSeriesStacked,
    )

    # Only register if they are not there yet
    if not INPUT_TYPE_REGISTRY:
        register_custom_xr_input_type(ShnitselDataset)
        register_custom_xr_input_type(DataSeries)
        register_custom_xr_input_type(Trajectory)
        register_custom_xr_input_type(Frames)
        register_custom_xr_input_type(PerState)
        register_custom_xr_input_type(InterState)
        register_custom_xr_input_type(MultiSeriesDataset)
        register_custom_xr_input_type(MultiSeriesLayered)
        register_custom_xr_input_type(MultiSeriesStacked)


setup_defaults()
