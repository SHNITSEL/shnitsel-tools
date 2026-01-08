from abc import abstractmethod
import logging
from typing import Protocol, TypeAlias, TypeVar
import xarray as xr

ResType = TypeVar("ResType")
MetaData: TypeAlias = dict[str, str]


class SupportsToXrConversion(Protocol):
    """Definition of the protocol to support conversion of a type into
    xarray Dataset structs mostly for io purposes
    """

    @abstractmethod
    def as_xr_dataset(self) -> tuple[str | None, xr.Dataset, MetaData]:
        """Base function to implement by classes supporting this protocol
        to allow for standardized conversion to a dataset

        Returns
        -------
        tuple[str, xr.Dataset, MetaData]
            A tuple of the `io_type_tag` under which the deserializer is registered
            with the Shnitsel Tools framework (or `None` if no
            deserialization is desired) then the xr.Dataset that is the result of the conversion
            and last a dict of metadata that might help with deserialization later on.

        Raises
        ------
        ValueError
            If the conversion failed for some reason.
        """
        raise NotImplementedError(
            "The class %s did not implement the `as_xr_dataset` method." % type(self)
        )


class SupportsFromXrConversion(Protocol):
    """Definition of the protocol to support instantiation from
    xarray dataset structs.
    """

    @classmethod
    @abstractmethod
    def from_xr_dataset(
        cls: type[ResType], dataset: xr.Dataset, metadata: MetaData
    ) -> ResType:
        """Class method to support standardized deserialization of arbitrary classes.
        Implemented as a class method to avoid need to construct instance for
        deserialization.

        Parameters
        ----------
        cls : type[ResType]
            The class executing the deserialization.
        dataset : xr.Dataset
            The dataset to be deserialized into the output type.
        metadata : MetaData
            Metdatata from the serialization process.

        Returns
        -------
        instance of cls
            The deserialized instance of the target class.

        Raises
        ------
        TypeError
            If deserialization of the object was not possible
        """
        raise NotImplementedError(
            "The class %s did not implement the `from_xr_dataset` method." % cls
        )


Readable = TypeVar("Readable", bound=SupportsFromXrConversion)

INPUT_TYPE_REGISTRY: dict[str, type[SupportsFromXrConversion]] = {}


def register_custom_xr_input_type(cls: type[Readable], io_type_tag: str) -> bool:
    """Function to register a custom type that can be parsed from an xarray Dataset
    using the `SupportsFromXrConversion` protocol definition to support input from
    xarray style netcdf files.

    Parameters
    ----------
    cls : type[Readable]
        A class supporting the `SupportsFromXrConversion` protocol to be invoked to deserialize an object from
        a `xr.Dataset` instance when reading a netcdf file.
    io_type_tag : str
        The string type tag to be used to mark this class as the executing instance for
        deserialization.

    Returns
    -------
    bool
        True if registration succeeded. False if there was a clash with the `io_type_tag` of an existing type.
    """
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
