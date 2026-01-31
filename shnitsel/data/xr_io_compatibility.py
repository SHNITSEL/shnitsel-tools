from abc import abstractmethod
import logging
from typing import Protocol, TypeAlias, TypeVar, runtime_checkable
import xarray as xr

ResType = TypeVar("ResType")
MetaData: TypeAlias = dict[str, str]


@runtime_checkable
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
            deserialization is desired/supported)/
            Then the `xr.Dataset that is the result of the conversion.
            And lastly a dict of metadata that might help with deserialization later on.

        Raises
        ------
        ValueError
            If the conversion failed for some reason.
        """
        raise NotImplementedError(
            "The class %s did not implement the `as_xr_dataset` method." % type(self)
        )


@runtime_checkable
class SupportsFromXrConversion(Protocol):
    """Definition of the protocol to support instantiation from
    xarray dataset structs.
    """

    @classmethod
    @abstractmethod
    def get_type_marker(cls) -> str:
        raise NotImplementedError(
            "The class %s did not implement the `get_type_marker()` method." % cls
        )

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
