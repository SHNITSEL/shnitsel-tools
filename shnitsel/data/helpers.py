import dataclasses
import logging
from typing import Dict, List, TypeVar

T = TypeVar("T")

def dataclass_from_dict(datatype: type[T], d: List | Dict | T) -> T:
    """Helper function to restore a Dataclass object from its dict representation.

    Mainly used for serialization or storage of data in the DataTree db structure.

    Args:
        datatype (Type[T]): The dataclass type to restore
        d (List|Dict|T): The datasource to convert back into the Dataclass instance.

    Raises:
        ValueError: If value decoding fails during reconstruction
        TypeError: If some type mismatch occurs between the provided dict and the target type

    Returns:
        T: A resonstructed `datatype` instance.
    """
    if isinstance(d, list):
        (inner,) = datatype.__args__ # type: ignore
        return [dataclass_from_dict(inner, i) for i in d] # type: ignore

    try:
        fieldtypes = {f.name: f.type for f in dataclasses.fields(datatype)} # type: ignore
        return datatype(**{f: dataclass_from_dict(fieldtypes[f], d[f]) for f in d}) # type: ignore
    except Exception as e:
        logging.exception(
            f"Failed to convert object {d} of type {type(d)} to class {datatype}: {e}"
        )
        return d  # type: ignore # Not a dataclass field
