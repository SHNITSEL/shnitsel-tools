import logging
from typing import Literal, Mapping, TypeAlias

from shnitsel.io.ase.format_reader import ASEFormatReader

from .format_reader_base import FormatReader
from .pyrai2md.format_reader import PyrAI2mdFormatReader
from .sharc.format_reader import SHARCFormatReader
from .shnitsel.format_reader import ShnitselFormatReader
from .newtonx.format_reader import NewtonXFormatReader


# TODO: FIXME: add ASE support
_newton_reader = NewtonXFormatReader()
REGISTERED_READERS: dict[str, FormatReader] = {
    "nx": _newton_reader,  # parse_newtonx,
    "newtonx": _newton_reader,  # parse_newtonx,
    "sharc": SHARCFormatReader(),  # parse_sharc,
    "pyrai2md": PyrAI2mdFormatReader(),
    "shnitsel": ShnitselFormatReader(),  # read_shnitsel_file,
}

FormatIdentifierType: TypeAlias = (
    Literal['nx', 'newtonx', 'sharc', 'pyrai2md', 'shnitsel'] | str
)


def register_format_reader(
    format_identifier: str, reader_instance: FormatReader
) -> bool:
    """Function to register a new input format implementation with the io registry.

    Parameters
    ----------
    format_identifier : str
        The string with which the format should be identified. E.g. 'newtonx' or 'sharc'. Must not already be in use. Identifiers already in use: 'nx', 'newtonx', 'sharc', 'pyrai2md', 'shnitsel'
    reader_instance : FormatReader
        An instance of an object of a format-specific subclass of `FormatReader` to handle io requests.

    Returns
    -------
    bool
        Returns True if registration was successful. False if there was a clash with an existing format declaration.
    """

    normalized_identifier = format_identifier.lower()

    if normalized_identifier in REGISTERED_READERS:
        logging.error(
            "Format %s already mapped to an io class. Cannot reassign.",
            format_identifier,
        )
        return False

    if reader_instance is None or not isinstance(reader_instance, FormatReader):
        logging.error(
            "Format reader is not of required type %s. Was of type %s. Cannot register format.",
            FormatReader,
            type(reader_instance),
        )
        return False

    REGISTERED_READERS[normalized_identifier] = reader_instance
    return True


def get_available_io_handlers() -> Mapping[str, FormatReader]:
    """Function to access the list of available input reader classes.

    Returns
    -------
    Mapping[str, FormatReader]
        The mapping from format identifier to actual io handler for the various supported formats.
    """
    return REGISTERED_READERS


# Register additional format readers
register_format_reader('ase', ASEFormatReader())
