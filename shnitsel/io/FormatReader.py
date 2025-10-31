from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
import pathlib
from typing import Dict

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.helpers import PathOptionsType


@dataclass
class FormatInformation:
    """Information to keep track of relevant information for"""

    format_name: str = "none"
    version: str = "none"
    path: pathlib.Path | None = None


class FormatReader(ABC):
    """Abstract base class for all input formats to define a unified input reader interface.

    Should be subclassed and the functions `check_path_for_format_info()` and `read_from_path()`
    overridden in the subclass
    """

    @abstractmethod
    def check_path_for_format_info(
        self, path: PathOptionsType, hints_or_settings: Dict | None = None
    ) -> FormatInformation:
        """Checks if a path is of a given format and returns a struct containing all relevant info for reading
        the format at this location. Additionally checks configured user settings provided in `hints_or_settings` whether they are
        consistent with the file format.

        Needs to be overridden by each format.

        Args:
            path (os.PathLike):
                The path to look for data from the respective method for.
                Depending on the format, this would need to point to a file or a directory containing the actual
                trajectory information
            hints (Dict|None, optional):
                Potential hints/configuration options provided by the user as input to the reader which can be
                checked for conflicts with the requirements of the format (i.e. requesting a static initial condition from a dynamic trajectory in SHARC).
                Defaults to None

        Raises:
            FileNotFoundError: If required files were not found, i.e. if the path does not actually constitute input data of the denoted format
            ValueError: If the hints/settings provided by the user conflict with the requirements of the format

        Returns:
            FormatInformation:
                A structure containing all of the information relevant to the interpretation or reading of the format.
                Can be used to differentiate different versions of the same format.
                Should be passed to the `read_from_path()` method of the same class.
        """
        return FormatInformation()

    @abstractmethod
    def read_from_path(
        self, path: PathOptionsType | None, format_info: FormatInformation | None = None
    ) -> Trajectory:
        """Method to read a path of the respective format (e.g. ) into a shnitsel-conform trajectory.

        The return value of type `Trajectory` is a wrapper for the raw `xarray.Dataset` read from the `path`.
        This allows provision of extra features like keeping track of the original data while post-processing is performed.

        Args:
            path (os.PathLike): Path to either the input file or input folder to be read.
            format_info (FormatInformation | None, optional): Format information previously constructed by `check_path_for_format_info()`. If None, will be constructed by calling `Self.check_path_for_format_info()` first. Defaults to None.

        Raises:
            FileNotFoundError: If required files were not found, i.e. if the path does not actually constitute input data of the denoted format
            ValueError: If the `format_info` provided by the user conflicts with the requirements of the format
            Valueerror: If neither `path` nor `format_info` are provided


        Returns:
            Trajectory: The parsed dataset as wrapper around `xarray.Dataset` to keep track of original data.
        """
        ...
