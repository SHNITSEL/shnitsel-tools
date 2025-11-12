from dataclasses import dataclass
import logging
import pathlib
import re
from typing import Dict, List, Tuple

from shnitsel.data.shnitsel_db_format import ShnitselDB
from shnitsel.data.trajectory_format import Trajectory
from shnitsel.io.helpers import LoadingParameters, PathOptionsType, make_uniform_path
from ..format_reader_base import FormatInformation, FormatReader
from .parse import read_shnitsel_file
from shnitsel.units.definitions import standard_shnitsel_units

import xarray as xr


@dataclass
class ShnitselFormatInformation(FormatInformation):
    pass


_shnitsel_default_pattern_regex = None
_shnitsel_default_pattern_glob = "*.nc"


class ShnitselFormatReader(FormatReader):
    """Class for providing the Shnitsel format reading functionality in the standardized `FormatReader` interface"""

    def find_candidates_in_directory(
        self, path: PathOptionsType
    ) -> List[pathlib.Path] | None:
        """Function to return a all potential matches for the current file format  within a provided directory at `path`.

        Returns:
            List[PathOptionsType] : A list of paths that should be checked in detail for whether they represent the format of this FormatReader.
            None: No potential candidate found
        """
        path_obj = make_uniform_path(path)
        res_entries = [
            e for e in path_obj.glob(_shnitsel_default_pattern_glob) if e.is_file()
        ]
        return None if len(res_entries) == 0 else res_entries

    def check_path_for_format_info(
        self, path: PathOptionsType, hints_or_settings: Dict | None = None
    ) -> FormatInformation:
        """Check if the `path` is a Shnitsel-style file

        Args:
            path (pathlib.Path): The path to check for shnitsel-style data
            hints_or_settings (Dict): Configuration options provided to the reader by the user

        Raises:
            FileNotFoundError: If the `path` is not a file.
            FileNotFoundError: If `path` is a file but not in the right format (i.e. not with `.nc` extension)

        Returns:
            FormatInformation: _description_
        """
        path_obj: pathlib.Path = make_uniform_path(path)

        _is_request_specific_to_shnitsel = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "shnitsel"
        )

        if not path_obj.exists() or not path_obj.is_file():
            message = f"Path `{path}` does not constitute a Shnitsel style trajectory file. Does not exist or is not a file."
            logging.debug(message)
            raise FileNotFoundError(message)

        if not path_obj.suffix.endswith(".nc"):
            message = f"Path `{path}` is not a NetCdf file (extension `.nc`)"

            logging.debug(message)
            raise FileNotFoundError(message)

        return ShnitselFormatInformation("shnitsel", "0.1", None, path_obj)

    def read_from_path(
        self,
        path: pathlib.Path,
        format_info: FormatInformation,
        loading_parameters: LoadingParameters | None = None,
    ) -> Trajectory | ShnitselDB:
        """Read a shnitsel-style file from `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (pathlib.Path): Path to a shnitsel-format `.nc` file.
            format_info (FormatInformation): Format information on the provided `path` that has been previously parsed.
            loading_parameters: (LoadingParameters|None, optional): Loading parameters to e.g. override default state names, units or configure the error reporting behavior

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            FileNotFoundError: Path was not found or was not of appropriate Shnitsel format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """
        try:
            loaded_dataset = read_shnitsel_file(
                path, loading_parameters=loading_parameters
            )
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading shnitsel file from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return loaded_dataset  # type: ignore # We know that the result of read_shnitsel_file is meant to be a ShnitselDB or single Trajectory

    def get_units_with_defaults(
        self, unit_overrides: Dict[str, str] | None = None
    ) -> Dict[str, str]:
        """Apply units to the default unit dictionary of the format SHNITSEL

        Args:
            unit_overrides (Dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            Dict[str, str]: The resulting, overridden default units
        """

        res_units = standard_shnitsel_units.copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
