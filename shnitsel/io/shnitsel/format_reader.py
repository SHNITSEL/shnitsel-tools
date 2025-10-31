from dataclasses import dataclass
import logging
import pathlib
from typing import Dict

from shnitsel.data.TrajectoryFormat import Trajectory
from ..FormatReader import FormatInformation, FormatReader
from .parse import read_shnitsel_file


@dataclass
class ShnitselFormatInformation(FormatInformation):
    pass


class ShnitselFormatReader(FormatReader):
    """Class for providing the Shnitsel format reading functionality in the standardized `FormatReader` interface"""

    def check_path_for_format_info(
        self, path: pathlib.Path, hints_or_settings: Dict | None = None
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

        is_request_specific_to_shnitsel = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "shnitsel"
        )

        if not path.exists() or not path.is_file():
            message = f"Path `{path}` does not constitute a Shnitsel style trajectory file. Does not exist or is not a file."
            if is_request_specific_to_shnitsel:
                logging.error(message)
            else:
                logging.debug(message)
            raise FileNotFoundError(message)

        if not path.suffix.endswith(".nc"):
            message = f"Path `{path}` is not a NetCdf file (extension `.nc`)"
            if is_request_specific_to_shnitsel:
                logging.error(message)
            else:
                logging.debug(message)
            raise FileNotFoundError(message)

        return ShnitselFormatInformation("shnitsel", "0.1", path)

    def read_from_path(
        self, path: pathlib.Path | None, format_info: FormatInformation | None = None
    ) -> Trajectory:
        """Read a shnitsel-style file from `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (pathlib.Path | None): Path to a shnitsel-format `.nc` file. If not provided explicitly, needs to be included in `format_info.path`
            format_info (FormatInformation | None, optional): Format information on the provided `path` if previously parsed. Will be parsed from `path` if not provided. Defaults to None.

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            FileNotFoundError: Path was not found or was not of appropriate Shnitsel format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """

        if path is not None and format_info is None:
            format_info = self.check_path_for_format_info(path)
        elif path is None and format_info is not None:
            path = format_info.path
        elif path is None and format_info is None:
            raise ValueError("Either `path` or `format_info` needs to be provided to ")

        if path is None:
            raise ValueError(
                "Not sufficient `path` information provided. Please set the `path` parameter"
            )

        try:
            loaded_dataset = read_shnitsel_file(path)
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading shnitsel file from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return Trajectory(loaded_dataset)
