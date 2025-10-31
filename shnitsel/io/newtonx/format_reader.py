from dataclasses import dataclass
from glob import glob
import logging
import os
import pathlib
from typing import Dict

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.helpers import PathOptionsType, make_uniform_path
from ..FormatReader import FormatInformation, FormatReader
from .parse import parse_newtonx


@dataclass
class NewtonXFormatInformation(FormatInformation):
    nx_log_path: pathlib.Path | None = None
    positions_file_path: pathlib.Path | None = None
    pass


class NewtonXFormatReader(FormatReader):
    """Class for providing the SHARC format reading functionality in the standardized `FormatReader` interface"""

    def check_path_for_format_info(
        self, path: PathOptionsType, hints_or_settings: Dict | None = None
    ) -> FormatInformation:
        """Check if the `path` is a SHARC-style output directory.

        Designed for a single input trajectory.

        Args:
            path (PathOptionsType): The path to check for SHARC data
            hints_or_settings (Dict): Configuration options provided to the reader by the user

        Raises:
            FileNotFoundError: If the `path` is not a directory.
            FileNotFoundError: If `path` is a directory but does not contain the required SHARC output files

        Returns:
            FormatInformation: _description_
        """
        path: pathlib.Path = make_uniform_path(path)

        is_request_specific_to_NewtonX = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and (
                hints_or_settings["kind"] == "nx"
                or hints_or_settings["kind"] == "newtonx"
            )
        )

        nx_log_path = path / "RESULTS" / "nx.log"
        nx_positions_path = path / "RESULTS" / "dyn.xyz"
        for file in [nx_log_path, nx_positions_path]:
            if not file.is_file():
                message = f"Input directory is missing {file}"
                if is_request_specific_to_NewtonX:
                    logging.error(message)
                else:
                    logging.debug(message)
                raise FileNotFoundError(message)

        return NewtonXFormatInformation(
            "NewtonX", "unkown", path, nx_log_path, nx_positions_path
        )

    def read_from_path(
        self, path: PathOptionsType | None, format_info: FormatInformation | None = None
    ) -> Trajectory:
        """Read a NewtonX-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (PathOptionsType | None): Path to a NewtonX-format directory. If not provided explicitly, needs to be included in `format_info.path`
            format_info (FormatInformation | None, optional): Format information on the provided `path` if previously parsed. Will be parsed from `path` if not provided. Defaults to None.

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            FileNotFoundError: Path was not found or was not of appropriate NewtonX format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """

        path: pathlib.Path = make_uniform_path(path)

        if path is not None and format_info is None:
            format_info = self.check_path_for_format_info(path)
        elif path is None and format_info is not None:
            path = format_info.path
        elif path is None and format_info is None:
            raise ValueError("Either `path` or `format_info` needs to be provided")

        if path is None:
            raise ValueError(
                "Not sufficient `path` information provided. Please set the `path` parameter"
            )

        try:
            loaded_dataset = parse_newtonx(path)
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading NewtonX trajectory from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return Trajectory(loaded_dataset)
