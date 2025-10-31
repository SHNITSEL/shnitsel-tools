from dataclasses import dataclass
import logging
import pathlib
from typing import Dict

from shnitsel.data.TrajectoryFormat import Trajectory
from ..FormatReader import FormatInformation, FormatReader
from .trajectory import read_traj


@dataclass
class SHARCFormatInformation(FormatInformation):
    pass


class SHARCFormatReader(FormatReader):
    """Class for providing the SHARC format reading functionality in the standardized `FormatReader` interface"""

    def check_path_for_format_info(
        self, path: pathlib.Path, hints_or_settings: Dict | None = None
    ) -> FormatInformation:
        """Check if the `path` is a SHARC-style output directory.

        Designed for a single input trajectory.

        Args:
            path (pathlib.Path): The path to check for SHARC data
            hints_or_settings (Dict): Configuration options provided to the reader by the user

        Raises:
            FileNotFoundError: If the `path` is not a directory.
            FileNotFoundError: If `path` is a directory but does not contain the required SHARC output files

        Returns:
            FormatInformation: _description_
        """

        is_request_specific_to_sharc = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "sharc"
        )

        if not path.exists() or not path.is_dir():
            message = f"Path `{path}` does not constitute a SHARC style output directory: Does not exist or is not a directory."
            if is_request_specific_to_sharc:
                logging.error(message)
            else:
                logging.debug(message)
            raise FileNotFoundError(message)
        
        input_file_path = path / "input"
        input_dat_path = path / "output.dat"
        input_xyz_path = path / "output.xyz"

        for file in [input_file_path, input_dat_path, input_xyz_path]:
            if not file.is_file():
                message = f"Input directory is missing {file}"
                if is_request_specific_to_sharc:
                    logging.error(message)
                else:
                    logging.debug(message)
                raise FileNotFoundError(message)

        return SHARCFormatInformation("SHARC", "unkown", path)

    def read_from_path(
        self, path: pathlib.Path | None, format_info: FormatInformation | None = None
    ) -> Trajectory:
        """Read a SHARC-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

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
            raise ValueError("Either `path` or `format_info` needs to be provided")

        if path is None:
            raise ValueError(
                "Not sufficient `path` information provided. Please set the `path` parameter"
            )

        try:
            loaded_dataset = read_traj(path)
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading SHARC trajectory from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return Trajectory(loaded_dataset)
