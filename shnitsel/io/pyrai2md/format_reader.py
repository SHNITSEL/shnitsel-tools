from dataclasses import dataclass
from glob import glob
import logging
import pathlib
from typing import Dict

from shnitsel.data.TrajectoryFormat import Trajectory
from ..FormatReader import FormatInformation, FormatReader
from .parse import parse_pyrai2md


@dataclass
class PyrAI2mdFormatInformation(FormatInformation):
    energy_file_path: pathlib.Path | None = None
    log_file_path: pathlib.Path | None = None
    pass


class PyrAI2mdFormatReader(FormatReader):
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

        is_request_specific_to_pyrai2md = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "pyrai2md"
        )

        md_energies_paths = glob(
            "*.md.energies",
            root_dir=path,
        )
        if (n := len(md_energies_paths)) != 1:
            message = (
                f"Path `{path}` does not constitute a PyrAI2md style output directory: Expected to find a single file ending with '.md.energies' "
                f"but found {n} files: {md_energies_paths}"
            )
            if is_request_specific_to_pyrai2md:
                logging.error(message)
            else:
                logging.debug(message)
            raise FileNotFoundError(message)

        energy_file_path = path / md_energies_paths[0]

        log_paths = glob(
            "*.log",
            root_dir=path,
        )
        if (n := len(md_energies_paths)) != 1:
            message = (
                "Path `{path}` does not constitute a PyrAI2md style output directory: Expected to find a single file ending with '.log' "
                f"but found {n} files: {log_paths}"
            )
            if is_request_specific_to_pyrai2md:
                logging.error(message)
            else:
                logging.debug(message)
            raise FileNotFoundError(message)

        log_file_path = path / log_paths[0]

        return PyrAI2mdFormatInformation(
            "SHARC", "unkown", path, energy_file_path, log_file_path
        )

    def read_from_path(
        self, path: pathlib.Path | None, format_info: FormatInformation | None = None
    ) -> Trajectory:
        """Read a PyrAI2md-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (pathlib.Path | None): Path to a PyrAI2md-format directory. If not provided explicitly, needs to be included in `format_info.path`
            format_info (FormatInformation | None, optional): Format information on the provided `path` if previously parsed. Will be parsed from `path` if not provided. Defaults to None.

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            FileNotFoundError: Path was not found or was not of appropriate PyrAI2md format

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
            loaded_dataset = parse_pyrai2md(path)
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading PyrAI2md trajectory from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return Trajectory(loaded_dataset)
