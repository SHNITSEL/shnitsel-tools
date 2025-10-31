from dataclasses import dataclass
import logging
import pathlib
from typing import Dict, List

from shnitsel.data.TrajectoryFormat import Trajectory
from ..FormatReader import FormatInformation, FormatReader
from .parse_trajectory import read_traj
from .parse_initial_conditions import list_iconds, dir_of_iconds


@dataclass
class SHARCDynamicFormatInformation(FormatInformation):
    pass


@dataclass
class SHARCInitialFormatInformation(FormatInformation):
    list_of_iconds: List | None = None


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

        # Check if dynamic SHARC format satisfied
        is_dynamic = False
        format_information: FormatInformation | None = None
        try:
            input_file_path = path / "input"
            input_dat_path = path / "output.dat"
            input_xyz_path = path / "output.xyz"

            for file in [input_file_path, input_dat_path, input_xyz_path]:
                if not file.is_file():
                    message = f"Input directory `{path}` is missing {file}"
                    if is_request_specific_to_sharc:
                        logging.error(message)
                    else:
                        logging.debug(message)
                    raise FileNotFoundError(message)
            is_dynamic = True
            SHARCDynamicFormatInformation("SHARC", "unkown", path)
        except Exception as e:
            dynamic_check_error = e

        # Check if static/initial condition SHARC format satisfied

        is_static = False
        try:
            list_of_initial_condition_paths = list_iconds(path)
            is_static = True
            SHARCInitialFormatInformation(
                "SHARC", "unkown", path, list_of_initial_condition_paths
            )
        except Exception as e:
            static_check_error = e

        if is_dynamic and is_static:
            message = (
                f"Input directory {path} contains both static initial conditions and dynamic trajectory data of type SHARC."
                f"Please only point to a directory containing exactly one of the two kinds of data"
            )
            if is_request_specific_to_sharc:
                logging.error(message)
            else:
                logging.debug(message)
            raise ValueError(message)
        if format_information is None:
            message = (
                f"Input directory {path} contains neither static initial conditions nor dynamic trajectory data of type SHARC."
                f"Please point to a directory containing exactly one of the two kinds of data"
            )
            if is_request_specific_to_sharc:
                logging.error(message)
            else:
                logging.debug(message)
            raise FileNotFoundError(message)

        return format_information

    def read_from_path(
        self, path: pathlib.Path | None, format_info: FormatInformation | None = None
    ) -> Trajectory:
        """Read a SHARC-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (pathlib.Path | None): Path to a shnitsel-format `.nc` file. If not provided explicitly, needs to be included in `format_info.path`
            format_info (FormatInformation | None, optional): Format information on the provided `path` if previously parsed. Will be parsed from `path` if not provided. Defaults to None.

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            ValueError: `format_info` was of a wrong non-SHARC type.
            FileNotFoundError: Path was not found or was not of appropriate Shnitsel format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """

        is_dynamic = False

        if path is not None and format_info is None:
            format_info = self.check_path_for_format_info(path)
        elif path is None and format_info is not None:
            path = format_info.path
        elif path is None and format_info is None:
            raise ValueError("Either `path` or `format_info` needs to be provided")

        if isinstance(format_info, SHARCDynamicFormatInformation):
            is_dynamic = True
        elif isinstance(format_info, SHARCInitialFormatInformation):
            is_dynamic = False
        else:
            raise ValueError("The provided `format_info` object is not SHARC-specific.")

        if path is None:
            raise ValueError(
                "Not sufficient `path` information provided. Please set the `path` parameter"
            )

        try:
            if is_dynamic:
                loaded_dataset = read_traj(path)
            else:
                loaded_dataset = dir_of_iconds(path)
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading SHARC trajectory from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return Trajectory(loaded_dataset)
