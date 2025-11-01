from dataclasses import dataclass
from glob import glob
import logging
import pathlib
import re
from typing import Dict, Tuple

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.helpers import LoadingParameters, PathOptionsType, make_uniform_path
from ..FormatReader import FormatInformation, FormatReader
from .parse import parse_pyrai2md


@dataclass
class PyrAI2mdFormatInformation(FormatInformation):
    energy_file_path: pathlib.Path | None = None
    log_file_path: pathlib.Path | None = None
    pass


class PyrAI2mdFormatReader(FormatReader):
    """Class for providing the SHARC format reading functionality in the standardized `FormatReader` interface"""

    def get_default_trajectory_pattern(self) -> Tuple[str, re.Pattern | None] | None:
        """Function to retrieve PyrAI2md specific naming convention for trajectory directories.

        Defaults to checking all subdirectories

        Returns:
            Tuple[str, re.Pattern | None] | None: _description_
        """
        return ("*", None)

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
        path_obj: pathlib.Path = make_uniform_path(path)
        is_request_specific_to_pyrai2md = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "pyrai2md"
        )

        md_energies_paths = glob(
            "*.md.energies",
            root_dir=path_obj,
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

        energy_file_path = path_obj / md_energies_paths[0]

        log_paths = glob(
            "*.log",
            root_dir=path_obj,
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

        log_file_path = path_obj / log_paths[0]

        return PyrAI2mdFormatInformation(
            "sharc", "unkown", path_obj, energy_file_path, log_file_path
        )

    def read_from_path(
        self, path: PathOptionsType | None,
        format_info: FormatInformation | None = None,
        loading_parameters: LoadingParameters | None = None
    ) -> Trajectory:
        """Read a PyrAI2md-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (PathOptionsType | None): Path to a PyrAI2md-format directory. If not provided explicitly, needs to be included in `format_info.path`
            format_info (FormatInformation | None, optional): Format information on the provided `path` if previously parsed. Will be parsed from `path` if not provided. Defaults to None.
            loading_parameters: (LoadingParameters|None, optional): Loading parameters to e.g. override default state names, units or configure the error reporting behavior

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            FileNotFoundError: Path was not found or was not of appropriate PyrAI2md format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """

        path_obj: pathlib.Path = make_uniform_path(path)

        if path_obj is not None and format_info is None:
            format_info = self.check_path_for_format_info(path_obj)
        elif path_obj is None and format_info is not None:
            path_obj = format_info.path
        elif path_obj is None and format_info is None:
            raise ValueError(
                "Either `path` or `format_info` needs to be provided")

        if path_obj is None:
            raise ValueError(
                "Not sufficient `path` information provided. Please set the `path` parameter"
            )

        try:
            loaded_dataset = parse_pyrai2md(
                path_obj, loading_parameters=self.get_loading_parameters_with_defaults(loading_parameters))
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading PyrAI2md trajectory from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return Trajectory(loaded_dataset)

    def get_units_with_defaults(self, unit_overrides: Dict[str, str] | None = None) -> Dict[str, str]:
        """Apply units to the default unit dictionary of the format PyrAI2md

        Args:
            unit_overrides (Dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            Dict[str, str]: The resulting, overridden default units
        """
        from shnitsel.units.definitions import standard_units_of_formats

        res_units = standard_units_of_formats['pyrai2md'].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)
            
        return res_units