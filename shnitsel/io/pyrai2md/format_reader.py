from dataclasses import dataclass
from glob import glob
import logging
import pathlib
import re
import sys
import traceback
from typing import Dict, List, Tuple

from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
from ..format_reader_base import FormatInformation, FormatReader
from .parse import parse_pyrai2md

import xarray as xr


@dataclass
class PyrAI2mdFormatInformation(FormatInformation):
    energy_file_path: pathlib.Path | None = None
    log_file_path: pathlib.Path | None = None
    pass


class PyrAI2mdFormatReader(FormatReader):
    """Class for providing the PyrAI2md format reading functionality in the standardized `FormatReader` interface"""

    def find_candidates_in_directory(
        self, path: PathOptionsType
    ) -> List[pathlib.Path] | None:
        """Function to return a all potential matches for the current file format  within a provided directory at `path`.

        Returns:
            List[PathOptionsType] : A list of paths that should be checked in detail for whether they represent the format of this FormatReader.
            None: No potential candidate found
        """
        # TODO: FIXME: Add option to specify if we want only file or only directory paths
        # TODO: FIXME: maybe just turn into a "filter" function and provide the paths?
        path_obj = make_uniform_path(path)

        res_entries = [e for e in path_obj.glob("*") if e.is_dir()]
        return None if len(res_entries) == 0 else res_entries

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
        path_obj: pathlib.Path = make_uniform_path(path)  # type: ignore
        base_format_info = super().check_path_for_format_info(
            path_obj, hints_or_settings
        )

        _is_request_specific_to_pyrai2md = (
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
            logging.debug(message)
            raise FileNotFoundError(message)

        log_file_path = path_obj / log_paths[0]

        return PyrAI2mdFormatInformation(
            "pyrai2md",
            "unkown",
            base_format_info.trajid,
            path_obj,
            energy_file_path,
            log_file_path,
        )

    def read_from_path(
        self,
        path: pathlib.Path,
        format_info: FormatInformation,
        loading_parameters: LoadingParameters | None = None,
    ) -> xr.Dataset:
        """Read a PyrAI2md-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (pathlib.Path): Path to a PyrAI2md-format directory.
            format_info (FormatInformation): Format information on the provided `path` that has been previously parsed.
            loading_parameters: (LoadingParameters|None, optional): Loading parameters to e.g. override default state names, units or configure the error reporting behavior

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            FileNotFoundError: Path was not found or was not of appropriate PyrAI2md format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """

        try:
            loaded_dataset = parse_pyrai2md(path, loading_parameters=loading_parameters)
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading PyrAI2md trajectory from path `{path}` failed because of original error: {v_e}.\n Trace: \n {traceback.format_exc()}"
            logging.error(message)
            raise FileNotFoundError(message)

        return loaded_dataset

    def get_units_with_defaults(
        self, unit_overrides: Dict[str, str] | None = None
    ) -> Dict[str, str]:
        """Apply units to the default unit dictionary of the format PyrAI2md

        Args:
            unit_overrides (Dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            Dict[str, str]: The resulting, overridden default units
        """
        from shnitsel.units.definitions import standard_units_of_formats

        res_units = standard_units_of_formats["pyrai2md"].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
