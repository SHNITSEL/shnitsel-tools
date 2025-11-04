from dataclasses import dataclass
from glob import glob
import logging
import os
import pathlib
import re
from typing import Dict, List, Tuple

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.helpers import LoadingParameters, PathOptionsType, make_uniform_path
from ..FormatReader import FormatInformation, FormatReader
from .parse import parse_newtonx


@dataclass
class NewtonXFormatInformation(FormatInformation):
    nx_log_path: pathlib.Path | None = None
    positions_file_path: pathlib.Path | None = None
    pass


_newtonx_default_pattern_regex = re.compile(r"TRAJ(?P<trajid>\d+)")
_newtonx_default_pattern_glob = r"TRAJ*"


class NewtonXFormatReader(FormatReader):
    """Class for providing the SHARC format reading functionality in the standardized `FormatReader` interface"""

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

        res_entries = [
            e
            for e in path_obj.glob(_newtonx_default_pattern_glob)
            if _newtonx_default_pattern_regex.match(e.name) and e.is_dir()
        ]
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
        path_obj: pathlib.Path = make_uniform_path(path)

        is_request_specific_to_NewtonX = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and (
                hints_or_settings["kind"] == "nx"
                or hints_or_settings["kind"] == "newtonx"
            )
        )

        nx_log_path = path_obj / "RESULTS" / "nx.log"
        nx_positions_path = path_obj / "RESULTS" / "dyn.xyz"
        for file in [nx_log_path, nx_positions_path]:
            if not file.is_file():
                message = f"Input directory is missing {file}"
                if is_request_specific_to_NewtonX:
                    logging.error(message)
                else:
                    logging.debug(message)
                raise FileNotFoundError(message)

        format_information = NewtonXFormatInformation(
            "newtonx", "unkown", path_obj, nx_log_path, nx_positions_path
        )

        # Try and extract a trajectory ID from the path name
        match_attempt = _newtonx_default_pattern_regex.match(path_obj.name)

        if match_attempt:
            path_based_trajid = match_attempt.group("trajid")
            format_information.trajid = int(path_based_trajid)

        return format_information

    def read_from_path(
        self,
        path: PathOptionsType | None,
        format_info: FormatInformation | None = None,
        loading_parameters: LoadingParameters | None = None,
    ) -> Trajectory:
        """Read a NewtonX-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (PathOptionsType | None): Path to a NewtonX-format directory. If not provided explicitly, needs to be included in `format_info.path`
            format_info (FormatInformation | None, optional): Format information on the provided `path` if previously parsed. Will be parsed from `path` if not provided. Defaults to None.
            loading_parameters: (LoadingParameters|None, optional): Loading parameters to e.g. override default state names, units or configure the error reporting behavior

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            FileNotFoundError: Path was not found or was not of appropriate NewtonX format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """

        path_obj: pathlib.Path = make_uniform_path(path)

        if path_obj is not None and format_info is None:
            format_info = self.check_path_for_format_info(path_obj)
        elif path_obj is None and format_info is not None:
            path_obj = format_info.path
        elif path_obj is None and format_info is None:
            raise ValueError("Either `path` or `format_info` needs to be provided")

        if path_obj is None:
            raise ValueError(
                "Not sufficient `path` information provided. Please set the `path` parameter"
            )

        try:
            loaded_dataset = parse_newtonx(
                path_obj,
                loading_parameters=self.get_loading_parameters_with_defaults(
                    loading_parameters
                ),
            )
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading NewtonX trajectory from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        # If trajid has been extracted from the input path, set it
        if format_info.trajid is not None:
            loaded_dataset.attrs["trajid"] = format_info.trajid
            logging.info(f"Assigning id {format_info.trajid} to trajectory")

        loaded_dataset.attrs["trajectory_input_path"] = format_info.path.as_posix()

        return Trajectory(loaded_dataset)

    def get_units_with_defaults(
        self, unit_overrides: Dict[str, str] | None = None
    ) -> Dict[str, str]:
        """Apply units to the default unit dictionary of the format NewtonX

        Args:
            unit_overrides (Dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            Dict[str, str]: The resulting, overridden default units
        """
        from shnitsel.units.definitions import standard_units_of_formats

        res_units = standard_units_of_formats["newtonx"].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
