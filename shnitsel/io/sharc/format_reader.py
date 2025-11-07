from dataclasses import dataclass
import logging
import pathlib
import re
from typing import Dict, List, Tuple

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.helpers import LoadingParameters, PathOptionsType, make_uniform_path
from ..format_reader_base import FormatInformation, FormatReader
from .parse_trajectory import read_traj
from .parse_initial_conditions import list_iconds, read_iconds_individual

import xarray as xr


@dataclass
class SHARCDynamicFormatInformation(FormatInformation):
    pass


@dataclass
class SHARCInitialFormatInformation(FormatInformation):
    pass


@dataclass
class SHARCMultiInitialFormatInformation(FormatInformation):
    list_of_iconds: List | None = None


_sharc_default_pattern_regex = re.compile(r"(?P<dynstat>TRAJ|ICOND)_(?P<trajid>\d+)")
_sharc_default_pattern_glob_traj = "TRAJ_*"
_sharc_default_pattern_glob_icond = "ICOND_*"


class SHARCFormatReader(FormatReader):
    """Class for providing the SHARC format reading functionality in the standardized `FormatReader` interface"""

    def find_candidates_in_directory(
        self, path: PathOptionsType
    ) -> List[pathlib.Path] | None:
        """Function to return a all potential matches for the current file format  within a provided directory at `path`.

        Returns:
            List[PathOptionsType] : A list of paths that should be checked in detail for whether they represent the format of this FormatReader.
            None: No potential candidate found
        """
        path_obj = make_uniform_path(path)

        tmp_entries_traj = [e for e in path_obj.glob(_sharc_default_pattern_glob_traj)]
        tmp_entries_icond = [
            e for e in path_obj.glob(_sharc_default_pattern_glob_icond)
        ]
        res_entries = [
            e
            for e in tmp_entries_traj + tmp_entries_icond
            if _sharc_default_pattern_regex.match(e.name) and e.is_dir()
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

        is_request_specific_to_sharc = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "sharc"
        )

        if not path_obj.exists() or not path_obj.is_dir():
            message = f"Path `{path}` does not constitute a SHARC style output directory: Does not exist or is not a directory."
            if is_request_specific_to_sharc:
                logging.error(message)
            else:
                logging.debug(message)
            raise FileNotFoundError(message)

        dontanalyze_file_path = path_obj / "DONT_ANALYZE"

        if dontanalyze_file_path.exists() and dontanalyze_file_path.is_file():
            message = f"The path {path} does contain a `DONT_ANALYZE` file and will therefore be skipped. Please remove that file if you want the directory to be read."
            logging.warning(message)
            raise FileNotFoundError(message)

        # Check if dynamic SHARC format satisfied
        is_dynamic = False
        format_information: FormatInformation | None = None
        try:
            input_file_path = path_obj / "input"
            input_dat_path = path_obj / "output.dat"
            input_xyz_path = path_obj / "output.xyz"

            for file in [input_file_path, input_dat_path, input_xyz_path]:
                if not file.is_file():
                    message = f"Input directory `{path}` is missing {file}"
                    if is_request_specific_to_sharc:
                        logging.error(message)
                    else:
                        logging.debug(message)
                    raise FileNotFoundError(message)
            is_dynamic = True
            format_information = SHARCDynamicFormatInformation(
                "sharc", "unkown", None, path_obj
            )
            logging.debug(
                f"Input directory `{path}` fulfils data requirements of dynamic SHARC trajectory"
            )
        except Exception as e:
            dynamic_check_error = e

        # Check if static/initial condition SHARC format satisfied

        is_static = False
        try:
            qm_out_path = path_obj / "QM.out"
            qm_log_path = path_obj / "QM.log"
            qm_in_path = path_obj / "QM.in"

            if not qm_out_path.is_file() or (
                not qm_log_path.is_file() and not qm_in_path.is_file()
            ):
                message = f"Input directory `{path}` is missing `QM.out` or both `QM.log` and `QM.in`"
                if is_request_specific_to_sharc:
                    logging.info(message)
                else:
                    logging.debug(message)
                raise FileNotFoundError(message)

            # list_of_initial_condition_paths = list_iconds(path_obj)
            is_static = True
            format_information = SHARCInitialFormatInformation(
                "sharc",
                "unkown",
                None,
                path_obj,  # , list_of_initial_condition_paths
            )
            logging.debug(
                f"Input directory `{path}` fulfils data requirements of SHARC Initial Conditions"
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

        # Try and extract a trajectory ID from the path name
        match_attempt = _sharc_default_pattern_regex.match(path.name)

        if match_attempt:
            path_based_trajid = match_attempt.group("trajid")
            format_information.trajid = int(path_based_trajid)

        return format_information

    def read_from_path(
        self,
        path: pathlib.Path,
        format_info: FormatInformation,
        loading_parameters: LoadingParameters | None = None,
    ) -> xr.Dataset:
        """Read a SHARC-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Args:
            path (pathlib.Path): Path to a SHARC-format directory.
            format_info (FormatInformation): Format information on the provided `path` that has been previously parsed.
            loading_parameters: (LoadingParameters|None, optional): Loading parameters to e.g. override default state names, units or configure the error reporting behavior

        Raises:
            ValueError: Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
            ValueError: `format_info` was of a wrong non-SHARC type.
            FileNotFoundError: Path was not found or was not of appropriate Shnitsel format

        Returns:
            Trajectory: The loaded Shnitsel-conforming trajectory
        """

        is_dynamic = False
        if isinstance(format_info, SHARCDynamicFormatInformation):
            is_dynamic = True
        elif isinstance(format_info, SHARCInitialFormatInformation):
            is_dynamic = False
        else:
            raise ValueError("The provided `format_info` object is not SHARC-specific.")

        try:
            if is_dynamic:
                loaded_dataset = read_traj(
                    path,
                    loading_parameters=self.get_loading_parameters_with_defaults(
                        loading_parameters
                    ),
                )
            else:
                loaded_dataset = read_iconds_individual(
                    path,
                    loading_parameters=self.get_loading_parameters_with_defaults(
                        loading_parameters
                    ),
                )
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = f"Attempt at reading SHARC trajectory from path `{path}` failed because of original error: {v_e}"
            logging.error(message)
            raise FileNotFoundError(message)

        return loaded_dataset

    def get_units_with_defaults(
        self, unit_overrides: Dict[str, str] | None = None
    ) -> Dict[str, str]:
        """Apply units to the default unit dictionary of the format SHARC

        Args:
            unit_overrides (Dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            Dict[str, str]: The resulting, overridden default units
        """
        from shnitsel.units.definitions import standard_units_of_formats

        # TODO: FIXME: Check if default units are the same for icond and traj
        res_units = standard_units_of_formats["sharc"].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
