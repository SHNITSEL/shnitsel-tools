from dataclasses import dataclass
import logging
import pathlib
import re
import traceback
from types import UnionType
from typing import Any, Sequence, TypeVar
from typing_extensions import TypeForm

from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
from shnitsel.data.xr_io_compatibility import SupportsFromXrConversion
from ..format_reader_base import FormatInformation, FormatReader
from .parse_trajectory import read_traj
from .parse_initial_conditions import read_iconds_individual

import xarray as xr


@dataclass
class SHARCDynamicFormatInformation(FormatInformation):
    pass


@dataclass
class SHARCInitialFormatInformation(FormatInformation):
    pass


@dataclass
class SHARCMultiInitialFormatInformation(FormatInformation):
    list_of_iconds: list | None = None


_sharc_default_pattern_regex = re.compile(r"(?P<dynstat>TRAJ|ICOND)_(?P<trajid>\d+)")
_sharc_default_pattern_glob_traj = "TRAJ_*"
_sharc_default_pattern_glob_icond = "ICOND_*"

DataType = TypeVar("DataType")


class SHARCFormatReader(FormatReader):
    """Class for providing the SHARC format reading functionality in the standardized `FormatReader` interface"""

    def find_candidates_in_directory(
        self, path: PathOptionsType
    ) -> list[pathlib.Path] | None:
        """
        Function to return a all potential matches for the current file format  within a provided directory at `path`.
        Parameters

        ----------
        path : PathOptionsType
            The path to a directory to check for potential candidate files or subdirectories

        Returns
        -------
        list[pathlib.Path]
            A list of paths that should be checked in detail for whether they represent the format of this FormatReader.
        None
            If no potential candidates were found
        """
        path_obj = make_uniform_path(path)

        base_format_info = super().check_path_for_format_info(path_obj)

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
        self, path: PathOptionsType, hints_or_settings: dict | None = None
    ) -> FormatInformation:
        """Check if the `path` is a SHARC-style output directory.

        Designed for a single input trajectory.

        Parameters
        ----------
        path : PathOptionsType
            The path to check for SHARC data
        hints_or_settings : dict | None, optional
            Configuration options provided to the reader by the user, by default None

        Returns
        -------
        FormatInformation
            The object holding all relevant format information for the path contents if
            it matches the SHARC format

        Raises
        ------
        FileNotFoundError
            If the `path` is not a directory.
        FileNotFoundError
            If `path` is a directory but does not contain the required SHARC output files
        """
        path_obj: pathlib.Path = make_uniform_path(path)

        base_format_info = super().check_path_for_format_info(
            path_obj, hints_or_settings
        )

        _is_request_specific_to_sharc = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "sharc"
        )

        if not path_obj.exists() or not path_obj.is_dir():
            message = "Path `%(path)s` does not constitute a SHARC style output directory: Does not exist or is not a directory."
            logging.debug(
                message,
                {"path": path},
            )
            raise FileNotFoundError(message % {"path": path})

        dontanalyze_file_path = path_obj / "DONT_ANALYZE"

        if dontanalyze_file_path.exists() and dontanalyze_file_path.is_file():
            message = "The path `%(path)s` does contain a `DONT_ANALYZE` file and will therefore be skipped. Please remove that file if you want the directory to be read."
            logging.warning(
                message,
                {"path": path},
            )
            raise FileNotFoundError(message % {"path": path})

        # Check if dynamic SHARC format satisfied
        is_dynamic = False
        format_information: FormatInformation | None = None
        try:
            input_file_path = path_obj / "input"
            input_dat_path = path_obj / "output.dat"
            input_xyz_path = path_obj / "output.xyz"

            for file in [input_file_path, input_dat_path, input_xyz_path]:
                if not file.is_file():
                    message = "Input directory `%(path)s` is missing `%(file)s`"

                    logging.debug(message, {"path": path, "file": file})
                    raise FileNotFoundError(message % {"path": path, "file": file})
            is_dynamic = True
            format_information = SHARCDynamicFormatInformation(
                "sharc", "unkown", None, path_obj
            )
            logging.debug(
                "Input directory `%(path)s` fulfils data requirements of dynamic SHARC trajectory",
                {"path": path},
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
                message = "Input directory `%(path)s` is missing `QM.out` or both `QM.log` and `QM.in`"
                logging.debug(message, {"path": path})
                raise FileNotFoundError(message % {"path": path})

            # list_of_initial_condition_paths = list_iconds(path_obj)
            is_static = True
            format_information = SHARCInitialFormatInformation(
                "sharc",
                "unkown",
                None,
                path_obj,  # , list_of_initial_condition_paths
            )
            logging.debug(
                "Input directory `%(path)s` fulfils data requirements of SHARC Initial Conditions",
                {"path": path},
            )
        except Exception as e:
            static_check_error = e

        if is_dynamic and is_static:
            message = (
                "Input directory `%(path)s` contains both static initial conditions and dynamic trajectory data of type SHARC."
                "Please only point to a directory containing exactly one of the two kinds of data"
            )
            logging.debug(message, {"path": path})
            raise ValueError(message % {"path": path})
        if format_information is None:
            message = (
                "Input directory `%(path)s` contains neither static initial conditions nor dynamic trajectory data of type SHARC."
                "Please point to a directory containing exactly one of the two kinds of data"
            )
            logging.debug(message, {"path": path})
            raise FileNotFoundError(message % {"path": path})

        # Try and extract a trajectory ID from the path name
        match_attempt = _sharc_default_pattern_regex.match(path_obj.name)

        if match_attempt:
            path_based_trajid = match_attempt.group("trajid")
            format_information.trajid = int(path_based_trajid)
        else:
            format_information.trajid = base_format_info.trajid

        return format_information

    def read_from_path(
        self,
        path: pathlib.Path,
        *,
        format_info: FormatInformation,
        loading_parameters: LoadingParameters | None = None,
        expect_dtype: type[DataType] | UnionType | None = None,
    ) -> (
        xr.Dataset
        | xr.DataArray
        | ShnitselDataset
        | SupportsFromXrConversion
        | TreeNode[
            Any, ShnitselDataset | SupportsFromXrConversion | xr.Dataset | xr.DataArray
        ]
        | TreeNode[Any, DataType]
        | Sequence[xr.Dataset | ShnitselDataset | SupportsFromXrConversion]
        | DataType
        | None
    ):
        """Read a SHARC-style trajcetory from path at `path`.
        Implements `FormatReader.read_from_path()`.

        Designed for a single input trajectory.

        Parameters
        ----------
        path : pathlib.Path
            Path to a SHARC-format directory.
        format_info : FormatInformation
            Format information on the provided `path` that has been previously parsed.
        loading_parameters : LoadingParameters | None, optional
            Loading parameters to e.g. override default state names, units or configure the error reporting behavior, by default None
        expect_dtype : type[DataType] | UnionType | None, optional
            An optional parameter to specify the return type.
            For this class, it should be `xr.Dataset`, `Trajectory` or `Frames`, by default None

        Returns
        -------
        xr.Dataset | Trajectory | Frames | None
            The loaded Shnitsel-conforming trajectory.

        Raises
        ------
        ValueError
            Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
        FileNotFoundError
            Path was not found or was not of appropriate SHARC format
        """

        if expect_dtype is not None:
            logging.debug(
                "The SHARC format reader only supports dtypes `Trajectory`, `Frames` or `xr.Dataset` "
            )

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
                    loading_parameters=loading_parameters,
                )
            else:
                loaded_dataset = read_iconds_individual(
                    path,
                    loading_parameters=loading_parameters,
                )
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = "Attempt at reading SHARC trajectory from path `%(path)s` failed because of original error: %(v_e)s.\n Trace: \n %(v_e)s"
            logging.error(
                message, {"path": path, "v_e": v_e, "tb": traceback.format_exc()}
            )
            raise FileNotFoundError(message)

        return loaded_dataset

    def get_units_with_defaults(
        self, unit_overrides: dict[str, str] | None = None
    ) -> dict[str, str]:
        """
        Apply units to the default unit dictionary of the format SHARC

        Parameters
        ----------
        unit_overrides : Dict[str, str] | None, optional
            Units denoted by the user to override format default settings., by default None

        Returns
        -------
        dict[str, str]
            The resulting, overridden default units
        """
        from shnitsel.units.definitions import standard_units_of_formats

        # TODO: FIXME: Check if default units are the same for icond and traj
        res_units = standard_units_of_formats["sharc"].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
