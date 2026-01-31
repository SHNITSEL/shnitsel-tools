from dataclasses import dataclass
import logging
import pathlib
import re
import traceback
from types import UnionType
from typing import Any, Sequence, TypeVar
from typing_extensions import TypeForm

from shnitsel.data.dataset_containers import Trajectory, Frames
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
from shnitsel.data.xr_io_compatibility import SupportsFromXrConversion
from ..format_reader_base import FormatInformation, FormatReader
from .parse import parse_newtonx

import xarray as xr


@dataclass
class NewtonXFormatInformation(FormatInformation):
    nx_log_path: pathlib.Path | None = None
    positions_file_path: pathlib.Path | None = None
    pass


DataType = TypeVar("DataType")

_newtonx_default_pattern_regex = re.compile(r"TRAJ(?P<trajid>\d+)")
_newtonx_default_pattern_glob = r"TRAJ*"


class NewtonXFormatReader(FormatReader):
    """Class for providing the NewtonX format reading functionality in the standardized `FormatReader` interface"""

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
        self, path: PathOptionsType, hints_or_settings: dict | None = None
    ) -> FormatInformation:
        """
        Check if the `path` is a NewtonX-style output directory.

        Designed for a single input trajectory.

        Parameters
        ----------
        path : PathOptionsType
            The path to check for NewtonX data
        hints_or_settings : dict | None, optional
            Configuration options provided to the reader by the user, by default None

        Returns
        -------
        FormatInformation
            The object holding all relevant format information for the path contents if
            it matches the NewtonX format

        Raises
        ------
        FileNotFoundError
            If the `path` is not a directory.
        FileNotFoundError
            If `path` is a directory but does not contain the required NewtonX output files
        """
        path_obj: pathlib.Path = make_uniform_path(path)

        base_format_info = super().check_path_for_format_info(
            path_obj, hints_or_settings
        )

        _is_request_specific_to_NewtonX = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and (
                hints_or_settings["kind"] == "nx"
                or hints_or_settings["kind"] == "newtonx"
            )
        )

        nx_log_path = path_obj / "RESULTS" / "nx.log"
        nx_dynxyz_path = path_obj / "RESULTS" / "dyn.xyz"
        nx_dynout_path = path_obj / "RESULTS" / "dyn.out"
        nx_endat_path = path_obj / "RESULTS" / "en.dat"

        num_hits = 0
        file_not_found_miss = None
        for file in [nx_log_path, nx_dynxyz_path, nx_dynout_path, nx_endat_path]:
            if not file.is_file():
                message = "Input directory is missing %(file)s"
                params = {"file": file}
                logging.debug(message, params)
                file_not_found_miss = FileNotFoundError(message % params)
            else:
                num_hits += 1

        if num_hits < 3 and file_not_found_miss is not None:
            raise file_not_found_miss

        format_information = NewtonXFormatInformation(
            "newtonx", "unkown", None, path_obj, nx_log_path, nx_dynxyz_path
        )

        # Try and extract a trajectory ID from the path name
        match_attempt = _newtonx_default_pattern_regex.match(path_obj.name)

        if match_attempt:
            path_based_trajid = match_attempt.group("trajid")
            format_information.trajid = int(path_based_trajid)
            logging.info(
                "Assigning id %(id)d to trajectory", {"id": format_information.trajid}
            )
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
        """Read a NewtonX-style trajcetory from path at `path`. Implements `FormatReader.read_from_path()`

        Parameters
        ----------
        path : pathlib.Path
            Path to a NewtonX-format directory.
        format_info : FormatInformation
            Format information on the provided `path` that has been previously parsed.
        loading_parameters : LoadingParameters | None, optional
            Loading parameters to e.g. override default state names, units or configure the error reporting behavior, by default None
        expect_dtype : type[DataType] | TypeForm[DataType] | None, optional
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
            Path was not found or was not of appropriate NewtonX format
        """

        if expect_dtype is not None:
            logging.debug(
                "The NewtonX format reader only supports dtypes `Trajectory`, `Frames` or `xr.Dataset` "
            )

        try:
            loaded_dataset = parse_newtonx(path, loading_parameters=loading_parameters)
        except FileNotFoundError:
            raise
        except ValueError as v_e:
            message = "Attempt at reading NewtonX trajectory from path `%(path)s` failed because of original error: %(v_e)s.\n Trace: \n %(v_e)s"
            logging.error(
                message, {"path": path, "v_e": v_e, "tb": traceback.format_exc()}
            )
            raise

        return loaded_dataset

    def get_units_with_defaults(
        self, unit_overrides: dict[str, str] | None = None
    ) -> dict[str, str]:
        """
        Apply units to the default unit dictionary of the format NewtonX

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

        res_units = standard_units_of_formats["newtonx"].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
