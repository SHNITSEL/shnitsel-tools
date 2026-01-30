from dataclasses import dataclass
from glob import glob
import logging
import pathlib
import re
import sys
import traceback
from types import UnionType
from typing import Any, Sequence, TypeVar
from typing_extensions import TypeForm
from shnitsel.data.dataset_containers import Trajectory, Frames

from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.node import TreeNode
from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
from shnitsel.data.xr_io_compatibility import SupportsFromXrConversion
from ..format_reader_base import FormatInformation, FormatReader
from .parse import parse_pyrai2md

import xarray as xr


@dataclass
class PyrAI2mdFormatInformation(FormatInformation):
    energy_file_path: pathlib.Path | None = None
    log_file_path: pathlib.Path | None = None
    pass


DataType = TypeVar("DataType")


class PyrAI2mdFormatReader(FormatReader):
    """Class for providing the PyrAI2md format reading functionality in the standardized `FormatReader` interface"""

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

        res_entries = [e for e in path_obj.glob("*") if e.is_dir()]
        return None if len(res_entries) == 0 else res_entries

    def check_path_for_format_info(
        self, path: PathOptionsType, hints_or_settings: dict | None = None
    ) -> FormatInformation:
        """Check if the `path` is a PyrAI2MD-style output directory.

        Designed for a single input trajectory.

        Parameters
        ----------
        path : PathOptionsType
            The path to check for PyrAI2MD data
        hints_or_settings : dict | None, optional
            Configuration options provided to the reader by the user, by default None

        Returns
        -------
        FormatInformation
            The object holding all relevant format information for the path contents if
            it matches the PyrAI2MD format

        Raises
        ------
        FileNotFoundError
            If the `path` is not a directory.
        FileNotFoundError
            If `path` is a directory but does not contain the required PyrAI2MD output files
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
                "Path `%(path)s` does not constitute a PyrAI2md style output directory: Expected to find a single file ending with '.md.energies' "
                "but found %(n)d files: %(md_energies_paths)s"
            )
            params = {"path": path, "n": n, "md_energies_paths": md_energies_paths}
            logging.debug(message, params)
            raise FileNotFoundError(message % params)

        energy_file_path = path_obj / md_energies_paths[0]

        log_paths = glob(
            "*.log",
            root_dir=path_obj,
        )
        if (n := len(md_energies_paths)) != 1:
            message = (
                "Path `%(path)s` does not constitute a PyrAI2md style output directory: Expected to find a single file ending with '.log' "
                "but found %(n)d files: %(log_paths)s"
            )
            params = {"path": path, "n": n, "log_paths": log_paths}
            logging.debug(message, params)
            raise FileNotFoundError(message % params)

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
        """Read a PyrAI2MD-style trajcetory from path at `path`.
        Implements `FormatReader.read_from_path()`.

        Designed for a single input trajectory.

        Parameters
        ----------
        path : pathlib.Path
            Path to a PyrAI2MD-format directory.
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
            Path was not found or was not of appropriate PyrAI2MD format
        """

        if expect_dtype is not None:
            logging.debug(
                "The NewtonX format reader only supports dtypes `Trajectory`, `Frames` or `xr.Dataset` "
            )

        try:
            loaded_dataset = parse_pyrai2md(path, loading_parameters=loading_parameters)
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = "Attempt at reading PyRAI2md trajectory from path `%(path)s` failed because of original error: %(v_e)s.\n Trace: \n %(v_e)s"
            logging.error(
                message, {"path": path, "v_e": v_e, "tb": traceback.format_exc()}
            )
            raise FileNotFoundError(message)

        return loaded_dataset

    def get_units_with_defaults(
        self, unit_overrides: dict[str, str] | None = None
    ) -> dict[str, str]:
        """
        Apply units to the default unit dictionary of the format PyrAI2md

        Parameters
        ----------
        unit_overrides : dict[str, str] | None, optional
            Units denoted by the user to override format default settings, by default None

        Returns
        -------
        dict[str, str]
            The resulting, overridden default units
        """
        from shnitsel.units.definitions import standard_units_of_formats

        res_units = standard_units_of_formats["pyrai2md"].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
