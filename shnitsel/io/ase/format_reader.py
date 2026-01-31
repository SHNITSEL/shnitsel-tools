from dataclasses import dataclass
import logging
import pathlib
import traceback
from types import UnionType
from typing import Any, Sequence, TypeVar
from typing_extensions import TypeForm

from shnitsel.data.dataset_containers import (
    Trajectory,
    Frames,
    InterState,
    PerState,
    wrap_dataset,
)
from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.dataset_containers.xr_conversion import xr_dataset_to_shnitsel_format
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.xr_conversion import xarray_datatree_to_shnitsel_tree
from shnitsel.data.helpers import is_assignable_to
from shnitsel.data.tree import (
    DataGroup,
    ShnitselDB,
    DataLeaf,
    CompoundGroup,
    complete_shnitsel_tree,
)
from shnitsel.io.shared.helpers import (
    LoadingParameters,
    PathOptionsType,
    make_uniform_path,
)
from shnitsel.io.shared.trajectory_finalization import normalize_dataset
from shnitsel.data.xr_io_compatibility import SupportsFromXrConversion
from ..format_reader_base import FormatInformation, FormatReader
from .parse import read_ase
import xarray as xr


@dataclass
class ASEFormatInformation(FormatInformation):
    pass


_ase_default_pattern_regex = "*.db"
_ase_default_pattern_glob = "*.db"

DataType = TypeVar("DataType")


class ASEFormatReader(FormatReader):
    """Class for providing the ASE format reading functionality in the standardized `FormatReader` interface"""

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
        res_entries = [
            e for e in path_obj.glob(_ase_default_pattern_glob) if e.is_file()
        ]
        return None if len(res_entries) == 0 else res_entries

    def check_path_for_format_info(
        self, path: PathOptionsType, hints_or_settings: dict | None = None
    ) -> FormatInformation:
        """Check if the `path` is a ASE-style `.db` file

        Designed for a single input file

        Parameters
        ----------
        path : PathOptionsType
            The path to check for being a ASE-style `.db` file
        hints_or_settings : dict | None, optional
            Configuration options provided to the reader by the user, by default None

        Returns
        -------
        FormatInformation
            The object holding all relevant format information for the path contents if it represents a ASE-style
            NetCDF file

        Raises
        ------
        FileNotFoundError
            If the `path` is not a directory.
        FileNotFoundError
            If `path` is a directory but does not satisfy the ASE file format input requirements
        """
        path_obj: pathlib.Path = make_uniform_path(path)

        _is_request_specific_to_ASE = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "ase"
        )

        if not path_obj.exists() or not path_obj.is_file():
            message = "Path `%(path)s` does not constitute a ASE style trajectory file. Does not exist or is not a file."
            logging.debug(message, {"path": path})
            raise FileNotFoundError(message % {"path": path})

        if not path_obj.suffix.endswith(".db"):
            message = "Path `%(path)s`` is not a NetCdf file (extension `.db`)"

            logging.debug(message, {"path": path})
            raise FileNotFoundError(message % {"path": path})

        return ASEFormatInformation("ASE", "1.0", None, path_obj)

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
        """
        Read an ASE-style file from `path`.

        Implements `FormatReader.read_from_path()`.
        Designed for a single input file.

        Parameters
        ----------
        path : pathlib.Path
            Path to an ASE-format `.db` file.
        format_info : FormatInformation
            Format information on the provided `path` that has been previously parsed.
        loading_parameters : LoadingParameters | None, optional
            Loading parameters to e.g. override default state names, units or configure the error reporting behavior, by default None
        expect_dtype : type[DataType] | UnionType | None, optional
            An optional parameter to specify the return type.
            For ASE-style db files, the return type can be either a single frameset (Frames), a trajectory (Trajectory) a collection/list
            of such sets or a blank xr.Dataset if conversion did not succeed.
            If for some reason a datatype that is convertible to a Shnitsel-type from arbitrary xr.Dataset containers (SupportsFromXrConversion)
            is contained within the database, it will be deserialized.
            This type should specify either the type of the full result or the type of data entries in a hierarchical Shnitsel tree structure.


        Returns
        -------
        xr.Dataset
        | ShnitselDataset
        | SupportsFromXrConversion
        | TreeNode[Any, ShnitselDataset | SupportsFromXrConversion | xr.Dataset]
        | TreeNode[Any, DataType]
        | DataType
        | None
            The data stored in the ASE-style file.
            Can be either a single frameset (Frames), a trajectory (Trajectory) a collection/list
            of such sets or a blank xr.Dataset if conversion did not succeed.
            If for some reason a datatype that is convertible to a Shnitsel-type from arbitrary xr.Dataset containers (SupportsFromXrConversion)
            is contained within the database, it will be deserialized.

        Raises
        ------
        ValueError
            Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
        FileNotFoundError
            Path was not found or was not of appropriate ASE-style `.db` format
        """
        try:
            # TODO: FIXME: Something is funky with type checks here.
            loaded_dataset_collection_or_tree = read_ase(
                path, loading_parameters=loading_parameters
            )
        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = "Attempt at reading ASE file from path `%(path)s` failed because of original error: %(v_e)s.\n Trace: \n %(v_e)s"
            logging.error(
                message, {"path": path, "v_e": v_e, "tb": traceback.format_exc()}
            )
            raise FileNotFoundError(
                message % {"path": path, "v_e": v_e, "tb": traceback.format_exc()}
            )

        return loaded_dataset_collection_or_tree  # type: ignore # We know that the result of read_shnitsel_file is meant to be a ShnitselDB or single Trajectory

    def get_units_with_defaults(
        self, unit_overrides: dict[str, str] | None = None
    ) -> dict[str, str]:
        """Apply units to the default unit dictionary of the format ASE

        Args:
            unit_overrides (dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            dict[str, str]: The resulting, overridden default units
        """
        from shnitsel.units.definitions import standard_units_of_formats

        res_units = standard_units_of_formats["ase"].copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
