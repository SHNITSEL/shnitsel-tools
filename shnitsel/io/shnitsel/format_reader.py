from dataclasses import dataclass
import logging
import pathlib
import traceback
from types import UnionType
from typing import Any, Sequence, TypeVar
from typing_extensions import TypeForm

from shnitsel.data.dataset_containers import Trajectory, Frames, InterState, PerState
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
from shnitsel.data.xr_io_compatibility import SupportsFromXrConversion
from ..format_reader_base import FormatInformation, FormatReader
from .parse import read_shnitsel_file
from shnitsel.units.definitions import standard_shnitsel_units
import xarray as xr


@dataclass
class ShnitselFormatInformation(FormatInformation):
    pass


_shnitsel_default_pattern_regex = None
_shnitsel_default_pattern_glob = "*.nc"

DataType = TypeVar("DataType")


class ShnitselFormatReader(FormatReader):
    """Class for providing the Shnitsel format reading functionality in the standardized `FormatReader` interface"""

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
            e for e in path_obj.glob(_shnitsel_default_pattern_glob) if e.is_file()
        ]
        return None if len(res_entries) == 0 else res_entries

    def check_path_for_format_info(
        self, path: PathOptionsType, hints_or_settings: dict | None = None
    ) -> FormatInformation:
        """Check if the `path` is a Shnitsel-style netcdf file

        Designed for a single input file

        Parameters
        ----------
        path : PathOptionsType
            The path to check for being a Shnitsel-style netcdf file
        hints_or_settings : dict | None, optional
            Configuration options provided to the reader by the user, by default None

        Returns
        -------
        FormatInformation
            The object holding all relevant format information for the path contents if it represents a Shnitsel-style
            NetCDF file

        Raises
        ------
        FileNotFoundError
            If the `path` is not a directory.
        FileNotFoundError
            If `path` is a directory but does not satisfy the shnitsel NetCDF requirements
        """
        path_obj: pathlib.Path = make_uniform_path(path)

        _is_request_specific_to_shnitsel = (
            hints_or_settings is not None
            and "kind" in hints_or_settings
            and hints_or_settings["kind"] == "shnitsel"
        )

        if not path_obj.exists() or not path_obj.is_file():
            message = "Path `%(path)s` does not constitute a Shnitsel style trajectory file. Does not exist or is not a file."
            logging.debug(message, {"path": path})
            raise FileNotFoundError(message % {"path": path})

        if not path_obj.suffix.endswith(".nc"):
            message = "Path `%(path)s`` is not a NetCdf file (extension `.nc`)"

            logging.debug(message, {"path": path})
            raise FileNotFoundError(message % {"path": path})

        return ShnitselFormatInformation("shnitsel", "0.1", None, path_obj)

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
        Read a shnitsel-style file from `path`.

        Implements `FormatReader.read_from_path()`.
        Designed for a single input file.

        Parameters
        ----------
        path : pathlib.Path
            Path to a shnitsel-format `.nc` file.
        format_info : FormatInformation
            Format information on the provided `path` that has been previously parsed.
        loading_parameters : LoadingParameters | None, optional
            Loading parameters to e.g. override default state names, units or configure the error reporting behavior, by default None
        expect_dtype : type[DataType] | TypeForm[DataType] | None, optional
            An optional parameter to specify the return type.
            For shnitsel-style NetCDF files, the return type can be pretty much arbitrary.
            This type should specify either the type of the full result or the type of data entries in a hierarchical
            ShnitselDB structure.


        Returns
        -------
        xr.Dataset
        | ShnitselDataset
        | SupportsFromXrConversion
        | TreeNode[Any, ShnitselDataset | SupportsFromXrConversion | xr.Dataset]
        | TreeNode[Any, DataType]
        | Sequence[xr.Dataset | ShnitselDataset | SupportsFromXrConversion]
        | DataType
        | None
            The data stored in the shnitsel-style file.
            Can be any kind of derived data, deserializable data from datasets, raw datasets or a deserialized hierarchy of the desired types.

        Raises
        ------
        ValueError
            Not enough loading information was provided via `path` and `format_info`, e.g. if both are None.
        FileNotFoundError
            Path was not found or was not of appropriate Shnitsel-style NetCDF format
        """
        try:
            # TODO: FIXME: Something is funky with type checks here.
            loaded_dataset_or_tree = read_shnitsel_file(
                path, loading_parameters=loading_parameters
            )

            if isinstance(loaded_dataset_or_tree, xr.Dataset):
                io_metadata = loaded_dataset_or_tree.attrs.get("_shnitsel_io_meta", {})
                try:
                    res = xr_dataset_to_shnitsel_format(
                        loaded_dataset_or_tree, io_metadata
                    )
                except:
                    return loaded_dataset_or_tree

                if expect_dtype is None:
                    return xr_dataset_to_shnitsel_format(
                        loaded_dataset_or_tree, io_metadata
                    )
                elif expect_dtype == xr.Dataset:
                    # Do not convert
                    return loaded_dataset_or_tree
                elif isinstance(loaded_dataset_or_tree, ShnitselDB):
                    return complete_shnitsel_tree(
                        loaded_dataset_or_tree, dtype=expect_dtype
                    )
                else:
                    if is_assignable_to(type(res), expect_dtype):
                        return res
                    else:
                        logging.error(
                            "Could not convert result of type %s to expected type %s. Returning bare xr.Dataset",
                            type(res),
                            expect_dtype,
                        )
                        return loaded_dataset_or_tree
            if isinstance(loaded_dataset_or_tree, xr.DataTree):
                try:
                    db = xarray_datatree_to_shnitsel_tree(
                        loaded_dataset_or_tree, dtype=expect_dtype
                    )
                    return db
                except:
                    logging.error(
                        "Could not convert result of type %s to expected type %s. Returning bare xr.Dataset",
                        type(loaded_dataset_or_tree),
                        expect_dtype,
                    )
                    raise

        except FileNotFoundError as fnf_e:
            raise fnf_e
        except ValueError as v_e:
            message = "Attempt at reading shnitsel file from path `%(path)s` failed because of original error: %(v_e)s.\n Trace: \n %(v_e)s"
            logging.error(
                message, {"path": path, "v_e": v_e, "tb": traceback.format_exc()}
            )
            raise FileNotFoundError(
                message % {"path": path, "v_e": v_e, "tb": traceback.format_exc()}
            )

        return loaded_dataset_or_tree  # type: ignore # We know that the result of read_shnitsel_file is meant to be a ShnitselDB or single Trajectory

    def get_units_with_defaults(
        self, unit_overrides: dict[str, str] | None = None
    ) -> dict[str, str]:
        """Apply units to the default unit dictionary of the format SHNITSEL

        Args:
            unit_overrides (dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            dict[str, str]: The resulting, overridden default units
        """

        res_units = standard_shnitsel_units.copy()

        if unit_overrides is not None:
            res_units.update(unit_overrides)

        return res_units
