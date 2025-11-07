from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
import pathlib
from typing import Callable, Dict, List

from shnitsel.data.TrajectoryFormat import Trajectory
from shnitsel.io.helpers import (
    LoadingParameters,
    PathOptionsType,
    default_state_name_assigner,
    default_state_type_assigner,
    make_uniform_path,
    random_trajid_assigner,
)

import xarray as xr

from shnitsel.io.shared.trajectory_finalization import finalize_loaded_trajectory
from shnitsel.io.shared.trajectory_setup import (
    OptionalTrajectorySettings,
    assign_optional_settings,
)
from shnitsel.io.shared.variable_flagging import mark_variable_assigned


@dataclass
class FormatInformation:
    """Information to keep track of relevant information for"""

    format_name: str = "none"
    version: str = "none"
    trajid: int | None = None
    path: pathlib.Path | None = None


class FormatReader(ABC):
    """Abstract base class for all input formats to define a unified input reader interface.

    Should be subclassed and the functions `check_path_for_format_info()` and `read_from_path()`
    overridden in the subclass
    """

    @abstractmethod
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
        return None

    @abstractmethod
    def check_path_for_format_info(
        self, path: PathOptionsType, hints_or_settings: Dict | None = None
    ) -> FormatInformation:
        """Checks if a path is of a given format and returns a struct containing all relevant info for reading
        the format at this location. Additionally checks configured user settings provided in `hints_or_settings` whether they are
        consistent with the file format.

        Needs to be overridden by each format.

        Args:
            path (os.PathLike):
                The path to look for data from the respective method for.
                Depending on the format, this would need to point to a file or a directory containing the actual
                trajectory information
            hints (Dict|None, optional):
                Potential hints/configuration options provided by the user as input to the reader which can be
                checked for conflicts with the requirements of the format (i.e. requesting a static initial condition from a dynamic trajectory in SHARC).
                Defaults to None

        Raises:
            FileNotFoundError: If required files were not found, i.e. if the path does not actually constitute input data of the denoted format
            ValueError: If the hints/settings provided by the user conflict with the requirements of the format

        Returns:
            FormatInformation:
                A structure containing all of the information relevant to the interpretation or reading of the format.
                Can be used to differentiate different versions of the same format.
                Should be passed to the `read_from_path()` method of the same class.
        """
        return FormatInformation()

    @abstractmethod
    def read_from_path(
        self,
        path: pathlib.Path,
        format_info: FormatInformation,
        loading_parameters: LoadingParameters | None = None,
    ) -> xr.Dataset | None:
        """Method to read a path of the respective format (e.g. ) into a shnitsel-conform trajectory.

        The return value of type `Trajectory` is a wrapper for the raw `xarray.Dataset` read from the `path`.
        This allows provision of extra features like keeping track of the original data while post-processing is performed.

        Args:
            path (pathlib.Path): Path to either the input file or input folder to be read.
            format_info (FormatInformation): Format information previously constructed by `check_path_for_format_info()`. If None, will be constructed by calling `Self.check_path_for_format_info()` first. Defaults to None.
            loading_parameters: (LoadingParameters|None, optional): Loading parameters to e.g. override default state names, units or configure the error reporting behavior

        Raises:
            FileNotFoundError: If required files were not found, i.e. if the path does not actually constitute input data of the denoted format
            ValueError: If the `format_info` provided by the user conflicts with the requirements of the format
            Valueerror: If neither `path` nor `format_info` are provided


        Returns:
            Trajectory: The parsed dataset as wrapper around `xarray.Dataset` to keep track of original data.
        """
        ...

    def read_trajectory(
        self,
        path: PathOptionsType | None,
        format_info: FormatInformation | None = None,
        loading_parameters: LoadingParameters | None = None,
    ) -> Trajectory | None:
        """Wrapper function to perform some potential initialization and finalization on the read trajectory objects.

        Uses the format-specific `self.read_from_path()` method to read the trajectory and then performs some standard post processing on it.


        Args:
            path (PathOptionsType, optional): Path to either the input file or input folder to be read.
            format_info (FormatInformation, optional): Format information previously constructed by `check_path_for_format_info()`. If None, will be constructed by calling `Self.check_path_for_format_info()` first. Defaults to None.
            loading_parameters: (LoadingParameters|None, optional): Loading parameters to e.g. override default state names, units or configure the error reporting behavior

        Raises:
            FileNotFoundError: If required files were not found, i.e. if the path does not actually constitute input data of the denoted format
            ValueError: If the `format_info` provided by the user conflicts with the requirements of the format
            Valueerror: If neither `path` nor `format_info` are provided

        Returns:
            Trajectory|None: Returns a wrapped Trajectory/xr.Dataset object with standard units, only assigned variables remaining and all
                                variables with appropriate attributes.
                                If no result was obtained by the call to `self.read_from_path()`, it will return `None`.
        """

        loading_parameters = self.get_loading_parameters_with_defaults(
            loading_parameters
        )

        path_obj: pathlib.Path = make_uniform_path(path)

        if path_obj is None:
            raise ValueError(
                "Not sufficient `path` information provided. Please set the `path` parameter"
            )

        if path_obj is not None and format_info is None:
            format_info = self.check_path_for_format_info(path_obj)
        if path_obj is None and format_info is not None:
            path_obj = format_info.path
        if path_obj is None or format_info is None:
            raise ValueError(
                "Either `path` or `format_info` needs to be provided and the other must be derivable from the other information. Not enough information provided for loading trajectory."
            )

        if not path_obj.exists():
            raise FileNotFoundError(f"Path at {path_obj} does not exist.")

        res = self.read_from_path(path_obj, format_info, loading_parameters)

        if res is not None:
            # Set some optional settings.
            optional_settings = OptionalTrajectorySettings()
            optional_settings.trajectory_input_path = path_obj.as_posix()
            
            assign_optional_settings(res, optional_settings)

            # If trajid has been extracted from the input path, set it
            if format_info is not None:
                # If trajid has been extracted from the input path, set it
                if loading_parameters.trajectory_id is not None:
                    # the trajectory_id assignment should have been transformed into a callable
                    traj_id_assigner: Callable[[pathlib.Path], int] = (
                        loading_parameters.trajectory_id
                    )  # type: ignore
                    optional_settings.trajid = traj_id_assigner(
                        format_info.path if format_info.path is not None else path_obj
                    )

                    if (
                        traj_id_assigner == random_trajid_assigner
                        and format_info.trajid is not None
                    ):
                        optional_settings.trajid = format_info.trajid
                elif format_info.trajid is not None:
                    optional_settings.trajid = format_info.trajid

                # Assign state types if provided
                if loading_parameters.state_types is not None:
                    state_types_assigner: Callable[[xr.Dataset], xr.Dataset] = (
                        loading_parameters.state_types
                    )  # type: ignore
                    res = state_types_assigner(res)

                # Assign state names if provided
                if loading_parameters.state_names is not None:
                    state_names_assigner: Callable[[xr.Dataset], xr.Dataset] = (
                        loading_parameters.state_names
                    )  # type: ignore
                    res = state_names_assigner(res)

            assign_optional_settings(res, optional_settings)

            return finalize_loaded_trajectory(res, loading_parameters)
        else:
            return res

    @abstractmethod
    def get_units_with_defaults(
        self, unit_overrides: Dict[str, str] | None = None
    ) -> Dict[str, str]:
        """Apply units to the default unit dictionary of the format

        Args:
            unit_overrides (Dict[str, str] | None, optional): Units denoted by the user to override format default settings. Defaults to None.

        Raises:
            NotImplementedError: The class does not provide this functionality yet

        Returns:
            Dict[str, str]: The resulting, overridden default units
        """
        raise NotImplementedError()

    def get_loading_parameters_with_defaults(
        self, base_loading_parameters: LoadingParameters | None
    ) -> LoadingParameters:
        """Populate loading parameters with default settings for this format

        Args:
            base_loading_parameters (LoadingParameters | None): User-provided parameter overrides

        Returns:
            LoadingParameters: The default parameters modified by user overrides
        """

        # Transform different options of input settings into a callable that assigns the values.
        get_state_name_callable: Callable[[xr.Dataset], xr.Dataset] | None = None
        get_state_types_callable: Callable[[xr.Dataset], xr.Dataset] | None = None
        get_traj_id_callable: Callable[[pathlib.Path], int] | None = None
        if base_loading_parameters is not None:
            state_names_override = base_loading_parameters.state_names
            if state_names_override is not None:
                if callable(state_names_override):
                    get_state_name_callable = state_names_override
                elif isinstance(state_names_override, list):

                    def tmp_state_assigner(dataset: xr.Dataset) -> xr.Dataset:
                        dataset = dataset.assign_coords(
                            {
                                "state_names": (
                                    "state",
                                    state_names_override,
                                    dataset.state_names.attrs,
                                )
                            }
                        )
                        mark_variable_assigned(dataset.state_names)
                        return dataset

                    get_state_name_callable = tmp_state_assigner
                else:
                    get_state_name_callable = default_state_name_assigner
            state_types_override = base_loading_parameters.state_types
            if state_types_override is not None:
                if callable(state_types_override):
                    get_state_types_callable = state_types_override
                elif isinstance(state_types_override, list):

                    def tmp_state_assigner(dataset: xr.Dataset) -> xr.Dataset:
                        dataset = dataset.assign_coords(
                            {
                                "state_types": (
                                    "state",
                                    state_types_override,
                                    dataset.state_types.attrs,
                                )
                            }
                        )
                        mark_variable_assigned(dataset.state_types)
                        return dataset

                    get_state_types_callable = tmp_state_assigner
                else:
                    get_state_types_callable = default_state_type_assigner
            trajid_override = base_loading_parameters.trajectory_id
            if trajid_override is not None:
                if callable(trajid_override):
                    get_traj_id_callable = trajid_override
                elif isinstance(trajid_override, dict):

                    def tmp_trajid_assigner(path: pathlib.Path) -> int:
                        path_str = path.absolute().as_posix()
                        if path_str in trajid_override:
                            return trajid_override[path_str]
                        else:
                            return -1

                    get_traj_id_callable = tmp_trajid_assigner
                else:
                    get_traj_id_callable = random_trajid_assigner

        return LoadingParameters(
            self.get_units_with_defaults(
                base_loading_parameters.input_units
                if base_loading_parameters is not None
                else None
            ),
            (
                base_loading_parameters.error_reporting
                if base_loading_parameters is not None
                else "raise"
            ),
            get_traj_id_callable,
            get_state_types_callable,
            get_state_name_callable,
        )
