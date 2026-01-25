import logging
from typing import Any, TypeVar, overload
import xarray as xr

from shnitsel.data.charge_helpers import guess_molecular_charge
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree import TreeNode

from shnitsel.data.dataset_containers import Trajectory, Frames, wrap_dataset
from shnitsel.io.shared.helpers import LoadingParameters
from shnitsel.data.state_helpers import (
    default_state_name_assigner,
    default_state_type_assigner,
)
from shnitsel.io.shared.variable_flagging import (
    clean_unassigned_variables,
    is_variable_assigned,
)
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults

NodeType = TypeVar("NodeType", bound=TreeNode)
DataType = TypeVar("DataType", bound=xr.Dataset | ShnitselDataset)


@overload
def finalize_loaded_trajectory(
    dataset: TreeNode[Any, DataType],
    loading_parameters: LoadingParameters | None,
) -> TreeNode[Any, DataType]: ...


@overload
def finalize_loaded_trajectory(
    dataset: DataType,
    loading_parameters: LoadingParameters | None,
) -> DataType: ...


@overload
def finalize_loaded_trajectory(
    dataset: None,
    loading_parameters: LoadingParameters | None,
) -> None: ...


def finalize_loaded_trajectory(
    dataset: xr.Dataset | ShnitselDataset | TreeNode[Any, DataType] | None,
    loading_parameters: LoadingParameters | None,
) -> xr.Dataset | ShnitselDataset | TreeNode[Any, xr.Dataset | ShnitselDataset] | None:
    """Function to apply some final postprocessing common to all input routines that allow reading of single trajectories from input formats.

    Parameters
    ----------
    dataset : xr.Dataset | Trajectory | Frames | NodeType | None
        The dataset to perform finalization on. Only updates Dataset, Trajectory and Frames data.
        All other data will be returned unchanged
    loading_parameters : LoadingParameters | None
        Parameters to set some defaults.

    Returns
    -------
    xr.Dataset | Trajectory | Frames | NodeType | None
        The same type as the original `dataset` parameter but potentially with some
        default values and conversions applied.
    """
    if dataset is not None:
        if isinstance(dataset, TreeNode):
            return dataset.map_data(
                lambda x: finalize_loaded_trajectory(
                    x, loading_parameters=loading_parameters
                ),
                keep_empty_branches=True,
            )
        # logging.debug(f"Finalizing: {repr(dataset)}")
        elif isinstance(dataset, (xr.Dataset, ShnitselDataset)):
            rebuild_type = None
            if isinstance(dataset, ShnitselDataset):
                rebuild_type = type(dataset)
                res_dataset = dataset.dataset
            else:
                res_dataset = dataset

            res_dataset = set_state_defaults(res_dataset, loading_parameters)
            # Clean up variables if the variables are not assigned yet.
            res_dataset = clean_unassigned_variables(res_dataset)
            res_dataset = convert_all_units_to_shnitsel_defaults(res_dataset)
            res_dataset = normalize_dataset(res_dataset)

            charge_guess: int | None = None

            # If charge is not set, guess it if we have positional data.
            if "atXYZ" in dataset and (
                "charge" not in dataset
                or isinstance(dataset, ShnitselDataset)
                and not dataset.has_coordinate("charge")
            ):
                try:
                    charge_guess = guess_molecular_charge(res_dataset)
                    if charge_guess != 0:
                        logging.info(
                            "A charge of %d e was guessed for the loaded trajectory with id %s. "
                            "If this was not correct, please set the charge manually using `.set_charge(<charge>)`",
                            charge_guess,
                            str(
                                dataset.trajectory_id
                                if isinstance(dataset, ShnitselDataset)
                                else dataset.attrs.get(
                                    "trajid",
                                    dataset.attrs.get(
                                        "trajectory_id", "<unidentified>"
                                    ),
                                )
                            ),
                        )
                except (RuntimeError, ValueError) as e:
                    logging.debug(
                        "During the guessing of molecular a charge, the following exception was encountered: %s",
                        e,
                    )
                    logging.warning(
                        "The loaded trajectory did not contain charge information and a self-consistent charge could not be guessed. "
                    )
            if rebuild_type:
                tmp_res = rebuild_type(res_dataset)
            else:
                tmp_res = wrap_dataset(res_dataset)

            # Set the charge guess result if we got a result.
            if charge_guess is not None:
                if isinstance(tmp_res, ShnitselDataset):
                    return tmp_res.set_charge(charge_guess)
                else:
                    tmp_res.attrs["charge"] = charge_guess
                    return tmp_res
            else:
                return tmp_res
        else:
            return dataset
    return None


def set_state_defaults(
    dataset: xr.Dataset, loading_parameters: LoadingParameters | None
) -> xr.Dataset:
    """Helper function to apply default settings to dataset variables
    for state names and state types if they have not been assigned at some point
    earlier during the configuration process.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset to set state name and state type information on.
    loading_parameters : LoadingParameters | None
        Currently unused settings to be applied to trajectory import.

    Returns
    -------
    xr.Dataset
        The dataset but with default values for state types and state names
    """
    # TODO: FIXME: apply configured names from loading_parameters

    if is_variable_assigned(dataset.state_types) and is_variable_assigned(
        dataset.state_names
    ):
        logging.debug(
            "Types and names of state already set for dataset in finalization."
        )

        # logging.debug(f"Types: {dataset.state_types}")
        # logging.debug(f"Names: {dataset.state_names}")
        return dataset

    logging.debug("Assigning default state names and/or types.")

    if not is_variable_assigned(dataset.state_types):
        dataset = default_state_type_assigner(dataset)
    if not is_variable_assigned(dataset.state_names):
        dataset = default_state_name_assigner(dataset)
    return dataset


def normalize_dataset(ds: xr.Dataset | ShnitselDataset) -> xr.Dataset:
    """Helper method to perform some standardized renaming operations as well as some
    restructuring, e.g. if a multi-index is missing.

    May also convert some legacy attributes to promoted dimensionless variables.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to normalize according to current Shnitsel standards

    Returns
    -------
    xr.Dataset
        The renamed dataset.
    """
    if isinstance(ds, ShnitselDataset):
        ds = ds.dataset

    if not isinstance(ds, xr.Dataset):
        logging.error(
            "Normalization only supports xr.Dataset or ShnitselDataset entries not %s",
            type(ds),
        )
        return ds

    # Rename legacy uses of `trajid` and `trajid_`
    if 'trajid' in ds.dims:
        ds = ds.rename_dims(trajid='trajectory')
    if 'trajid_' in ds.dims:
        ds = ds.rename_dims(trajid_='trajectory')

    if 'trajid' in ds.coords:
        # print(ds.time)
        if (
            'trajectory' not in ds.coords['trajid'].dims
            or "frame" in ds.coords['trajid'].dims
        ):
            ds = ds.rename_vars(trajid="atrajectory")
        else:
            ds = ds.rename_vars(trajid="trajectory")
        # print(ds.time)

    if 'trajid_' in ds.coords:
        if (
            'trajectory' not in ds.coords['trajid_'].dims
            or "frame" in ds.coords['trajid_'].dims
        ):
            ds = ds.rename_vars(trajid_="atrajectory")
        else:
            ds = ds.rename_vars(trajid_="trajectory")

    if (
        "trajectory" in ds.dims
        and "trajectory" in ds.coords
        and "trajectory" not in ds.indexes
    ):
        ds = ds.set_xindex("trajectory")

    # print(f"Normalization: multiindex {ds.time}")
    # Check if frameset has a multi-index
    if 'frame' in ds.dims and 'time' in ds.coords and 'frame' not in ds.xindexes:
        if 'frame' in ds.coords['time'].dims:
            # Make frame a multi-index
            if "atrajectory" in ds.coords and "frame" in ds.atrajectory.dims:
                ds = ds.set_xindex(["time", "atrajectory"])
            else:
                ds = ds.set_xindex('frame')

    # Turn delta_t and t_max and max_ts into a variable
    if 'time' in ds.coords:
        if 'units' not in ds.time.attrs:
            logging.warning(
                "The imported dataset does not have a time unit set. We will assume the unit to be `fs`."
            )
        time_unit = ds.time.attrs.get('units', 'fs')
        ds.time.attrs['units'] = time_unit
        ds.time.attrs['unitdim'] = "time"
    else:
        time_unit = 'fs'

    for var, has_unit, var_type in [
        ('delta_t', True, float),
        ('max_ts', False, int),
        ('t_max', True, float),
    ]:
        if var not in ds.coords:
            if var in ds.data_vars:
                ds = ds.set_coords(var)
            elif var in ds.attrs:
                new_array = xr.DataArray(
                    var_type(ds.attrs.get(var, -1)), dims=(), name=var
                ).astype(var_type)
                ds = ds.assign_coords({var: new_array})

            if has_unit and var in ds and 'units' not in ds[var].attrs:
                ds[var].attrs['unitdim'] = 'time'
                ds[var].attrs['units'] = time_unit

    if 'astate' in ds.data_vars:
        ds = ds.set_coords("astate")
    if 'sdiag' in ds.data_vars:
        ds = ds.set_coords("sdiag")

    return ds
