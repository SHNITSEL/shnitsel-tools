import logging
from typing import Any, List, Set, Tuple

from shnitsel.data.TrajectoryFormat import Trajectory
import xarray as xr
from shnitsel.units import standard_shnitsel_units
from shnitsel.units.definitions import unit_dimensions

_required_shnitsel_variables = ["energy", "forces", "time"]
_optional_shnitsel_variables = ["atXYZ", "nacs", "dip_perm", "dip_trans", "socs", "state_names",
                                "state_types", "states", "astate", "sdiag", "state2", "statecomb", "phases", "atNames", "atNums", "e_kin", "from", "to", "atom", "state", "direction"]
_required_shnitsel_attributes = ["input_format", "input_type", "completed", "max_ts", "delta_t", "num_singlets", "num_doublets", "num_triplets", "t_max"]
_optional_shnitsel_attributes = ["has_forces", "trajectory_input_path", "trajid", "__original_dataset"]


def check_shnitsel_trajectory_data(trajectory: Trajectory | xr.Dataset, report: bool = False) -> Tuple[Set[str], Set[str], Set[str], Set[str]] | None:
    """Function to check whether all required and only denoted optional variables and meta information is available on a shnitsel-loaded trajectory.

    Args:
        trajectory (Trajectory | xr.Dataset): The trajectory to check for the presence of variables and settings
        report (bool, optional): Whether to raise an error if discrepancies were found. Defaults to False.

    Raises:
        ValueError: If missing variables or attributes or unexpected variables or attributes are encountered in the provided trajectory and ``report=True`` is set.

    Returns:
        Tuple[Set[str], Set[str], Set[str], Set[str]]: The sets of missing required variables, present unexpected variables, missing required attributes and present unexpected attributes.
        None: No discrepancies found.
    """
    missing_required_vars = _required_shnitsel_variables.copy()
    unexpected_vars = []
    missing_required_attrs = _required_shnitsel_attributes.copy()
    unexpected_attrs = []

    for var_name in trajectory.variables:
        # TODO: FIXME: Should we check if the variables are actually set and not just all the same? Should we filter out coordinates separately?
        if var_name in missing_required_vars:
            missing_required_vars.remove(var_name)
        elif var_name not in _optional_shnitsel_variables:
            unexpected_vars.append(var_name)
    # TODO: We should have an even more sophisticated check recursing into the attributes on variables.
    for attr_name in trajectory.attrs:
        if attr_name in missing_required_attrs:
            missing_required_attrs.remove(attr_name)
        elif attr_name not in _optional_shnitsel_attributes:
            unexpected_attrs.append(attr_name)

    # All is well
    if len(missing_required_attrs) == 0 and len(missing_required_vars) == 0 and len(unexpected_vars) == 0 and len(unexpected_attrs) == 0:
        return None

    if report:
        message = "Encountered errors while checking validity of trajectory:"
        if len(missing_required_attrs) > 0:
            message += f"\n Trajectory was missing required variables: {missing_required_vars}"
        if len(missing_required_attrs) > 0:
            message += f"\n Trajectory was missing required attributes: {missing_required_attrs}"
        if len(unexpected_vars) > 0:
            message += f"\n Trajectory had unexpected variables set: {unexpected_vars}"
        if len(unexpected_attrs) > 0:
            message += f"\n Trajectory had unexpected attributes set: {unexpected_attrs}"

        logging.error(message)
        raise ValueError(message)

    return (set(missing_required_vars), set(unexpected_vars), set(missing_required_attrs), set(unexpected_attrs))


def verify_trajectory_format(
    obj: Any, asserted_properties: List[str] | None = None
) -> bool:
    return is_permitted_traj_result(obj) and has_required_properties(
        obj, asserted_properties
    )


def is_permitted_traj_result(obj: Any) -> bool:
    logging.debug(f"Obj has type: {type(obj)}")
    type_check_result = (
        isinstance(obj, xr.Dataset)
        or isinstance(obj, Trajectory)
        or (isinstance(obj, list) and (len(obj) == 0 or isinstance(obj[0], Trajectory)))
    )
    if not type_check_result:
        logging.error(
            "Object is not of type {xr.Dataset}, {Trajectory} or List of {Trajectory}."
        )
    return type_check_result


def has_required_properties(
    traj: List[Trajectory] | Trajectory, asserted_properties: List[str] | None = None
) -> bool:
    res = True
    if isinstance(traj, list):
        for i_traj in traj:
            res = res and has_required_properties(i_traj)
        return res
    else:

        print(traj.variables.keys())
        print(traj["atXYZ"].attrs.keys())
        assert check_shnitsel_trajectory_data(traj, report=True) is None

        check_prop_units = {
            "atXYZ": {"unit": standard_shnitsel_units[unit_dimensions.length]},
            "forces": {"unit": standard_shnitsel_units[unit_dimensions.force]},
            "energy": {"unit": standard_shnitsel_units[unit_dimensions.energy]},
            "dip_trans": {"unit": standard_shnitsel_units[unit_dimensions.dipole]},
            "dip_perm": {"unit": standard_shnitsel_units[unit_dimensions.dipole]},
            "state": {},
            "time": {"unit": standard_shnitsel_units[unit_dimensions.time]},
        }

        available_keys = list(traj.variables.keys())

        if asserted_properties is not None:
            for prop in asserted_properties:
                assert (
                    prop in traj.variables.keys() or prop in traj.coords.keys()
                ), f"Asserted property {prop} of format is missing in resulting trajectory. Only has: {available_keys}"

        for prop in check_prop_units:
            if asserted_properties is not None and prop in asserted_properties:
                assert (
                    prop in traj.variables.keys() or prop in traj.coords.keys()
                ), f"Property {prop} is missing in resulting trajectory. Only has: {available_keys}"
                if "unit" in check_prop_units[prop]:
                    assert (
                        "units" in traj[prop].attrs
                    ), f"Property {prop} has no unit set in resulting trajectory"
                    required_unit = check_prop_units[prop]["unit"]
                    actual_unit = traj[prop].attrs["units"]
                    assert (
                        actual_unit == required_unit
                    ), f"Property {prop} has unit {actual_unit} instead of required unit {required_unit}"

        for var in traj:
            if "unitdim" in traj[var].attrs:
                assert (
                    "units" in traj[var].attrs
                ), f"Variable {var} has property `unitdim` but no `units` set."
                unit_dim = traj[var].attrs["unitdim"]
                actual_unit = traj[var].attrs["units"]
                required_unit = standard_shnitsel_units[unit_dim]
                assert (
                    required_unit == "1" or actual_unit == required_unit
                ), f"Variable {var} has unit {actual_unit} instead of required unit {required_unit}."

        return True
