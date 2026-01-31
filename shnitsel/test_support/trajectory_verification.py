import logging
from typing import Any, List, Set, Tuple

from shnitsel.data.dataset_containers import Frames, Trajectory, wrap_dataset
import xarray as xr
from shnitsel.data.tree.node import TreeNode
from shnitsel.units import standard_shnitsel_units
from shnitsel.units.definitions import unit_dimensions

# TODO: FIXME: Some attributes are turned into variables on merged trajectories.
_required_shnitsel_variables = ["energy", "time", "delta_t", "t_max", "max_ts"]
_optional_shnitsel_variables = [
    "atXYZ",
    "nacs",
    "dip_perm",
    "dip_trans",
    "forces",
    "socs",
    "state_names",
    "state_types",
    "state_charges",
    "state_magnetic_number",
    "state_degeneracy_group",
    "astate",
    "sdiag",
    "state2",
    "phases",
    "atNames",
    "atNums",
    "e_kin",
    "velocities",
    "statecomb",
    "from",
    "to",
    "atom",
    "state",
    "direction",
    "trajid",
    'full_statecomb',
    'full_statecomb_from',
    'full_statecomb_to',
    'charge'
]
_required_shnitsel_attributes = [
    "input_format",
    "input_format_version",
    "input_type",
    "completed",
    "num_singlets",
    "num_doublets",
    "num_triplets",
]
_optional_shnitsel_attributes = [
    "has_forces",
    "DataTree_Level",
    "trajectory_input_path",
    "trajid",
    "__original_dataset",
    "trajectory_id",
    "_shnitsel_setup_for_cleanup",
    "is_multi_trajectory",
    "misc_input_settings",
    "theory_basis_set",
    "est_level",
    "max_ts",
    "delta_t",
    "t_max",
]


def check_shnitsel_trajectory_data(
    trajectory: Frames | Trajectory | xr.Dataset, report: bool = False
) -> tuple[set[str], set[str], set[str], set[str]] | None:
    """Function to check whether all required and only denoted optional variables and meta information is available 
    on a shnitsel-loaded trajectory.

    Parameters
    ----------
    trajectory : Frames | Trajectory | xr.Dataset
        The trajectory to check for the presence of variables and settings
    report : bool, optional
        Whether to raise an error if discrepancies were found. Defaults to False.

    Raises
    ------
    ValueError
        If missing variables or attributes or unexpected variables or attributes are encountered in the provided 
        trajectory and ``report=True`` is set.

    Returns
    -------
    tuple[set[str], set[str], set[str], set[str]]
        The sets of missing required variables, present unexpected variables, missing required attributes and present 
        unexpected attributes.
    None
        No discrepancies found.
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
        elif (
            not attr_name.startswith("__")
            and attr_name not in _optional_shnitsel_attributes
        ):
            unexpected_attrs.append(attr_name)

    # All is well
    if (
        len(missing_required_attrs) == 0
        and len(missing_required_vars) == 0
        and len(unexpected_vars) == 0
        and len(unexpected_attrs) == 0
    ):
        return None

    if report:
        message = "Encountered errors while checking validity of trajectory:"
        if len(missing_required_vars) > 0:
            message += (
                f"\n Trajectory was missing required variables: {missing_required_vars}"
            )
        if len(missing_required_attrs) > 0:
            message += f"\n Trajectory was missing required attributes: {missing_required_attrs}"
        if len(unexpected_vars) > 0:
            message += f"\n Trajectory had unexpected variables set: {unexpected_vars}"
        if len(unexpected_attrs) > 0:
            message += (
                f"\n Trajectory had unexpected attributes set: {unexpected_attrs}"
            )

        logging.error(message)
        raise ValueError(message)

    return (
        set(missing_required_vars),
        set(unexpected_vars),
        set(missing_required_attrs),
        set(unexpected_attrs),
    )


def verify_trajectory_format(
    obj: Trajectory
    | Frames
    | xr.Dataset
    | List[Trajectory | Frames | xr.Dataset]
    | TreeNode[Any, Trajectory | Frames | xr.Dataset],
    asserted_properties: List[str] | None = None,
) -> bool:
    """Verify whether the data in `obj` has all required properties of a ShnitselTools trajectory.

    Parameters
    ----------
    obj : Any
        Either a read trajectory (Trajectory|Frames|xr.Dataset) or a Tree structure with the appropriate
        entries in its data nodes.
    asserted_properties : List[str] | None, optional
        The list of keys that are required to be in the trajectory, by default None

    Returns
    -------
    bool
        True if the format is as asserted.
    """
    from shnitsel.data.tree import TreeNode

    if isinstance(obj, TreeNode):
        mapped_checks = obj.map_data(
            lambda x: verify_trajectory_format(
                x, asserted_properties=asserted_properties
            ),
            keep_empty_branches=True,
            dtype=bool,
        )
        return all(mapped_checks.collect_data())
    else:
        return is_permitted_traj_result(obj) and has_required_properties(
            obj, asserted_properties
        )


def is_permitted_traj_result(obj: Any) -> bool:
    """Check the type returned by the `read()` function.

    Parameters
    ----------
    obj : Any
        The object returned by `read`.

    Returns
    -------
    bool
        Whether the type is in agreement with the intended type model
    """
    from shnitsel.data.tree import ShnitselDB

    logging.debug(f"Obj has type: {type(obj)}")
    type_check_result = isinstance(
        obj, (xr.Dataset, Trajectory, Frames, ShnitselDB)
    ) or (
        isinstance(obj, list)
        and (
            len(obj) == 0
            or all(isinstance(d, (Trajectory, Frames, xr.Dataset)) for d in obj)
        )
    )

    if not type_check_result:
        logging.error(
            "Object is not of type {xr.Dataset}, {Trajectory}, {Frames} or List thereof."
        )
    return type_check_result


def has_required_properties(
    traj: List[Trajectory | Frames | xr.Dataset] | Trajectory | Frames | xr.Dataset,
    asserted_properties: List[str] | None = None,
) -> bool:
    """Function to check whether a certain set of properties is set on the trajectory

    Parameters
    ----------
    traj : List[Trajectory] | Trajectory
        The parsed Trajectory object(s)
    asserted_properties : List[str] | None, optional
        The keys of properties that are required to be set on the trajectory, by default None

    Returns
    -------
    bool
        Whether all properties were present.
    """
    res = True
    if isinstance(traj, list):
        return all(has_required_properties(i_traj) for i_traj in traj)
    else:
        wrapped_ds = wrap_dataset(traj, (Trajectory | Frames))
        # print(traj.variables.keys())
        # print(traj["atXYZ"].attrs.keys())
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
                assert prop in traj.variables.keys() or prop in traj.coords.keys(), (
                    f"Asserted property {prop} of format is missing in resulting trajectory. Only has: {available_keys}"
                )

        for prop in check_prop_units:
            if asserted_properties is not None and prop in asserted_properties:
                assert prop in traj.variables.keys() or prop in traj.coords.keys(), (
                    f"Property {prop} is missing in resulting trajectory. Only has: {available_keys}"
                )
                if "unit" in check_prop_units[prop]:
                    assert "units" in traj[prop].attrs, (
                        f"Property {prop} has no unit set in resulting trajectory"
                    )
                    required_unit = check_prop_units[prop]["unit"]
                    actual_unit = traj[prop].attrs["units"]
                    assert actual_unit == required_unit, (
                        f"Property {prop} has unit {actual_unit} instead of required unit {required_unit}"
                    )
        for varset in [traj.data_vars, traj.coords]:
            for var in varset:
                if "unitdim" in traj[var].attrs:
                    assert "units" in traj[var].attrs, (
                        f"Variable {var} has property `unitdim` but no `units` set."
                    )
                    unit_dim = traj[var].attrs["unitdim"]
                    actual_unit = traj[var].attrs["units"]
                    required_unit = standard_shnitsel_units[unit_dim]
                    assert required_unit == "1" or actual_unit == required_unit, (
                        f"Variable {var} has unit {actual_unit} instead of required unit {required_unit}."
                    )

        return True
