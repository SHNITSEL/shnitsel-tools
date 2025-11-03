import logging
from typing import Any, List

from shnitsel.data.TrajectoryFormat import Trajectory
import xarray as xr
from shnitsel.units import standard_shnitsel_units
from shnitsel.units.definitions import unit_dimensions


def verify_trajectory_format(obj: Any) -> bool:
    return is_permitted_traj_result(obj) and has_required_properties(obj)


def is_permitted_traj_result(obj: Any) -> bool:
    logging.debug(f"Obj has type: {type(obj)}")
    return (
        isinstance(obj, xr.Dataset)
        or isinstance(obj, Trajectory)
        or (isinstance(obj, list) and (len(obj) == 0 or isinstance(obj[0], Trajectory)))
    )


def has_required_properties(traj: List[Trajectory] | Trajectory) -> bool:
    res = True
    if isinstance(traj, list):
        for i_traj in traj:
            res = res and has_required_properties(i_traj)
        return res
    else:

        print(traj.variables.keys())
        print(traj["atXYZ"].attrs.keys())

        check_prop_units = {
            "atXYZ": {"unit": standard_shnitsel_units[unit_dimensions.length]},
            "forces": {"unit": standard_shnitsel_units[unit_dimensions.force]},
            "energy": {"unit": standard_shnitsel_units[unit_dimensions.energy]},
            "dip_trans": {"unit": standard_shnitsel_units[unit_dimensions.dipole]},
            "dip_perm": {"unit": standard_shnitsel_units[unit_dimensions.dipole]},
            "states": {},
        }

        for prop in check_prop_units:
            assert (
                prop in traj.variables.keys()
            ), f"Property {prop} is missing in resulting trajectory"
            if "unit" in check_prop_units[prop]:
                assert (
                    "units" in traj["atXYZ"].attrs
                ), f"Property {prop} has no unit set in resulting trajectory"
                required_unit = check_prop_units[prop]["unit"]
                actual_unit = traj["atXYZ"].attrs["units"]
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
