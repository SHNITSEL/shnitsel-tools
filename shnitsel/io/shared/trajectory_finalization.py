import xarray as xr

from shnitsel.data.TrajectoryFormat import Trajectory, wrap_trajectory
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults


def finalize_loaded_trajectory(dataset: xr.Dataset | None) -> Trajectory | None:
    if dataset is not None:
        unset_vars = []
        for var in dataset.variables:
            if "__assigned" in dataset[var].attrs:
                # Remove tags
                del dataset[var].attrs["__assigned"]
            else:
                unset_vars.append(var)

        dataset.drop_vars(unset_vars)

        return wrap_trajectory(convert_all_units_to_shnitsel_defaults(dataset))
    
    return None
