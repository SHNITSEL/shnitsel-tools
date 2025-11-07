import logging
import xarray as xr

from shnitsel.data.TrajectoryFormat import Trajectory, wrap_trajectory
from shnitsel.io.helpers import LoadingParameters
from shnitsel.units.conversion import convert_all_units_to_shnitsel_defaults


def finalize_loaded_trajectory(dataset: xr.Dataset | None, loading_parameters: LoadingParameters | None) -> Trajectory | None:
    if dataset is not None:
        # TODO: FIXME: use loading_parameters to configure state names
        # TODO: FIXME: use loading_parameters to configure state names
        unset_vars = []
        for var in dataset.variables:
            if "__assigned" in dataset[var].attrs:
                # Remove tags
                del dataset[var].attrs["__assigned"]
            else:
                unset_vars.append(var)

        dataset.drop_vars(unset_vars)

        dataset = set_state_defaults(dataset, loading_parameters)

        return wrap_trajectory(convert_all_units_to_shnitsel_defaults(dataset))

    return None


def set_state_defaults(dataset: xr.Dataset, loading_parameters: LoadingParameters | None) -> xr.Dataset:

    # TODO: FIXME: apply configured names from loading_parameters
    nsinglets = dataset.attrs["num_singlets"]
    ndoublets = dataset.attrs["num_doublets"]
    ntriplets = dataset.attrs["num_triplets"]

    if nsinglets >= 0 and ndoublets >= 0 and ntriplets >= 0:
        logging.warning("We made a best-effort guess for the names and types/multiplicities of the individual states. "
                        "Please provide a list of state types or a function ot assign the state types to have the correct values assigned.")
        if nsinglets > 0:
            dataset.state_types[:nsinglets] = 1
            dataset.state_names[:nsinglets] = [
                f"S{i}" for i in range(nsinglets)]
        if ndoublets > 0:
            dataset.state_types[nsinglets: nsinglets + 2 * ndoublets] = 2
            dataset.state_names[nsinglets: nsinglets + 2 * ndoublets] = [
                f"D{i}" for i in range(2 * ndoublets)
            ]
        if ntriplets > 0:
            dataset.state_types[nsinglets + 2 * ndoublets:] = 3
            dataset.state_names[nsinglets + 2 *
                                ndoublets:] = [f"T{i}" for i in range(3 * ntriplets)]
    else:
        logging.error("Could not determine state multiplicities and names")

    return dataset
