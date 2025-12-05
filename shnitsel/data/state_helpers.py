import logging

import numpy as np
import xarray as xr

from shnitsel.io.shared.variable_flagging import (
    is_variable_assigned,
    mark_variable_assigned,
)


def default_state_type_assigner(dataset: xr.Dataset) -> xr.Dataset:
    """Function to assign default state types to states.

    Args:
        dataset (xr.Dataset): The dataset to assign the states to

    Returns:
        xr.Dataset: The dataset after the assignment
    """
    # If state types have already been set, do not touch them
    if is_variable_assigned(dataset.state_types):
        return dataset

    # Try and extract the state types from the number of different states
    nsinglets = dataset.attrs["num_singlets"]
    ndoublets = dataset.attrs["num_doublets"]
    ntriplets = dataset.attrs["num_triplets"]

    if nsinglets >= 0 and ndoublets >= 0 and ntriplets >= 0:
        # logging.debug(f"S/D/T = {nsinglets}/{ndoublets}/{ntriplets}")
        logging.warning(
            "We made a best-effort guess for the types/multiplicities of the individual states. "
            "Please provide a list of state types or a function to assign the state types to have the correct values assigned."
        )
        if nsinglets > 0:
            dataset.state_types[:nsinglets] = 1
        if ndoublets > 0:
            dataset.state_types[nsinglets : nsinglets + ndoublets] = 2
        if ntriplets > 0:
            dataset.state_types[nsinglets + ndoublets :] = 3
        keep_attr = dataset.state_types.attrs

        dataset = dataset.reindex({"state_types": dataset.state_types.values})
        dataset.state_types.attrs.update(keep_attr)

        mark_variable_assigned(dataset.state_types)
    return dataset


def default_state_name_assigner(dataset: xr.Dataset) -> xr.Dataset:
    """Function to assign default state names to states.

    Args:
        dataset (xr.Dataset): The dataset to assign the states to

    Returns:
        xr.Dataset: The dataset after the assignment
    """
    # Do not touch previously set names
    if is_variable_assigned(dataset.state_names):
        logging.info("State names already assigned")
        return dataset

    if is_variable_assigned(dataset.state_types):
        counters = np.array([0, 0, 0], dtype=int)
        type_prefix = np.array(["S", "D", "T"])
        type_values = dataset.state_types.values

        res_names = []
        for i in range(len(type_values)):
            type_index = int(round(type_values[i]))
            assert type_index >= 1 and type_index <= 3, (
                f"Found invalid state multiplicity for default naming: {type_index} (must be 1,2 or 3)"
            )
            # logging.debug(
            #     f"{i}, {type_index}, {type_prefix[type_index - 1]}, {counters[type_index - 1]}"
            # )
            res_names.append(
                type_prefix[type_index - 1] + f"{counters[type_index - 1]:d}"
            )
            counters[type_index - 1] += 1

        # logging.info(
        #    "State names assigned based on types: {type_values} -> {res_names}"
        # )
        dataset = dataset.assign_coords(
            {"state_names": ("state", res_names, dataset.state_names.attrs)}
        )

        mark_variable_assigned(dataset.state_names)
        # logging.debug(f"Default name set on type basis: {repr(dataset)}")
    else:
        nsinglets = dataset.attrs["num_singlets"]
        ndoublets = dataset.attrs["num_doublets"]
        ntriplets = dataset.attrs["num_triplets"]

        if nsinglets >= 0 and ndoublets >= 0 and ntriplets >= 0:
            logging.warning(
                "We made a best-effort guess for the names of the individual states. "
                "Please provide a list of state names or a function ot assign the state names to have the correct values assigned."
            )
            new_name_values = dataset.state_names
            if nsinglets > 0:
                new_name_values[:nsinglets] = [f"S{i}" for i in range(nsinglets)]
            if ndoublets > 0:
                new_name_values[nsinglets : nsinglets + ndoublets] = [
                    f"D{i}" for i in range(ndoublets)
                ]
            if ntriplets > 0:
                new_name_values[nsinglets + ndoublets :] = [
                    f"T{i}" for i in range(ntriplets)
                ]
            dataset = dataset.assign_coords(
                {"state_names": ("state", new_name_values, dataset.state_names.attrs)}
            )

            mark_variable_assigned(dataset.state_names)

    return dataset
