from functools import lru_cache
import logging
import re

import numpy as np
import xarray as xr

from shnitsel.io.shared.variable_flagging import (
    is_variable_assigned,
    mark_variable_assigned,
)
from typing import List, Optional

higher_order_state_pattern_re = re.compile(
    r"^S\[(?P<multiplicity>\d+)\](?P<state_mult_index>\d+)::(?P<magnetic_number>-?\d+)$"
)
triplet_state_pattern_re = re.compile(r"^T(?P<state_mult_index>\d+)(?P<suffix>[-\+]?)$")
doublet_state_pattern_re = re.compile(r"^D(?P<state_mult_index>\d+)(?P<suffix>[-\+])$")
singlet_state_pattern_re = re.compile(r"^S(?P<state_mult_index>\d+)$")


def default_state_type_assigner(dataset: xr.Dataset) -> xr.Dataset:
    """Function to assign default state types to states independent of the format.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset to assign the states to

    Returns
    -------
    xr.Dataset
        The dataset after the assignment
    """
    # If state types have already been set, do not touch them
    if is_variable_assigned(dataset.state_types):
        return dataset

    # Try and extract the state types from the number of different states
    nsinglets = dataset.attrs.get("num_singlets", -1)
    ndoublets = dataset.attrs.get("num_doublets", -1)
    ntriplets = dataset.attrs.get("num_triplets", -1)

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

        if ndoublets == 0 and ntriplets == 0:
            # Special case if we only have singlets:

            if not is_variable_assigned(dataset.state_magnetic_number):
                dataset.state_magnetic_number[:nsinglets] = 0.0
                mark_variable_assigned(dataset.state_magnetic_number)
            if not is_variable_assigned(dataset.state_degeneracy_group):
                dataset.state_degeneracy_group[:nsinglets] = dataset.state[:nsinglets]
                mark_variable_assigned(dataset.state_degeneracy_group)
    return dataset


def default_state_name_assigner(dataset: xr.Dataset) -> xr.Dataset:
    """Function to assign default state names to states.

    State names for Singlets are (S0, S1, S2, S3, S4...) higher-order multiplicities start with index 1 (no suffix for momentum in default naming due to lack of information).
    Prefixes for singlets, doublets and triplets are `S`, `D`, `T`. Higher-order states are not considered.

    Parameters
    ----------
    dataset : xr.Dataset
        The dataset to assign the states to

    Returns
    -------
    xr.Dataset
        The dataset after the assignment
    """
    # Do not touch previously set names
    if is_variable_assigned(dataset.state_names):
        logging.info("State names already assigned")
        return dataset

    if is_variable_assigned(dataset.state_types):
        counters = np.array([0, 1, 1], dtype=int)
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
        nsinglets = dataset.attrs.get("num_singlets",-1)
        ndoublets = dataset.attrs.get("num_doublets",-1)
        ntriplets = dataset.attrs.get("num_triplets",-1)

        if nsinglets >= 0 and ndoublets >= 0 and ntriplets >= 0:
            logging.warning(
                "We made a best-effort guess for the names of the individual states. "
                "Please provide a list of state names or a function ot assign the state names to have the correct values assigned."
            )
            new_name_values = dataset.state_names
            if nsinglets > 0:
                new_name_values[:nsinglets] = [f"S{i}" for i in range(nsinglets)]
            if ndoublets > 0:
                # We skip label 0 for higher-order states
                new_name_values[nsinglets : nsinglets + ndoublets] = [
                    f"D{i + 1}" for i in range(ndoublets)
                ]
            if ntriplets > 0:
                # We skip label 0 for higher-order states
                new_name_values[nsinglets + ndoublets :] = [
                    f"T{i + 1}" for i in range(ntriplets)
                ]
            dataset = dataset.assign_coords(
                {"state_names": ("state", new_name_values, dataset.state_names.attrs)}
            )

            mark_variable_assigned(dataset.state_names)

    return dataset


def set_sharc_state_type_and_name_defaults(
    dataset: xr.Dataset,
    multiplicity_counts: List[int] | int,
    multiplicity_charges: Optional[List | int | float] = None,
) -> xr.Dataset:
    """Apply default sharc naming scheme to a dataset and set the state order appropriately. This is more specific than the general naming convention. Enumerates spin numbers on top of common plain enumeration of dublets and triplets.

    Can also be used to set the charges per state.
    State names for Singlets are (S0, S1, S2, S3, S4...) higher-order multiplicities start with index 1 and have a suffix depending on the angular momentum (+,none, -).
    Prefixes for singlets, doublets and triplets are `S`, `D`, `T`. Higher-order states are named with the pattern `S[<multiplicity>]<label index in multiplicity>::<angular momentum index>`.

    Parameters
    ----------
    dataset : xr.Dataset
        The input dataset to set the states on
    multiplicity_counts : List[int] | int
        The list of amount of states of different multiplicities or the number of singlet states
    multiplicity_charges : List | int | float, optional
        The list of charges of different states or the charge to apply to all states.
            If not set, no charge will be set for all states

    Returns
    -------
    xr.Dataset
        The dataset with state types, names and charges applied.
    """
    if not isinstance(multiplicity_counts, list):
        multiplicity_counts = [multiplicity_counts]

    max_mult = len(multiplicity_counts)

    if multiplicity_charges is None:
        multiplicity_charges = [0] * max_mult
    elif isinstance(multiplicity_charges, list):
        len_charges = len(multiplicity_charges)
        if len_charges != max_mult:
            logging.warning(
                f"Length of charge and multiplicity arrays differ: {max_mult} vs. {len(multiplicity_charges)}. Padding with zeroes."
            )
            if max_mult > len_charges:
                multiplicity_charges = multiplicity_charges + [0] * (
                    max_mult - len_charges
                )
    else:
        multiplicity_charges = [multiplicity_charges] * max_mult

    curr_index = 0
    degeneracy_index = 0

    if max_mult >= 1:
        for i in range(multiplicity_counts[0]):
            dataset.state_types[curr_index] = 1
            dataset.state_names[curr_index] = f"S{i}"
            dataset.state_charges[curr_index] = multiplicity_charges[0]
            dataset.state_magnetic_number[curr_index] = 0.0  # For singlets
            dataset.state_degeneracy_group[curr_index] = degeneracy_index + i
            curr_index += 1
        degeneracy_index += multiplicity_counts[0]

    if max_mult >= 2:
        curr_mult = 2
        suffix = ["-", "+"]
        magnetic_nums = [-0.5, 0.5]
        charge = multiplicity_charges[1]
        for m in range(0, curr_mult):
            # We skip label 0 for higher-order states
            for i in range(1, 1 + multiplicity_counts[1]):
                dataset.state_types[curr_index] = curr_mult
                dataset.state_names[curr_index] = f"D{i}{suffix[m]}"
                dataset.state_charges[curr_index] = charge
                dataset.state_magnetic_number[curr_index] = magnetic_nums[m]
                dataset.state_degeneracy_group[curr_index] = degeneracy_index + i  #
                curr_index += 1

        degeneracy_index += multiplicity_counts[1]

    if max_mult >= 3:
        curr_mult = 3
        suffix = ["-", "", "+"]
        magnetic_nums = [-1.0, 0.0, 1.0]
        charge = multiplicity_charges[2]
        for m in range(0, curr_mult):
            # We skip label 0 for higher-order states
            for i in range(1, 1 + multiplicity_counts[2]):
                dataset.state_types[curr_index] = curr_mult
                dataset.state_names[curr_index] = f"T{i}{suffix[m]}"
                dataset.state_charges[curr_index] = charge
                dataset.state_magnetic_number[curr_index] = magnetic_nums[m]
                dataset.state_degeneracy_group[curr_index] = degeneracy_index + i  #
                curr_index += 1
        degeneracy_index += multiplicity_counts[2]

    if max_mult > 3:
        for curr_mult in range(4, max_mult + 1):
            charge = multiplicity_charges[curr_mult - 1]
            center_mag_offset = float(curr_mult) / 2.0
            for m in range(0, curr_mult):
                mag_num = float(m) - center_mag_offset

                # We skip label 0 for higher-order states
                for i in range(1, 1 + multiplicity_counts[curr_mult - 1]):
                    dataset.state_types[curr_index] = curr_mult
                    dataset.state_names[curr_index] = (
                        f"S[{curr_mult}]{i}::{mag_num:.1f}"
                    )
                    dataset.state_charges[curr_index] = charge
                    dataset.state_magnetic_number[curr_index] = mag_num
                    dataset.state_degeneracy_group[curr_index] = degeneracy_index + i  #
                    curr_index += 1
            degeneracy_index += multiplicity_counts[curr_mult - 1]

    mark_variable_assigned(dataset.state_types)
    mark_variable_assigned(dataset.state_names)
    mark_variable_assigned(dataset.state_charges)
    mark_variable_assigned(dataset.state_magnetic_number)
    mark_variable_assigned(dataset.state_degeneracy_group)

    return dataset


@lru_cache
def state_name_to_tex_label(statename: str) -> str:
    """Function to translate default state naming conventions into a general latex-subcscrip/-superscript label.

    Parameters
    ----------
    statename : str
        Statename as per Shnitsel default convention.

    Returns
    -------
    str
        A LaTeX representation of the state label
    """
    singlet_match = singlet_state_pattern_re.match(statename)
    if singlet_match:
        singlet_index = singlet_match.group("state_mult_index")
        return r"S_{" + singlet_index + r"}"

    triplet_match = triplet_state_pattern_re.match(statename)
    if triplet_match:
        triplet_index = triplet_match.group("state_mult_index")
        suffix = triplet_match.group("suffix")
        return r"T_{" + triplet_index + r"}" + ("^{" + suffix + "}" if suffix else "")

    doublet_match = doublet_state_pattern_re.match(statename)
    if doublet_match:
        doublet_index = doublet_match.group("state_mult_index")
        suffix = doublet_match.group("suffix")
        return r"D_{" + doublet_index + r"}" + ("^{" + suffix + "}" if suffix else "")

    higher_match = higher_order_state_pattern_re.match(statename)
    if higher_match:
        higher_mult = higher_match.group("multiplicity")
        higher_index = higher_match.group("state_mult_index")
        higher_m = higher_match.group("magnetic_number")
        return r"S_{m=" + higher_mult + ",i=" + higher_index + ",j=" + higher_m + r"}"

    logging.info(f"Failed to translate state name to label for state {statename}")
    return statename
