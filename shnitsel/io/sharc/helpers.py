import logging
import pathlib
from typing import Any, Callable, Dict, List, Optional
import numpy as np
import xarray as xr

from shnitsel.io.shared.variable_flagging import mark_variable_assigned


def set_sharc_state_type_and_name_defaults(
    dataset: xr.Dataset,
    multiplicity_counts: List[int] | int,
    multiplicity_charges: Optional[List | int | float] = None,
) -> xr.Dataset:
    """Apply default sharc naming scheme to a dataset and set the state order appropriately.

    Can also be used to set the charges per state.

    Args:
        dataset (xr.Dataset): The input dataset to set the states on
        multiplicity_counts (List[int] | int): The list of amount of states of different multiplicities or the number of singlet states
        multiplicity_charges (List | int | float, optional): The list of charges of different states or the charge to apply to all states.
                If not set, no charge will be set for all states

    Returns:
        xr.Dataset: The dataset with state types, names and charges applied.
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

    if max_mult >= 1:
        for i in range(multiplicity_counts[0]):
            dataset.state_types[curr_index] = 1
            dataset.state_names[curr_index] = f"S{i}"
            dataset.state_charges[curr_index] = multiplicity_charges[0]
            curr_index += 1
    if max_mult >= 2:
        curr_mult = 2
        suffix = ["-", "+"]
        charge = multiplicity_charges[1]
        for i in range(multiplicity_counts[1]):
            for m in range(curr_mult):
                dataset.state_types[curr_index] = curr_mult
                dataset.state_names[curr_index] = f"D{i}{suffix[m]}"
                dataset.state_charges[curr_index] = charge
                curr_index += 1
    if max_mult >= 3:
        curr_mult = 3
        suffix = ["-", "", "+"]
        charge = multiplicity_charges[2]
        for i in range(multiplicity_counts[2]):
            for m in range(curr_mult):
                dataset.state_types[curr_index] = curr_mult
                dataset.state_names[curr_index] = f"T{i}{suffix[m]}"
                dataset.state_charges[curr_index] = charge
                curr_index += 1
    if max_mult > 3:
        for curr_mult in range(4, max_mult + 1):
            charge = multiplicity_charges[curr_mult - 1]
            lower_mult = -int(curr_mult // 2)
            upper_mult = int(curr_mult // 2)
            for i in range(multiplicity_counts[curr_mult - 1]):
                for m in range(lower_mult, upper_mult + 1):
                    dataset.state_types[curr_index] = curr_mult
                    dataset.state_names[curr_index] = f"S[{curr_mult}]{i}::{m}"
                    dataset.state_charges[curr_index] = charge
                    curr_index += 1

    mark_variable_assigned(dataset.state_types)
    mark_variable_assigned(dataset.state_names)
    mark_variable_assigned(dataset.state_charges)

    return dataset


def read_molcas_qm_info(dataset: xr.Dataset, qm_path: pathlib.Path) -> Dict[str, Any]:
    """An attempt to read QM-interface info from SHARC 2.x and 3.0 MOLCAS interface data.

    Used to determine the charge of the molecule.

    Args:
        dataset (xr.Dataset): Input set to access certain features of the setup
        qm_path (pathlib.Path): Path to the QM/ folder in the simulation directory

    Returns:
        Dict[str, Any]: Resulting settings.
    """
    resources_file = qm_path / "MOLCAS.resources"
    template_file = qm_path / "MOLCAS.template"

    res = {}
    if resources_file.is_file():
        resource_res = {}
        with open(resources_file) as mol_res:
            lines = mol_res.readlines()

            for line in lines:
                line = line.strip()
                if len(line) == 0 or line.startswith("#"):
                    continue
                parts = [x.strip() for x in line.split(maxsplit=1)]

                if len(parts) == 1:
                    resource_res[parts[0]] = True
                else:
                    resource_res[parts[0]] = parts[1]
        res["MOLCAS.resources"] = resource_res

    if template_file.is_file():
        template_res = {}
        with open(resources_file) as mol_res:
            lines = mol_res.readlines()

            for line in lines:
                line = line.strip()
                if len(line) == 0 or line.startswith("#"):
                    continue
                parts = [x.strip() for x in line.split(maxsplit=1)]

                if len(parts) == 1:
                    template_res[parts[0]] = True
                else:
                    template_res[parts[0]] = parts[1]
        res["MOLCAS.template"] = template_res

        if "basis" in template_res:
            res["theory_basis"] = template_res["basis"]
        if "method" in template_res:
            res["est_level"] = template_res["method"]

        num_active_electrons = 0
        num_inactive_orbitals = 0
        if "nactel" in template_res:
            num_active_electrons = int(template_res["nactel"])
        if "inactive" in template_res:
            num_inactive_orbitals = int(template_res["inactive"])

        total_simulated_electrons = num_active_electrons + 2 * num_inactive_orbitals

        if total_simulated_electrons > 0:
            total_required_electrons = np.sum(dataset.atNums.values)

            total_charge = total_required_electrons - total_simulated_electrons

            res["charge"] = total_charge

    return res


def read_columbus_qm_info(dataset: xr.Dataset, qm_path: pathlib.Path) -> Dict[str, Any]:
    """An attempt to read QM-interface info from SHARC 2.x and 3.0 COLUMBUS interface data.

    Used to determine the charge of the molecule.

    Args:
        dataset (xr.Dataset): Input set to access certain features of the setup
        qm_path (pathlib.Path): Path to the QM/ folder in the simulation directory

    Returns:
        Dict[str, Any]: Resulting settings.
    """
    resources_file = qm_path / "COLUMBUS.resources"
    template_dir: Optional[pathlib.Path] = None

    res = {}
    template_state_dirs = {}
    if resources_file.is_file():
        resource_res = {}
        with open(resources_file) as mol_res:
            lines = mol_res.readlines()

            for line in lines:
                line = line.strip()
                if len(line) == 0 or line.startswith("#"):
                    continue
                parts = [x.strip() for x in line.split(maxsplit=1)]

                if len(parts) == 1:
                    resource_res[parts[0]] = True
                else:
                    if parts[0] == "DIR":
                        s = parts[1].split()
                        state_num, state_dir = int(s[0]), s[1]
                        template_state_dirs[state_num] = state_dir

                    resource_res[parts[0]] = parts[1]
        res["COLUMBUS.resources"] = lines

        # The template directory must be provided in the resources file
        template_dir = pathlib.Path(resource_res["template"])

    if template_dir is not None and template_dir.is_file():
        total_required_electrons = np.sum(dataset.atNums.values)

        state_charges = []

        for multiplicity in sorted([int(x) for x in template_state_dirs.keys()]):
            multiplicity_dir = template_dir / str(
                template_state_dirs[str(multiplicity)]
            )

            total_charge = None

            for path in multiplicity_dir.glob("cidtrin*"):
                with open(path) as cidtrin_f:
                    for line in cidtrin_f.read():
                        line = line.strip()

                        if line.find("total number of electrons") != -1:
                            total_num_electrons = int(line.split("/")[0])

                            total_charge = (
                                total_required_electrons - total_num_electrons
                            )
                            break
                if total_charge is not None:
                    break

            state_charges.append(0 if total_charge is None else total_charge)

        res["charge"] = state_charges
    return res


INTERFACE_READERS: Dict[str, Callable[[xr.Dataset, pathlib.Path], Dict[str, Any]]] = {
    "molcas": read_molcas_qm_info,
    "columbus": read_columbus_qm_info,
}
