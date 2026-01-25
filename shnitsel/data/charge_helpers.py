import logging
from typing import Any
import xarray as xr

from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.node import TreeNode
from shnitsel.bridges import construct_default_mol


def guess_molecular_charge(
    atXYZ_source: xr.DataArray
    | xr.Dataset
    | ShnitselDataset
    | TreeNode[Any, xr.DataArray | xr.Dataset | ShnitselDataset],
    max_range_abs: int = 5,
) -> int:
    """Helper function to guess the charge of a system configuration.

    Will attempt to make RDKit generate a molecule with various molecular charges starting with the smallest possible
    absolute values (0) and working upwards to the maximum absolute value of `max_range_abs`.
    It will try both negative and positive values.

    Parameters
    ----------
    atXYZ_source : xr.DataArray | xr.Dataset | ShnitselDataset | TreeNode[Any, xr.DataArray  |  xr.Dataset  |  ShnitselDataset]
        The source for positional data.
        Can either be only the positional data as a DataArray or a full (Shnitsel) dataset.
    max_range_abs : int

    Raises
    ------
    RuntimeError
        If no charge that makes the system self-consistent could be found in the range.

    Returns
    -------
    int
        The charge with the smallest absolute value within the specified range for which `rdkit` returned a self-consistent
        molecular state.
    """

    try:
        initial_attempt = construct_default_mol(
            atXYZ_source, charge=0, silent_mode=True
        )
        logging.info(
            "RDKit returned a self-consistent molecular structure for charge=0. Assuming molecule to be uncharged."
        )

        # Set mol as cached on the dataset
        if isinstance(atXYZ_source, xr.Dataset):
            atXYZ_source.attrs["__mol"] = initial_attempt
        elif isinstance(atXYZ_source, ShnitselDataset):
            atXYZ_source.dataset.attrs["__mol"] = initial_attempt

        return 0
    except:
        logging.info(
            "RDKit was not able to construct a molecular structure without a charge.\n Trying to guess molecular charge."
        )

        for abs_val in range(1, abs(max_range_abs) + 1, 1):
            for sign in [1, -1]:
                charge_guess = sign * abs_val
                try:
                    charged_attempt = construct_default_mol(
                        atXYZ_source, charge=charge_guess, silent_mode=True
                    )
                    logging.info(
                        "RDKit returned a self-consistent molecular structure for charge = %d. Assuming that molecular charge. Use `.set_charge(<charge>)` to set the charge if the value is not what you expected.",
                        charge_guess,
                    )

                    # Set mol as cached on the dataset
                    if isinstance(atXYZ_source, xr.Dataset):
                        atXYZ_source.attrs["__mol"] = charged_attempt
                    elif isinstance(atXYZ_source, ShnitselDataset):
                        atXYZ_source.dataset.attrs["__mol"] = charged_attempt

                    return charge_guess
                except Exception as e:
                    logging.debug(
                        "Caught exception during charge guess with charge=%d : %s",
                        charge_guess,
                        e,
                    )
                    continue

        raise RuntimeError(
            "No molecular charge up to a maximum absolute value of %d yielded a self-consistent molecular structure."
            % (abs(max_range_abs))
        )
