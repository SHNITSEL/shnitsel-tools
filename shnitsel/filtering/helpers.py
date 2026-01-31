from typing import Sequence
import rdkit.Chem as rc
import xarray as xr
from shnitsel.bridges import construct_default_mol
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.filtering.state_selection import StateSelection, StateSelectionDescriptor
from .structure_selection import (
    FeatureLevelOptions,
    StructureSelection,
    StructureSelectionDescriptor,
)


def _get_default_structure_selection(
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    mol: rc.Mol | None = None,
    atXYZ_source: xr.Dataset | xr.DataArray | ShnitselDataset | None = None,
    charge_info: int | None = None,
    default_levels: Sequence[FeatureLevelOptions] = ['atoms', 'bonds'],
) -> StructureSelection:
    """Get a default structure selection object from any accessible data if possible.

    Parameters
    ----------
    structure_selection : StructureSelection | StructureSelectionDescriptor | None, optional
        A potential already provided structure/feature selection.
        Alternatively a description of the selection to be applied to a dataset as, e.g.,
        with a SMARTS selection, a selection of indices to consider, etc. Defaults to None.
        If provided, will be returned back.
    mol : rc.Mol | None, optional
        An optional instance of an RDKit molecule. Used to construct a StructureSelection instance if `structure_selection` is None. Defaults to None.
    atXYZ : xr.Dataset | xr.DataArray | ShnitselDataset | None, optional
        An xr.DataArray holding positional data of a molecule or a different source for the positional data. Only the first frame will be used to create a default mol if `structure_selection` and `mol` were not provided. Defaults to None.
    default_levels : Sequence[FeatureLevelOptions], optional
        the desired default levels included in the selection that may be recreated if none was provided. Defaults to 'atoms' and 'bonds'.

    Raises
    ------
    ValueError
        Not enough data provided to construct a default selection from.

    Returns
    -------
    StructureSelection
        The initialized default structure selection.
    """

    if structure_selection is not None:
        if isinstance(structure_selection, StructureSelection):
            return structure_selection
        else:
            # try and construct structure selection:
            try:
                if mol is not None:
                    sel = StructureSelection.init_from_mol(mol).select_all()
                else:
                    sel = StructureSelection.init_from_dataset(
                        atXYZ_source
                    ).select_all()

                tmp_res = sel.derive_other_from_descriptor(structure_selection)
            except Exception as e:
                raise e

            if tmp_res is None:
                raise ValueError(
                    "Could not construct StructureSelection from provided rdkit.Mol or dataset combined and StructureSelection descriptor"
                )

            return tmp_res

    if mol is not None and isinstance(mol, rc.Mol):
        return StructureSelection.init_from_mol(mol, default_selection=default_levels)
    elif mol is None and atXYZ_source is not None:
        mol = construct_default_mol(atXYZ_source, charge=charge_info)
    elif mol is None:
        raise ValueError(
            "You did not provide sufficient data to construct a default feature selection. Please provide your own StructureSelection object."
        )

    return StructureSelection.init_from_mol(mol, default_selection=default_levels)


def _get_default_state_selection(
    state_selection: StateSelection | StateSelectionDescriptor | None = None,
    state_source: xr.Dataset | xr.DataArray | ShnitselDataset | None = None,
) -> StateSelection:
    """Get a default state selection object from any accessible data if possible.

    Parameters
    ----------
    state_selection : StateSelection | StateSelectionDescriptor | None, optional
        A potential already provided state selection.
        Alternatively a description of the selection to be applied to a dataset as, e.g.,
        with a description of state types selection, a selection of indices to consider, etc. Defaults to None.
        If provided, will be returned back.
    state_source : xr.Dataset | xr.DataArray | ShnitselDataset | None, optional
        An xr.DataArray or a kind of dataset that may hold state information of a system. Defaults to None.

    Raises
    ------
    ValueError
        Not enough data provided to construct a default selection from.

    Returns
    -------
    StateSelection
        The initialized default structure selection.
    """

    if state_selection is not None:
        if isinstance(state_selection, StateSelection):
            return state_selection
        else:
            # try and construct structure selection:
            try:
                sel = StateSelection.init_from_dataset(state_source)

                tmp_res = sel.select(state_selection)
            except:
                try:
                    tmp_res = StateSelection.init_from_descriptor(state_selection)
                except:
                    raise

            if tmp_res is None:
                raise ValueError(
                    "Could not construct StateSelection from provided dataset combined with the StateSelection descriptor"
                )

            return tmp_res

    return StateSelection.init_from_dataset(state_source)
