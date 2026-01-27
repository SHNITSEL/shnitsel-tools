from typing import Any, overload
from shnitsel._contracts import needs
import xarray as xr

from shnitsel.core.typedefs import AtXYZ
from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.node import TreeNode
from shnitsel.filtering.structure_selection import (
    BondDescriptor,
    FeatureTypeLabel,
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.geo.geocalc_.distances import get_distances
from shnitsel.geo.geocalc_.helpers import (
    _assign_descriptor_coords,
    _empty_descriptor_results,
)
from shnitsel.filtering.helpers import _get_default_structure_selection


@overload
def get_max_chromophor_BLA(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray],
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    SMARTS: str | None = None,
    num_double_bonds: int | None = None,
    allowed_chain_elements: str = "#6,#7,#8,#15,#16",
    max_considered_BLA_double_bonds: int = 50,
) -> TreeNode[Any, xr.DataArray]: ...


@overload
def get_max_chromophor_BLA(
    atXYZ_source: Trajectory | Frames | xr.Dataset | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    SMARTS: str | None = None,
    num_double_bonds: int | None = None,
    allowed_chain_elements: str = "#6,#7,#8,#15,#16",
    max_considered_BLA_double_bonds: int = 50,
) -> xr.DataArray: ...


@needs(dims={'atom', 'direction'})
def get_max_chromophor_BLA(
    atXYZ_source: TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray]
    | Trajectory
    | Frames
    | xr.Dataset
    | xr.DataArray,
    structure_selection: StructureSelection
    | StructureSelectionDescriptor
    | None = None,
    SMARTS: str | None = None,
    num_double_bonds: int | None = None,
    allowed_chain_elements: str = "#6,#7,#8,#15,#16",
    max_considered_BLA_double_bonds: int = 50,
) -> TreeNode[Any, xr.DataArray] | xr.DataArray:
    """Calculate bond length alternation value (BLA) for the maximum chromophor
    in the provided `structure_selection` or the maximum chromophor in the structure
    represented by `atXYZ`.

    Parameters
    ----------
    atXYZ_source : TreeNode[Any, Trajectory | Frames | xr.Dataset | xr.DataArray] | Trajectory | Frames | xr.Dataset | xr.DataArray
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions ``atom`` and
        ``direction`` or another source of positional data like a trajectory, a frameset,
        a dataset representing either of those or a tree structure holding such data.
    structure_selection | StructureSelectionDescriptor, optional
        Object encapsulating feature selection on the structure whose positional information is provided in `atXYZ`.
        If this argument is omitted altogether, a default selection for all bonds and atoms within the structure is created.
    SMARTS : str, optional
        SMARTS string to match for the maximum chromophor.
    num_double_bonds : int, optional
        The specified number of double bonds for the maximum chromophor.
    allowed_chain_elements : str, default="#6,#7,#8,#15,#16"
        SMARTS atomic specification, i.e. comma-separated list of element descriptors (default: C,N,O,P,S represented as '#6,#7,#8,#15,#16').
    max_considered_BLA_double_bonds: int, optional
        Maximum number of double bonds in a BLA chromophor if automatic maximum size detection is performed. Defaults to 50.

    Returns
    -------
        An :py:class:`xarray.DataArray` of the BLA for the maximum-length chromophor (alternating double bonds)

    Raises
    -------
    ValueError
        If the maximum chromophor within the provided selection or the entire molecule is not unique.

    """

    if isinstance(atXYZ_source, TreeNode):
        return atXYZ_source.map_data(
            lambda x: get_max_chromophor_BLA(
                x,
                structure_selection=structure_selection,
                SMARTS=SMARTS,
                num_double_bonds=num_double_bonds,
                allowed_chain_elements=allowed_chain_elements,
                max_considered_BLA_double_bonds=max_considered_BLA_double_bonds,
            ),
            keep_empty_branches=True,
            dtype=xr.DataArray,
        )

    position_data: xr.DataArray
    charge_info: int | None
    if isinstance(atXYZ_source, xr.DataArray):
        position_data = atXYZ_source
        charge_info = None
    else:
        wrapped_ds = wrap_dataset(atXYZ_source, (Trajectory | Frames))
        position_data = wrapped_ds.atXYZ
        charge_info = int(wrapped_ds.charge)

    structure_selection = _get_default_structure_selection(
        structure_selection,
        atXYZ_source=position_data,
        default_levels=['atoms', 'bonds'],
        charge_info=charge_info,
    )

    # Get the maximum BLA chromophor
    BLA_selection = structure_selection.select_BLA_chromophor(
        BLA_smarts=SMARTS,
        num_double_bonds=num_double_bonds,
        allowed_chain_elements=allowed_chain_elements,
        max_considered_BLA_double_bonds=max_considered_BLA_double_bonds,
    )

    if len(BLA_selection.bonds_selected) == 0:
        return _empty_descriptor_results(position_data)

    bond_lengths = get_distances(position_data, BLA_selection)

    single_idxs: list[BondDescriptor] = []
    double_idxs: list[BondDescriptor] = []
    for bond_id in BLA_selection.bonds_selected:
        bond_type = BLA_selection.get_bond_type(bond_id)
        # We eliminate the aromatic bonds by Kekulizing the atom in BLA selection.
        if bond_type > 1.1:
            double_idxs.append(bond_id)
        else:
            single_idxs.append(bond_id)

    single_bond_lengths = bond_lengths.sel(feature_indices=single_idxs)
    double_bond_lengths = bond_lengths.sel(feature_indices=double_idxs)

    BLA_res = single_bond_lengths.mean('descriptor') - double_bond_lengths.mean(
        'descriptor'
    )

    # Add back in the dimension we lost due to the `mean()`call
    BLA_res = BLA_res.expand_dims('descriptor')

    # Set metadatata
    BLA_res.name = "BLA"
    BLA_res.attrs.update(
        {
            "units": position_data.attrs["units"],
            "unitdim": "length",
            "long_name": "Bond length alternation (BLA)",
        }
    )

    descriptor_tex = [r"BLA"]
    descriptor_name = [r'BLA']
    descriptor_type: list[FeatureTypeLabel] = ['bla'] * len(descriptor_tex)
    return _assign_descriptor_coords(
        BLA_res,
        feature_descriptors=[(1, -1)],
        feature_type=descriptor_type,
        feature_tex_label=descriptor_tex,
        feature_name=descriptor_name,
    )
