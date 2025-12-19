from shnitsel._contracts import needs
import xarray as xr

from shnitsel.core.typedefs import AtXYZ
from shnitsel.filtering.structure_selection import BondDescriptor, StructureSelection
from shnitsel.geo.geocalc_.distances import get_distances
from shnitsel.geo.geocalc_.helpers import (
    _empty_descriptor_results,
    _get_default_selection,
)


@needs(dims={'atom', 'direction'})
def get_max_chromophor_BLA(
    atXYZ: AtXYZ,
    structure_selection: StructureSelection | None = None,
    SMARTS: str | None = None,
    num_double_bonds: int | None = None,
    allowed_chain_elements: str = "#6,#7,#8,#15,#16",
    max_considered_BLA_double_bonds: int = 50,
) -> xr.DataArray:
    """Calculate bond length alternation value (BLA) for the maximum chromophor
    in the provided `structure_selection` or the maximum chromophor in the structure
    represented by `atXYZ`.

    Parameters
    ----------
    atXYZ
        An :py:class:`xarray.DataArray` of molecular coordinates, with dimensions
        at least`atom` and `direction`
    structure_selection, optional
        Object encapsulating feature selection on the structure whose positional information is provided in `atXYZ`.
        If this argument is omitted altogether, a default selection for all bonds and atoms within the structure is created.
    SMARTS, optional
        SMARTS string to match for the maximum chromophor.
    num_double_bonds, optional
        The specified number of double bonds for the maximum chromophor.
    allowed_chain_elements, str
        SMARTS atomic specification, i.e. comma-separated list of element descriptors (default: C,N,O,P,S represented as '#6,#7,#8,#15,#16').
    max_considered_BLA_double_bonds, optional
        Maximum number of double bonds in a BLA chromophor if automatic maximum size detection is performed. Defaults to 50.

    Returns
    -------
        An :py:class:`xarray.DataArray` of the BLA for the maximum-length chromophor (alternating double bonds)

    Raises
    -------
    ValueError
        If the maximum chromophor within the provided selection or the entire molecule is not unique.

    """

    structure_selection = _get_default_selection(
        structure_selection, atXYZ=atXYZ, default_levels=['atoms', 'bonds']
    )

    BLA_selection = structure_selection.select_BLA_chromophor(
        BLA_smarts=SMARTS,
        num_double_bonds=num_double_bonds,
        allowed_chain_elements=allowed_chain_elements,
        max_considered_BLA_double_bonds=max_considered_BLA_double_bonds,
    )

    if len(BLA_selection.bonds_selected) == 0:
        return _empty_descriptor_results(atXYZ)

    bond_lengths = get_distances(atXYZ, BLA_selection)

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
    return BLA_res
    # data = data.expand_dims('descriptor')

    # # joined_atom_idxs = tuple(reduce(lambda x, y: set(x).union(y), atom_idxs))
    # joined_atom_idxs = tuple(set(sum(bond_idxs, ())))
    # format_str = 'BLA$_{' + ','.join(['%d'] * len(joined_atom_idxs)) + '}$'
    # return _assign_descriptor_coords(
    #     data,
    #     # *std_args,
    #     [joined_atom_idxs],
    #     [sum(bond_idxs, ())],
    #     [sum(bond_types, ())],
    #     [rc.Mol()],
    #     format_str,
    #     per_atom_coords=False,
    # )
