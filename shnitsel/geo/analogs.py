from dataclasses import dataclass
from functools import reduce
from itertools import combinations
import logging
from operator import and_
from typing import Any, Hashable, Iterable, Mapping, Sequence, TypeVar, overload

import numpy as np
import rdkit.Chem as rc
import xarray as xr

from shnitsel.data.dataset_containers import wrap_dataset
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.multi_indices import expand_midx
from shnitsel.bridges import construct_default_mol, set_atom_props
from shnitsel.data.tree.data_group import DataGroup
from shnitsel.data.tree.data_leaf import DataLeaf
from shnitsel.data.tree.node import TreeNode
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.filtering.structure_selection import SMARTSstring

DatasetOrArray = TypeVar(
    "DatasetOrArray", bound=ShnitselDataset | xr.Dataset | xr.DataArray
)


@dataclass
class StructureMapping:
    _submol: rc.Mol
    _orig_mol: rc.Mol
    _full_mapping: Mapping[int, int]

    def __init__(
        self, original_mol: rc.Mol, res_mol: rc.Mol, mapping: Mapping[int, int]
    ):
        self._submol = res_mol
        self._orig_mol = original_mol
        self._full_mapping = {k: v for k, v in mapping.items() if v >= 0}

    @overload
    def __call__(
        self, ds_or_da: TreeNode[Any, DatasetOrArray]
    ) -> TreeNode[Any, DatasetOrArray]: ...

    @overload
    def __call__(self, ds_or_da: DatasetOrArray) -> DatasetOrArray: ...

    def __call__(
        self, ds_or_da: DatasetOrArray | TreeNode[Any, DatasetOrArray]
    ) -> DatasetOrArray | TreeNode[Any, DatasetOrArray]:
        return self.apply(ds_or_da=ds_or_da)

    @overload
    def apply(
        self, ds_or_da: TreeNode[Any, DatasetOrArray]
    ) -> TreeNode[Any, DatasetOrArray]: ...

    @overload
    def apply(self, ds_or_da: DatasetOrArray) -> DatasetOrArray: ...

    def apply(
        self, ds_or_da: DatasetOrArray | TreeNode[Any, DatasetOrArray]
    ) -> DatasetOrArray | TreeNode[Any, DatasetOrArray]:
        if isinstance(ds_or_da, TreeNode):

            def map_ds(ds: DatasetOrArray) -> DatasetOrArray:
                return self.apply(ds)

            return ds_or_da.map_data(map_ds)
        else:
            orig_ids, new_ids = zip(*list(self._full_mapping.items()))
            mol = rc.Mol(self._submol)
            return (
                ds_or_da.isel(atom=list(orig_ids))
                .assign_coords(atom=("atom", list(new_ids)))
                .sortby("atom")
                .assign_attrs(__mol=mol)
                .assign_coords(
                    __mol=xr.DataArray(mol),
                    charge=xr.DataArray(rc.GetFormalCharge(mol)),
                )
            )


def _find_atom_pairs(mol: rc.Mol, atoms: Sequence[int]) -> list[int]:
    """Method to find all atom pairs that constitute a bond in the molecule
    and return the associated bond ids as a list.

    The list of bond ids can serve as a path within a molecule to support the
    extraction of a submolecule that contains also bonds at the edge of a
    SMARTS string structure.

    Parameters
    ----------
    mol : rc.Mol
        The molecule to find the bonds within
    atoms : Sequence[int]
        The set of atoms among which we want to try to find the bonds
        to construct a path of bonds.

    Returns
    -------
    list[int]
        The list of the bond ids of all bonds within `mol` between the `atoms`.
    """
    # TODO: Might we need this elsewhere?
    res = []
    # For each pair of atoms, check whether we have the respective bond
    for i, j in combinations(atoms, 2):
        # Check if this is a bond
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is not None:
            # Add the bond id if it was a bond
            res.append(bond.GetIdx())
    return res


def _substruct_match_to_submol(mol: rc.Mol, substruct_match: tuple[int, ...]) -> rc.Mol:
    """Build a sub-mol from a substructure match.

    This is used for analogs mapping.

    Parameters
    ----------
    mol : rc.Mol
        The mol on which the substructure has matched
    substruct_match : tuple[int,...]
        The indices of the atoms that have matched.

    Returns
    -------
    rc.Mol
        The mol object of the relevant substructure match.
    """
    # Store consistent atom order before extracting
    pattern_match_list = [-1] * mol.GetNumAtoms()
    for idx_in_pattern, idx_in_mol in enumerate(substruct_match):
        pattern_match_list[idx_in_mol] = idx_in_pattern
    mol_with_pattern_ids = set_atom_props(mol, pattern_idx=pattern_match_list)

    # Extract submol
    bond_path = _find_atom_pairs(mol_with_pattern_ids, substruct_match)
    res_mol = rc.PathToSubmol(mol_with_pattern_ids, bond_path)

    # Renumber atoms using stored order
    res_map = [-1] * res_mol.GetNumAtoms()
    for a in res_mol.GetAtoms():
        res_map[a.GetIntProp("pattern_idx")] = a.GetIdx()
    res_mol = rc.RenumberAtoms(res_mol, res_map)
    return res_mol


def _substruct_match_to_mapping(
    mol: rc.Mol, substruct_match: tuple[int, ...]
) -> tuple[rc.Mol, StructureMapping]:
    """Build a mapping of the mol ids to the ids of the substructure.

    This allows for selection of the atoms within a dataset based on a substructure match.

    Parameters
    ----------
    mol : rc.Mol
        The original mol on which the substructure match was obtained
    substruct_match : tuple[int,...]
        The index list of the substructure match.

    Returns
    -------
    rc.Mol
        The matched submol object.
    StructureMapping
        The mapping of original mol atom indices to submol indices or -1 if no longer present.
        Atoms no longer present in the final submol may either not be in the keys of the mapping or have a negative value
        associated with their index.
        Has an `.apply()` function to apply the mapping to datasets or data arrays
    """
    # Store consistent atom order before extracting
    res_map: dict[int, int] = dict()
    pattern_match_list = [-1] * mol.GetNumAtoms()
    for idx_in_pattern, idx_in_mol in enumerate(substruct_match):
        pattern_match_list[idx_in_mol] = idx_in_pattern
        res_map[idx_in_mol] = idx_in_pattern

    mol_with_pattern_ids = set_atom_props(mol, pattern_idx=pattern_match_list)

    # Extract submol
    bond_path = _find_atom_pairs(mol_with_pattern_ids, substruct_match)
    res_mol = rc.PathToSubmol(mol_with_pattern_ids, bond_path)

    # Renumber atoms using stored order
    res_list = [-1] * res_mol.GetNumAtoms()
    for a in res_mol.GetAtoms():
        res_list[a.GetIntProp("pattern_idx")] = a.GetIdx()
    res_mol = rc.RenumberAtoms(res_mol, res_list)

    # Clear the pattern idx markers
    set_atom_props(res_mol, inplace=True, pattern_idx=False)
    return res_mol, StructureMapping(mol_with_pattern_ids, res_mol, res_map)


def get_MCS_smarts(mols: Iterable[rc.Mol]) -> SMARTSstring:
    """Helper function to get the maximum common substructure (MCS) SMARTS string
    for a sequence of Molecular structures for further processing.

    Parameters
    ----------
    mols : Iterable[rc.Mol]
        The molecular structures to get the maximum common substructre between.

    Returns
    -------
    SMARTSstring
        The MCS SMARTS string.
    """
    from rdkit.Chem import rdFMCS

    mcs_settings = rdFMCS.MCSParameters()
    # TODO: FIXME: Do we want this `any heavy atom` match?
    # mcs_settings.AtomCompareParameters = rdFMCS.AtomCompare.CompareAnyHeavyAtom
    # # TODO: FIXME: Do we want bonds to always compare EXACT bond order?
    # mcs_settings.BondCompareParameters = rdFMCS.BondCompare.CompareOrder

    # TODO: We should probably consider the default config for FindMCS, which includes chains matching rings among others.
    substructure_smarts = rdFMCS.FindMCS(
        list(mols),
        parameters=mcs_settings,
    ).smartsString

    return substructure_smarts


def identify_analogs_mappings(
    mols: Mapping[Hashable | int, rc.Mol],
    smarts: SMARTSstring = "",
) -> tuple[SMARTSstring, Mapping[Hashable | int, StructureMapping]]:
    """Helper function to generate a maximum common substructure match and
    from that extract substructure mappings for each of the provided molecules.

    If provided a `smarts` string, the MCS will be skipped and instead an attempt will be made to match
    `smarts` against all provided structures.

    Parameters
    ----------
    mols : Mapping[Hashable  |  int, rc.Mol]
        The molecular structures to use as a basis for the MCS
        analysis or for finding the `smarts` string in if provided.

    Returns
    -------
    tuple[SMARTSstring, Mapping[Hashable | int, StructureMapping]]
        First the resulting (or used) SMARTS string for the structure.
        Then the Mapping between original keys and the resulting `StructureMapping` object
        that can be applied to the original data.
    """
    if not smarts:
        # TODO: We should probably consider the default config for FindMCS, which includes chains matching rings among others.
        substructure_smarts = get_MCS_smarts(mols.values())
    else:
        substructure_smarts = smarts

    # Get the mol match from the found MCS
    search_substructure = rc.MolFromSmarts(substructure_smarts)

    # Generate mappings
    results: dict[Hashable | int, StructureMapping] = {}
    for key, mol in mols.items():
        substruct_matches: tuple[int, ...] = mol.GetSubstructMatch(search_substructure)
        res_mol, res_mapping = _substruct_match_to_mapping(mol, substruct_matches)
        set_atom_props(res_mol, inplace=True, atomNote=True)
        results[key] = res_mapping

    logging.info(f"MCS SMARTS is: {substructure_smarts}")
    return substructure_smarts, results


# TODO: FIXME: This should accept a smarts of elements that should be considered equivalent for substructs.
def _list_analogs(
    ensembles: Mapping[Hashable | int, xr.DataArray],
    smarts: SMARTSstring = "",
    vis: bool = False,
) -> Mapping[Hashable | int, xr.DataArray]:
    """Extract a common moiety from a selection of ensembles.

    By default, this attempts to find the largest possible match using equivalence of any heavy atoms.
    H-atoms can only match other H-atoms.

    Parameters
    ----------
    ensembles : Mapping[Hashable | int, xr.DataArray]
        An ``Iterable`` of ``xr.DataArray``s, each containing the geometries of an ensemble of
        trajectories for a different compound; they
    smarts : SMARTSstring, optional
        A SMARTS-string indicating the moiety to cut out of each compound;
        in each case, the match returned by :py:func:`rdkit.Chem.Mol.GetSubstrucMatch`
        (not necessarily the only possible match) will be used;
        if no SMARTS is provided, a minimal common submol will be extracted using
        ``rdFMCS.FindMCS``
    vis : bool, default=False
        Whether to display a visual indication of the matches.

    Returns
    -------
        An ``Iterable`` of ``xr.DataArray``s
    """
    from rdkit.Chem import rdFMCS

    # The visualization should be separated
    if vis:
        from IPython.display import display
    else:
        display = None

    mcs_settings = rdFMCS.MCSParameters()
    # TODO: FIXME: Do we want this `any heavy atom` match?
    mcs_settings.AtomCompareParameters = rdFMCS.AtomCompare.CompareAnyHeavyAtom
    # TODO: FIXME: Do we want bonds to always compare EXACT bond order?
    mcs_settings.BondCompareParameters = rdFMCS.BondCompare.CompareOrder

    # Get mol representation per ensemble
    # This incorporates the atNames property, so
    # equivalence would need to be applied beforehand
    mols = {k: construct_default_mol(x) for k, x in ensembles.items()}
    if not smarts:
        # TODO: We should probably consider the default config for FindMCS, which includes chains matching rings among others.
        smarts = rdFMCS.FindMCS(
            mols.values(),
            parameters=mcs_settings,
        ).smartsString

    # Get the mol match either from the provided SMARTS or the found MCS
    search = rc.MolFromSmarts(smarts)

    results = {}
    mol_grid = []
    for key in ensembles.keys():
        mol = mols[key]
        compound_geo = ensembles[key]
        idxs = list(mol.GetSubstructMatch(search))
        res_mol = _substruct_match_to_submol(mol, idxs)
        set_atom_props(res_mol, inplace=True, atomNote=True)

        if vis:
            atom_labels = [""] * mol.GetNumAtoms()
            for patt_idx, mol_idx in enumerate(idxs):
                atom_labels[mol_idx] = f"{mol_idx}:{patt_idx}"
            vis_orig = rc.Mol(mol)  # avoid mutating original
            set_atom_props(vis_orig, inplace=True, atomNote=atom_labels)

            atom_labels = [
                f"{mol_idx}:{patt_idx}" for patt_idx, mol_idx in enumerate(idxs)
            ]
            vis_patt = rc.Mol(search)  # avoid mutating original
            set_atom_props(vis_patt, inplace=True, atomNote=atom_labels)

            mol_grid.append([vis_orig, vis_patt, res_mol])

        range_ = range(len(idxs))
        results[key] = (
            compound_geo.isel(atom=idxs)
            .assign_coords(atom=range_)
            .sortby("atom")
            .assign_attrs(__mol=res_mol)
        )

    if vis and display is not None:
        display(rc.Draw.MolsMatrixToGridImage(mol_grid))

    return results


# def _combine_compounds_unstacked(compounds, names=None, concat_kws=None):
#     if concat_kws is None:
#         concat_kws = {}

#     coord_names = [set(x.coords) for x in compounds]
#     coords_shared = reduce(and_, coord_names)
#     compounds = [
#         x.drop_vars(set(x.coords).difference(coords_shared)) for x in compounds
#     ]
#     if names is None:
#         names = range(len(compounds))
#     compounds = [
#         x.assign_coords(
#             {
#                 'compound': ('trajid', np.full(x.sizes['trajid'], name)),
#                 'traj': x.trajid,
#             }
#         )
#         .reset_index('trajid')
#         .set_xindex(['compound', 'traj'])
#         for x, name in zip(compounds, names)
#     ]

#     return xr.concat(compounds, dim='trajid', **concat_kws)


# def _combine_compounds_stacked(compounds, names=None, concat_kws=None):
#     if concat_kws is None:
#         concat_kws = {}

#     concat_dim = 'frame'

#     coord_names = [set(x.coords) for x in compounds]
#     coords_shared = reduce(and_, coord_names)
#     compounds = [
#         x.drop_vars(set(x.coords).difference(coords_shared)) for x in compounds
#     ]

#     if names is None:
#         names = range(len(compounds))

#     # Which coords are unique on a compound level? So far:
#     c_per_compound = ['atNames', 'atNums']

#     per_compound = {
#         crd: xr.concat(
#             [obj[crd] for obj in compounds],
#             dim='compound_',
#         )
#         for crd in c_per_compound
#         if all(crd in obj.coords for obj in compounds)
#     }

#     compounds = [
#         expand_midx(x, 'frame', 'compound', name)
#         .drop_dims('trajid_')
#         .drop_vars(c_per_compound, errors='ignore')
#         for x, name in zip(compounds, names)
#     ]

#     res = xr.concat(
#         compounds, dim=concat_dim, **({'combine_attrs': 'drop_conflicts'} | concat_kws)
#     )
#     res = res.assign_coords(per_compound)
#     res = res.assign_coords(compound_=names)

#     if any('time_' in x.dims for x in compounds):
#         time_only = xr.concat(
#             [obj.drop_dims(['frame', 'trajid_'], errors='ignore') for obj in compounds],
#             dim='time_',
#             **concat_kws,
#         )
#         res = res.assign(time_only)

#     # TODO: consider using MajMinIndex
#     return res


@overload
def extract_analogs(
    ensembles: TreeNode[Any, DatasetOrArray],
    smarts: SMARTSstring = "",
    vis: bool = False,
    *,
    concat_kws: dict[str, Any] | None = None,
) -> TreeNode[Any, DatasetOrArray] | None: ...


@overload
def extract_analogs(
    ensembles: Mapping[Hashable | int, DatasetOrArray],
    smarts: SMARTSstring = "",
    vis: bool = False,
    *,
    concat_kws: dict[str, Any] | None = None,
) -> Mapping[Hashable | int, DatasetOrArray] | None: ...


@overload
def extract_analogs(
    ensembles: Sequence[DatasetOrArray],
    smarts: SMARTSstring = "",
    vis: bool = False,
    *,
    concat_kws: dict[str, Any] | None = None,
) -> Sequence[DatasetOrArray] | None: ...


# TODO: FIXME: We should add a method that simply punches out a substructure from a match.
def extract_analogs(
    ensembles: (
        TreeNode[Any, DatasetOrArray]
        | Mapping[Hashable | int, DatasetOrArray]
        | Sequence[DatasetOrArray]
    ),
    # TODO: FIXME: Do we need a second smarts to restrict the search?
    smarts: SMARTSstring = "",
    vis: bool = False,
    *,
    concat_kws: dict[str, Any] | None = None,
) -> (
    TreeNode[Any, DatasetOrArray]
    | Mapping[Hashable | int, DatasetOrArray]
    | Sequence[DatasetOrArray]
    | None
):
    """Combine ensembles for different compounds by finding the
    moieties they have in common

    Parameters
    ----------
    ensembles : TreeNode[Any, DatasetOrArray] | Mapping[Hashable | int, DatasetOrArray] | Sequence[DatasetOrArray]
        Input of Datasets or DataArrays or Shnitsel Wrappers optionally in a tree structure,
        each containing the geometries of an ensemble of trajectories for a different compound or structure.
            - If the ensemble is provided as a tree, the result will be a tree of a mostly identical structure.
                A grouping operation may be performed beforehand to avoid different structures to be in the same group.
            - If the input is a mapping, the keys will be preserved and the mappings will be applied to each entry
            - If a sequence is provided, the order of inputs will be preserved and the mapping will be applied to each entry in order.
    smarts : SMARTSstring
        A SMARTS-string indicating the moiety to cut out of each compound;
        in each case, the match returned by :py:func:`rdkit.Chem.Mol.GetSubstructMatch`
        (not necessarily the only possible match) will be used;
        if no SMARTS is provided, a minimal common submol will be extracted using
        ``rdFMCS.FindMCS``
    vis : bool, default=False
        Deprecated; Whether to display a visual indication of the match, by default False
    **concat_kws
        Deprecated; Keyword arguments for internal calls to ``xr.concat``

    Returns
    -------
    TreeNode[Any, DatasetOrArray]:
        An tree holding the analog substructures in its leaves.
    Mapping[Hashable | int, DatasetOrArray]
    | Sequence[DatasetOrArray]
        Either a mapping or a sequence of `xr.Dataset` or `xr.DataArray` of trajectories,
        holding the mapped inputs from `ensembles`.

    Raises
    ------
    ValueError
        If the ensembles provided could not be brought into agreement.
    AssertionError
        If the tree is not of a suppported format.
    """
    if isinstance(ensembles, TreeNode):
        assert isinstance(ensembles, ShnitselDB), (
            "Analogs matching currently only supported on full shnitsel DB trees with a `CompoundGroup` level underneath them"
        )

        grouped_tree = ensembles.group_data_by_metadata()
        assert grouped_tree is not None, (
            "Tree could not be grouped by metadata. You may have provided a tree with wrong data inside."
        )

        # Get mol and path from the subtree
        def path_and_mol_from_group(
            group: TreeNode[Any, DatasetOrArray],
        ) -> TreeNode[Any, tuple[str, rc.Mol]] | None:
            if not isinstance(group, DataGroup):
                # Something went wrong, let's keep
                logging.warning(
                    "Mapping over groups yielded non-group tree node: %s", group.path
                )
                return None

            for child in group.subleaves.values():
                if child.has_data:
                    mol = construct_default_mol(child.data)
                    if mol:
                        return group.construct_copy(
                            children={
                                "_agg": DataLeaf(
                                    name="_agg" + str(child.name),
                                    data=(group.path, mol),
                                )
                            }
                        )
            return None

        # Collect all mols for groups from the tree with their path for later lookup
        path_mol_pairs = grouped_tree.map_filtered_nodes(
            filter_func=lambda x: isinstance(x, DataGroup) and len(x.subleaves) > 0,
            map_func=path_and_mol_from_group,
        ).collect_data()

        path_mol_map: dict[Hashable | int, rc.Mol] = {
            path: mol for path, mol in path_mol_pairs
        }

        res_SMARTS, res_mappings = identify_analogs_mappings(
            path_mol_map, smarts=smarts
        )

        # print(res_SMARTS)
        if not smarts:
            logging.info("Substructure matching resulted in SMARTS: %s", res_SMARTS)

        def patch_groups(
            group: TreeNode[Any, DatasetOrArray],
        ) -> TreeNode[Any, DatasetOrArray] | None:
            # print(f"Considering: {group.path} in {res_mappings}")
            if group.path in res_mappings:
                # print(f"Mapping group: {group.path}")
                struct_map = res_mappings[group.path]
                # print(struct_map._submol)

                def map_traj(x: DatasetOrArray) -> DatasetOrArray:
                    return wrap_dataset(struct_map.apply(x))

                return group.map_data(map_traj)
            return group.map_data(wrap_dataset)

        # Map the tree with the appropriate structure mapping:
        return grouped_tree.map_filtered_nodes(
            filter_func=lambda x: isinstance(x, DataGroup) and len(x.subleaves) > 0,
            map_func=patch_groups,
        )

    input_mol_map: dict[Hashable | int, rc.Mol]
    if isinstance(ensembles, Sequence):
        input_mol_map = {
            i: construct_default_mol(ensemble) for i, ensemble in enumerate(ensembles)
        }
    else:
        input_mol_map = {
            k: construct_default_mol(ensemble) for k, ensemble in ensembles.items()
        }
    res_smarts, res_mappings = identify_analogs_mappings(input_mol_map, smarts=smarts)

    # analogs = list_analogs(ensembles, smarts=smarts, vis=vis)
    # if all(is_stacked(x) for x in analogs):
    #     res = _combine_compounds_stacked(analogs, names=names, concat_kws=concat_kws)
    # elif not any(is_stacked(x) for x in analogs):
    #     res = _combine_compounds_unstacked(analogs, names=names, concat_kws=concat_kws)
    # else:
    #     raise ValueError("Inconsistent formats")

    # mols = [x.attrs['mol'] for x in analogs]
    # mol = mols[0]  # TODO: Try replacing with search pattern object
    # res = res.assign_attrs(mol=mol).assign_coords(mols=('compound_', mols))
    # if 'atXYZ' in res:
    #     res['atXYZ'].attrs['__mol'] = mol
    # return res

    # Map results
    if isinstance(ensembles, Sequence):
        return [
            wrap_dataset(res_mappings[i].apply(ensemble))
            for i, ensemble in enumerate(ensembles)
        ]
    else:
        return {
            k: wrap_dataset(res_mappings[k].apply(ensemble))
            for k, ensemble in ensembles.items()
        }
