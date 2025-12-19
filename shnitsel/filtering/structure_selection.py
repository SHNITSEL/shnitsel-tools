from dataclasses import dataclass
from itertools import combinations
import logging
from typing import Iterator, Literal, Self, Sequence, TypeAlias
import rdkit
import xarray as xr

from rdkit.Chem.rdchem import Mol

from shnitsel.bridges import to_mol
from shnitsel.vis.colormaps import hex2rgb, st_yellow
from IPython.display import SVG


AtomDescriptor: TypeAlias = int
BondDescriptor: TypeAlias = tuple[int, int]
AngleDescriptor: TypeAlias = tuple[int, int, int]
DihedralDescriptor: TypeAlias = tuple[int, int, int, int]
PyramidsDescriptor: TypeAlias = tuple[
    int, tuple[int, int, int]
]  # Center of the pyramid is the separate int

FeatureDescriptor: TypeAlias = (
    AtomDescriptor
    | BondDescriptor
    | AngleDescriptor
    | DihedralDescriptor
    | PyramidsDescriptor
)

FeatureList: TypeAlias = list[FeatureDescriptor]

ActiveFlag: TypeAlias = bool

FeatureLevelType: TypeAlias = Literal[
    'atoms', 'bonds', 'angles', 'dihedrals', 'pyramids'
]
FeatureLevelOptions: TypeAlias = FeatureLevelType | Literal[1, 2, 3, 4, 5]

FEATURE_LEVELS: list[FeatureLevelType] = [
    'atoms',
    'bonds',
    'angles',
    'dihedrals',
    'pyramids',
]

FeatureTypeLabel: TypeAlias = Literal['pyr', 'pos', 'dist', 'angle', 'dih', 'cos', 'sin']


@dataclass
class StructureSelection:
    """Class to keep track of a (sub-)selection of structural features."""

    mol: Mol | None

    atoms: set[AtomDescriptor]
    atoms_types: dict[AtomDescriptor, str]
    atoms_selected: set[AtomDescriptor]

    bonds: set[BondDescriptor]
    bonds_types: dict[BondDescriptor, float]
    bonds_selected: set[BondDescriptor]

    angles: set[AngleDescriptor]
    angles_types: dict[AngleDescriptor, bool]
    angles_selected: set[AngleDescriptor]

    dihedrals: set[DihedralDescriptor]
    dihedrals_types: dict[DihedralDescriptor, bool]
    dihedrals_selected: set[DihedralDescriptor]

    pyramids: set[PyramidsDescriptor]
    pyramids_types: dict[PyramidsDescriptor, bool]
    pyramids_selected: set[PyramidsDescriptor]

    def copy_or_update(
        self,
        mol: Mol | None = None,
        atoms: set[AtomDescriptor] | None = None,
        atoms_selected: set[AtomDescriptor] | None = None,
        atoms_types: dict[AtomDescriptor, str] | None = None,
        bonds: set[BondDescriptor] | None = None,
        bonds_selected: set[BondDescriptor] | None = None,
        bonds_types: dict[BondDescriptor, float] | None = None,
        angles: set[AngleDescriptor] | None = None,
        angles_selected: set[AngleDescriptor] | None = None,
        angles_types: dict[AngleDescriptor, bool] | None = None,
        dihedrals: set[DihedralDescriptor] | None = None,
        dihedrals_selected: set[DihedralDescriptor] | None = None,
        dihedrals_types: dict[DihedralDescriptor, bool] | None = None,
        pyramids: set[PyramidsDescriptor] | None = None,
        pyramids_selected: set[PyramidsDescriptor] | None = None,
        pyramids_types: dict[PyramidsDescriptor, bool] | None = None,
        inplace: bool = False,
    ) -> Self:
        """Function to create a copy with replaced member values.

        Meant as a helper for the `Frozen` logic of the selection, i.e. method calls return a new instance
        instead of updating the existing instance.


        Args:
            mol (Mol | None, optional): The new RDKit Mol object to assign to this selection. Should usually not be updated, but can be for completeness. Defaults to None.
            atoms (set[AtomDescriptor] | None, optional): New indices of atoms. Defaults to None.
            atoms_selected (set[AtomDescriptor] | None, optional): New indices of selected atoms. Defaults to None.
            atoms_types (dict[AtomDescriptor, str] | None, optional): New metadata dict for atom indices. Defaults to None.
            bonds (set[BondDescriptor] | None, optional): Bond indices set. Defaults to None.
            bonds_selected (set[BondDescriptor] | None, optional): Set of bond indices that have been selected. Defaults to None.
            bonds_types (dict[BondDescriptor, float] | None, optional): Dict with metadata for the bonds. Defaults to None.
            angles (set[AngleDescriptor] | None, optional): Set of all indices of angles. Defaults to None.
            angles_selected (set[AngleDescriptor] | None, optional): Set of selected indices of angles. Defaults to None.
            angles_types (dict[AngleDescriptor, bool] | None, optional): Dict with metadata about the angles. Defaults to None.
            dihedrals (set[DihedralDescriptor] | None, optional): Set of all indices of dihedrals in the structure. Defaults to None.
            dihedrals_selected (set[DihedralDescriptor] | None, optional): Set of selected indices of dihedrals. Defaults to None.
            dihedrals_types (dict[DihedralDescriptor, bool] | None, optional): Dict with metadata about dihedrals. Defaults to None.
            pyramids (set[PyramidsDescriptor] | None, optional): Set of all indices of pyramids in the structure. Defaults to None.
            pyramids_selected (set[PyramidsDescriptor] | None, optional): Set of selected indices of pyramids. Defaults to None.
            pyramids_types (dict[PyramidsDescriptor, bool] | None, optional): Dict with metadata about pyramids. Defaults to None.
            inplace (bool, optional): Flag to allow for in-place updates instead of returning a new cop. Defaults to False.

        Returns:
            StructureSelection: The selection update with the new members set. Can either be a copy if `inplace=False` or the old instance with updated members otherwise.
        """
        if inplace:
            # Update and create
            if mol is not None:
                self.mol = mol

            if atoms is not None:
                self.atoms = atoms
            if atoms_selected is not None:
                self.atoms_selected = atoms_selected
            if atoms_types is not None:
                self.atoms_types = atoms_types

            if bonds is not None:
                self.bonds = bonds
            if bonds_selected is not None:
                self.bonds_selected = bonds_selected
            if bonds_types is not None:
                self.bonds_types = bonds_types

            if angles is not None:
                self.angles = angles
            if angles_selected is not None:
                self.angles_selected = angles_selected
            if angles_types is not None:
                self.angles_types = angles_types

            if dihedrals is not None:
                self.dihedrals = dihedrals
            if dihedrals_selected is not None:
                self.dihedrals_selected = dihedrals_selected
            if dihedrals_types is not None:
                self.dihedrals_types = dihedrals_types

            if pyramids is not None:
                self.pyramids = pyramids
            if pyramids_selected is not None:
                self.pyramids_selected = pyramids_selected
            if pyramids_types is not None:
                self.pyramids_types = pyramids_types

            return self
        else:
            if mol is None:
                mol = self.mol

            if atoms is None:
                atoms = self.atoms
            if atoms_selected is None:
                atoms_selected = self.atoms_selected
            if atoms_types is None:
                atoms_types = self.atoms_types

            if bonds is None:
                bonds = self.bonds
            if bonds_selected is None:
                bonds_selected = self.bonds_selected
            if bonds_types is None:
                bonds_types = self.bonds_types

            if angles is None:
                angles = self.angles
            if angles_selected is None:
                angles_selected = self.angles_selected
            if angles_types is None:
                angles_types = self.angles_types

            if dihedrals is None:
                dihedrals = self.dihedrals
            if dihedrals_selected is None:
                dihedrals_selected = self.dihedrals_selected
            if dihedrals_types is None:
                dihedrals_types = self.dihedrals_types

            if pyramids is None:
                pyramids = self.pyramids
            if pyramids_selected is None:
                pyramids_selected = self.pyramids_selected
            if pyramids_types is None:
                pyramids_types = self.pyramids_types

            return type(self)(
                mol=mol,
                atoms=atoms,
                atoms_selected=atoms_selected,
                atoms_types=atoms_types,
                bonds=bonds,
                bonds_selected=bonds_selected,
                bonds_types=bonds_types,
                angles=angles,
                angles_selected=angles_selected,
                angles_types=angles_types,
                dihedrals=dihedrals,
                dihedrals_selected=dihedrals_selected,
                dihedrals_types=dihedrals_types,
                pyramids=pyramids,
                pyramids_selected=pyramids_selected,
                pyramids_types=pyramids_types,
            )

    @classmethod
    def init_from_dataset(
        cls: type[Self],
        dataset: xr.Dataset,
        default_selection: Sequence[FeatureLevelOptions] = [
            'atoms',
            'bonds',
        ],
        to2D: bool = True,
    ) -> Self:
        """Alternative constructor that creates an initial StructureSelection object from a dataset using the entire structural information in it.

        Args:
            cls (type[StructureSelection]): The type of this StructureSelection so that we can create instances of it.
            dataset (xr.Dataset): The dataset to extract the structure information out of.
                Must have at least `atXYZ` variable and a `atom` dimension.
                Ideally, an `atom` coordinate for feature selection is also provided.
                Should only represent a single frame of data.
            default_selection (Sequence[FeatureLevelOptions], optional): List of features to activate as selected by default. Defaults to [ 'atoms', 'bonds', ].
            to2D (bool, optional): Flag to control whether a mol representation is converted to a 2d projection before use for visualization.

        Raises:
            ValueError: If no structural information could be extracted from the dataset

        Returns:
            StructureSelection: A structure selection object initially covering all atoms and structural features.
        """
        filtered_dataset = dataset.squeeze()

        if 'frame' in filtered_dataset.dims or 'time' in filtered_dataset.dims:
            raise ValueError(
                "The dataset should not contain frame or data but represent a single frame of data. \n"
                f"Had dimensions : {filtered_dataset.dims}"
            )

        # TODO: FIXME: Consider the charges needing to be set from the dataset settings.s
        mol = to_mol(
            filtered_dataset.atXYZ,
            to2D=to2D,
        )
        # Create an initial state selection
        return cls.init_from_mol(mol, default_selection=default_selection)

    @classmethod
    def init_from_mol(
        cls: type[Self],
        mol: Mol,
        default_selection: Sequence[FeatureLevelOptions] = [
            'atoms',
            'bonds',
        ],
    ) -> Self:
        """Alternative constructor that creates an initial StructureSelection object from an RDKit Mol object

        Args:
            cls (type[StructureSelection]): The type of this StructureSelection so that we can create instances of it.
            mol (rdkit.rdchem.Mol): The RDKit Mol object to extract all initial structural information out of
            default_selection (Sequence[FeatureLevelOptions], optional): List of features to activate as selected by default. Defaults to [ 'atoms', 'bonds', ].

        Raises:
            ValueError: If no structural information could be extracted from the dataset.

        Returns:
            StructureSelection: A structure selection object initially covering all atoms and structural features.
        """
        # TODO: FIXME: Implement actual feature selection with geomatch

        atoms = set()
        atoms_selected = set()
        atoms_types = dict()
        bonds = set()
        bonds_selected = set()
        bonds_types = dict()
        angles = set()
        angles_selected = set()
        angles_types = dict()
        dihedrals = set()
        dihedrals_selected = set()
        dihedrals_types = dict()
        pyramids: set[PyramidsDescriptor] = set()
        pyramids_selected: set[PyramidsDescriptor] = set()
        pyramids_types = dict()

        default_selection = [
            StructureSelection._to_feature_level_str(x) for x in default_selection
        ]

        are_atoms_selected = 'atoms' in default_selection
        are_bonds_selected = 'bonds' in default_selection
        are_angles_selected = 'angles' in default_selection
        are_dihedrals_selected = 'dihedrals' in default_selection
        are_pyramids_selected = 'pyramids' in default_selection

        for atom in mol.GetAtoms():
            atomid = atom.GetIdx()
            atom_type = atom.GetSymbol()
            atoms.add(atomid)
            atoms_types[atomid] = atom_type

        if are_atoms_selected:
            atoms_selected.update(atoms)

        for bond in mol.GetBonds():
            beginIdx = bond.GetBeginAtomIdx()
            endIdx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondTypeAsDouble()
            bondId = (beginIdx, endIdx)
            bonds.add(bondId)
            bonds_types[bondId] = bond_type

        if are_bonds_selected:
            bonds_selected.update(bonds)

        for bond_j in mol.GetBonds():
            j = bond_j.GetBeginAtomIdx()
            k = bond_j.GetEndAtomIdx()

            # find atoms bonded to j (except k)
            neighbors_j = [
                nbr.GetIdx()
                for nbr in mol.GetAtomWithIdx(j).GetNeighbors()
                if nbr.GetIdx() != k
            ]

            # angles are (i, j, k)
            for i in neighbors_j:
                angles.add((i, j, k))

            # also angles (k, j, i) by symmetry
            neighbors_k = [
                nbr.GetIdx()
                for nbr in mol.GetAtomWithIdx(k).GetNeighbors()
                if nbr.GetIdx() != j
            ]
            for i in neighbors_k:
                angles.add((i, k, j))

        if are_angles_selected:
            angles_selected.update(angles)

        for bond_jk in mol.GetBonds():
            j = bond_jk.GetBeginAtomIdx()
            k = bond_jk.GetEndAtomIdx()

            # atoms bonded to j (except k)
            neighbors_j = [
                nbr.GetIdx()
                for nbr in mol.GetAtomWithIdx(j).GetNeighbors()
                if nbr.GetIdx() != k
            ]

            # atoms bonded to k (except j)
            neighbors_k = [
                nbr.GetIdx()
                for nbr in mol.GetAtomWithIdx(k).GetNeighbors()
                if nbr.GetIdx() != j
            ]

            # form dihedrals (i, j, k, l)
            for i in neighbors_j:
                for l in neighbors_k:
                    # TODO: FIXME: check if we want to exclude potential i=l
                    dihedrals.add((i, j, k, l))

                    # also handle reversed central bond direction (k–j)
                    # giving quadruples (i, k, j, l)
                    dihedrals.add((l, k, j, i))

        if are_dihedrals_selected:
            dihedrals_selected.update(dihedrals)

        for atom in mol.GetAtoms():
            i = atom.GetIdx()

            # get all neighbors of atom i
            neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]

            # a pyramid requires exactly (or at least) 3 bonded neighbors
            if len(neighbors) < 3:
                continue

            # choose all unique triplets of neighbors
            for j, k, l in combinations(neighbors, 3):
                pyramids.add((i, (j, k, l)))

        if are_pyramids_selected:
            pyramids_selected.update(pyramids)

        # Create an initial state selection
        return cls(
            mol=mol,
            atoms=atoms,
            atoms_selected=atoms_selected,
            atoms_types=atoms_types,
            bonds=bonds,
            bonds_selected=bonds_selected,
            bonds_types=bonds_types,
            angles=angles,
            angles_selected=angles_selected,
            angles_types=angles_types,
            dihedrals=dihedrals,
            dihedrals_selected=dihedrals_selected,
            dihedrals_types=dihedrals_types,
            pyramids=pyramids,
            pyramids_selected=pyramids_selected,
            pyramids_types=pyramids_types,
        )

    def select_all(
        self,
        feature_level: FeatureLevelOptions
        | Sequence[FeatureLevelOptions] = FEATURE_LEVELS,
        inplace: bool = False,
    ) -> Self:
        """Helper function to set all entries of a certain feature level to selected.

        By default marks all features as selected.

        Args:
            feature_level (FeatureLevelOptions | Sequence[FeatureLevelOptions], optional): The set of feature levels to mark as within the selection. Defaults to all FEATURE_LEVELS.
            inplace (bool, optional): Whether to update the selection in-place. Defaults to False.

        Returns:
            Self: The updated selection
        """
        if isinstance(feature_level, str) or not isinstance(feature_level, Sequence):
            feature_level = [feature_level]

        feature_levels = [self._to_feature_level_str(x) for x in feature_level]
        return self.copy_or_update(
            atoms_selected=self.atoms if 'atoms' in feature_levels else None,
            bonds_selected=self.bonds if 'bonds' in feature_levels else None,
            angles_selected=self.angles if 'angles' in feature_levels else None,
            dihedrals_selected=self.dihedrals
            if 'dihedrals' in feature_levels
            else None,
            pyramids_selected=self.pyramids if 'pyramids' in feature_levels else None,
            inplace=inplace,
        )

    def select_atoms(
        self,
        smarts_or_selection: str
        | Sequence[str]
        | AtomDescriptor
        | Sequence[AtomDescriptor]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        selection_list: Sequence
        if isinstance(smarts_or_selection, str) or isinstance(
            smarts_or_selection, AtomDescriptor
        ):
            selection_list = [smarts_or_selection]
        elif smarts_or_selection is None:
            selection_list = list(self.atoms.copy())
        else:
            selection_list = smarts_or_selection

        idx_set = []

        for entry in selection_list:
            if isinstance(entry, str):
                # Handle SMARTS parameter.
                if self.mol is None:
                    raise ValueError(
                        "No Molecule set in selection. Cannot match SMARTS."
                    )
                else:
                    matches = self.__match_pattern(self.mol, entry)

                    idx_set.extend(self._flatten(matches))
            else:
                # Handle atom descriptor
                idx_set.append(entry)

        return self.select_atoms_idx(idx_set, inplace=inplace)

    def select_atoms_idx(
        self,
        selection: AtomDescriptor | Sequence[AtomDescriptor] | None = None,
        extend_selection: bool = False,
        inplace: bool = False,
    ) -> Self:
        """Function to update the selection of atoms by specifying atom indices directly.

        Args:
            selection (AtomDescriptor | Sequence[AtomDescriptor] | None, optional): Either a single atom selector or a sequence of atom selectors. Defaults to None, which means that all available atoms will be considered.
            extend_selection (bool, optional): If set to True, the selection will be extended by the atoms denoted by `selection`. Otherwise, the new selection will be the intersection between the old selection and `selection`. Defaults to False.
            inplace (bool, optional): Whether to update this selection in-place. Defaults to False.

        Returns:
            StructureSelection: The updated selection.
        """
        if isinstance(selection, AtomDescriptor):
            selection_set = {selection}
        elif selection is None:
            selection_set = self.atoms.copy()
        else:
            selection_set = set(selection)

        get_only_available = self.atoms.intersection(selection_set)
        if extend_selection:
            new_selection = self.atoms_selected.union(get_only_available)
        else:
            new_selection = self.atoms_selected.intersection(get_only_available)

        return self.copy_or_update(atoms_selected=new_selection, inplace=inplace)

    def select_bonds(
        self,
        selection: str
        | Sequence[str]
        | BondDescriptor
        | Sequence[BondDescriptor]
        | Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        """Restrict the current selection of bonds by either specifying a SMARTS string (or sequence thereof) to specify substructures of
        the molecule to consider bonds in, or by providing one or more bond desciptor tuples or by providing a list of atoms that can be
        passed to the atom-based selection in `self.select_bonds_by_atoms()`.

        Args:
            selection (str | Sequence[str] | BondDescriptor | Sequence[BondDescriptor] | Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): The criterion or criteria by which to retain bonds in the selection. Defaults to None meaning that all bonds will be added back into the selection.
            inplace (bool, optional): Whether to update the selection in-place or return an updated copy. Defaults to False.

        Raises:
            ValueError: If no `self.mol` object is set and a SMARTS match is attempted

        Returns:
            Self: The updated selection object.
        """
        new_selection = set()
        if isinstance(selection, str) or isinstance(selection, tuple):
            selection_list = [selection]
        elif selection is None:
            new_selection = self.bonds.copy()
            selection_list = []
        else:
            selection_list = selection

        for entry in selection_list:
            if isinstance(entry, str):
                # Found smarts match:
                if self.mol is None:
                    raise ValueError(
                        "No Molecule set in selection. Cannot match SMARTS."
                    )

                matches = self.__match_pattern(self.mol, entry)
                from_matches = self._new_bond_selection_from_atoms(matches)
                new_selection.update(from_matches)
            elif isinstance(entry, AtomDescriptor):
                # We have an atom selection list or list thereof. Consume it and stop iteration.
                new_selection.update(
                    self._new_bond_selection_from_atoms(selection_list)
                )
                break
            elif isinstance(entry, tuple):
                new_selection.add(entry)
            else:
                # We have a sequence of sequences of atoms:
                try:
                    if isinstance(entry[0], AtomDescriptor):
                        new_selection.update(self._new_bond_selection_from_atoms(entry))
                except:
                    pass

        return self.copy_or_update(bonds_selected=new_selection, inplace=inplace)

    def select_bonds_idx(
        self,
        selection: BondDescriptor | Sequence[BondDescriptor],
        inplace: bool = False,
    ) -> Self:
        """Helper function to select bonds by specifying the explicit Bond descriptors/tuples.

        Restricts the selection further to this set.

        Args:
            selection (BondDescriptor | Sequence[BondDescriptor]): Either an individual bond selector or a sequence of bonds to select.
            inplace (bool, optional): Whether the selection should be updated in-place. Defaults to False.

        Returns:
            StructureSelection: The updated selections
        """
        new_selection = set()

        if isinstance(selection, tuple):
            if selection in self.bonds_selected:
                new_selection.add(selection)
        else:
            new_selection = self.bonds_selected.intersection(selection)

        return self.copy_or_update(bonds_selected=new_selection, inplace=inplace)

    def select_bonds_by_atoms(
        self,
        atoms: Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        """Helper function to select bonds by specifying a subset of atoms to consider for bonds between them.

        Allows provision of a single list of atoms or multiple such lists and will iterate over them as needed.

        Args:
            atoms (Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): Either a single set of atoms to keep bonds between or multiple sets within which the bonds should be kept. Defaults to None.
            inplace (bool, optional): Whether the selection should be updated in-place. Defaults to False.

        Returns:
            StructureSelection: The updated selections
        """
        new_selection = self._new_bond_selection_from_atoms(atoms)

        return self.copy_or_update(bonds_selected=new_selection, inplace=inplace)

    def _new_bond_selection_from_atoms(
        self,
        atoms: Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        consider_all: bool = False,
    ) -> set[BondDescriptor]:
        """Internal helper to get bond selection instead of directly updating the selection

        Args:
            atoms (Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional):  Either a single set of atoms to keep bonds between or multiple sets within which the bonds should be kept. Defaults to None.
            consider_all (bool, optional): Whether to use the entire set of features in the whole molecule as basis or just the selected set. Defaults to using only the currently selected set.
        Returns:
            set[BondDescriptor]: The set of bond descriptors in current selection fully covered by these atoms
        """

        basis_set = self.bonds if consider_all else self.bonds_selected
        new_selection: set[BondDescriptor]
        if atoms is None:
            new_selection = self.bonds.copy()
        else:
            new_selection = set()

            for entry in atoms:
                filter_set = None
                # Flag to allow for breaking out of the loop if the atoms array should have been used as filter instead.
                break_after = False
                if isinstance(entry, AtomDescriptor):
                    # atoms is a sequence of atoms to select from
                    filter_set = atoms
                    break_after = True
                else:
                    filter_set = entry

                for bond in basis_set:
                    if bond[0] in filter_set and bond[1] in filter_set:
                        new_selection.add(bond)

                if break_after:
                    break
        return new_selection

    def select_angles(
        self,
        selection: str
        | Sequence[str]
        | AngleDescriptor
        | Sequence[AngleDescriptor]
        | Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        """Function to restrict the angles selection by providing either providing SMARTS strings or explicit angles descriptors
        or sets of atoms between which to retain angles.

        Args:
            selection (str | Sequence[str] | AngleDescriptor | Sequence[AngleDescriptor] | Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): The criterion or criteria by which to retain angles in the selection. Defaults to None meaning that all angles will be added back into the selection.
            inplace (bool, optional): Whether the selection should be updated in-place. Defaults to False.

        Returns:
            Self: The updated selection
        """
        new_selection = set()
        if isinstance(selection, str) or isinstance(selection, tuple):
            selection_list = [selection]
        elif selection is None:
            new_selection = self.angles.copy()
            selection_list = []
        else:
            selection_list = selection

        for entry in selection_list:
            if isinstance(entry, str):
                # Found smarts match:
                matches = self.__match_pattern(self.mol, entry)
                new_selection.update(self._new_angle_selection_from_atoms(matches))
            elif isinstance(entry, AtomDescriptor):
                # We have an atom selection list or list thereof. Consume it and stop iteration.
                new_selection.update(
                    self._new_angle_selection_from_atoms(selection_list)
                )
                break
            elif isinstance(entry, tuple):
                new_selection.add(entry)
            else:
                # We have a sequence of sequences of atoms:
                try:
                    if isinstance(entry[0], AtomDescriptor):
                        new_selection.update(
                            self._new_angle_selection_from_atoms(entry)
                        )
                except:
                    pass

        return self.copy_or_update(angles_selected=new_selection, inplace=inplace)

    def select_angles_by_atoms(
        self,
        atoms: Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        """Helper function to select angles by specifying a subset of atoms to consider angles between.

        Allows provision of a single list of atoms or multiple such lists and will iterate over them as needed.

        Args:
            atoms (Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): Either a single set of atoms to keep angles between or multiple sets within which the bonds should be kept. Defaults to None.
            inplace (bool, optional): Whether the selection should be updated in-place. Defaults to False.

        Returns:
            StructureSelection: The updated selection.
        """
        new_selection = self._new_angle_selection_from_atoms(atoms)

        return self.copy_or_update(angles_selected=new_selection, inplace=inplace)

    def _new_angle_selection_from_atoms(
        self,
        atoms: Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        consider_all: bool = False,
    ) -> set[AngleDescriptor]:
        """Internal helper to get angle selection instead of directly updating the selection

        Args:
            atoms (Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): Either a single set of atoms to keep angles between or multiple sets within which the bonds should be kept. Defaults to None.
            consider_all (bool, optional): Whether to use the entire set of features in the whole molecule as basis or just the selected set. Defaults to using only the currently selected set.
        Returns:
            set[AngleDescriptor]: The set of agnle descriptors in current selection fully covered by these `atoms`.
        """

        basis_set = self.angles if consider_all else self.angles_selected
        new_selection: set[AngleDescriptor]
        if atoms is None:
            new_selection = self.angles.copy()
        else:
            new_selection = set()

            for entry in atoms:
                filter_set = None
                # Flag to allow for breaking out of the loop if the atoms array should have been used as filter instead.
                break_after = False
                if isinstance(entry, AtomDescriptor):
                    # atoms is a sequence of atoms to select from
                    filter_set = atoms
                    break_after = True
                else:
                    filter_set = entry

                for angle in basis_set:
                    if (
                        angle[0] in filter_set
                        and angle[1] in filter_set
                        and angle[2] in filter_set
                    ):
                        new_selection.add(angle)

                if break_after:
                    break
        return new_selection

    def select_dihedrals(
        self,
        selection: str
        | Sequence[str]
        | DihedralDescriptor
        | Sequence[DihedralDescriptor]
        | Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        """Function to restrict the dihedral selection by providing either providing SMARTS strings or explicit dihedral descriptors or sets of atoms between which to retain dihedrals.

        Args:
            selection (str | Sequence[str] | DihedralDescriptor | Sequence[DihedralDescriptor] | Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): The criterion or criteria by which to retain dihedrals in the selection. Defaults to None meaning that all dihedrals will be added back into the selection.
            inplace (bool, optional): Whether the selection should be updated in-place. Defaults to False.

        Returns:
            Self: The updated selection
        """
        new_selection = set()
        if isinstance(selection, str) or isinstance(selection, tuple):
            selection_list = [selection]
        elif selection is None:
            new_selection = self.dihedrals.copy()
            selection_list = []
        else:
            selection_list = selection

        for entry in selection_list:
            if isinstance(entry, str):
                # Found smarts match:
                matches = self.__match_pattern(self.mol, entry)
                new_selection.update(self._new_dihedral_selection_from_atoms(matches))
            elif isinstance(entry, AtomDescriptor):
                # We have an atom selection list or list thereof. Consume it and stop iteration.
                new_selection.update(
                    self._new_dihedral_selection_from_atoms(selection_list)
                )
                break
            elif isinstance(entry, tuple):
                new_selection.add(entry)
            else:
                # We have a sequence of sequences of atoms:
                try:
                    if isinstance(entry[0], AtomDescriptor):
                        new_selection.update(
                            self._new_dihedral_selection_from_atoms(entry)
                        )
                except:
                    pass

        return self.copy_or_update(dihedrals_selected=new_selection, inplace=inplace)

    def _new_dihedral_selection_from_atoms(
        self,
        atoms: Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        consider_all: bool = False,
    ) -> set[DihedralDescriptor]:
        """Internal helper to get dihedral selection instead of directly updating the selection

        Args:
            atoms (Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): Either a single set of atoms to keep dihedrals between or multiple sets within which the dihedrals should be kept. Defaults to None.
            consider_all (bool, optional): Whether to use the entire set of features in the whole molecule as basis or just the selected set. Defaults to using only the currently selected set.
        Returns:
            set[DihedralDescriptor]: The set of dihedral descriptors in current selection fully covered by these atoms
        """
        basis_set = self.dihedrals if consider_all else self.dihedrals_selected
        new_selection: set[DihedralDescriptor]
        if atoms is None:
            new_selection = self.dihedrals.copy()
        else:
            new_selection = set()

            for entry in atoms:
                filter_set = None
                # Flag to allow for breaking out of the loop if the atoms array should have been used as filter instead.
                break_after = False
                if isinstance(entry, AtomDescriptor):
                    # atoms is a sequence of atoms to select from
                    filter_set = atoms
                    break_after = True
                else:
                    filter_set = entry

                for dihedral in basis_set:
                    if all(x in filter_set for x in dihedral):
                        new_selection.add(dihedral)

                if break_after:
                    break
        return new_selection

    def select_pyramids(
        self,
        selection: str
        | Sequence[str]
        | PyramidsDescriptor
        | Sequence[PyramidsDescriptor]
        | Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        """Function to restrict the pyramid selection by providing either SMARTS strings or explicit pyramids descriptors or sets of atoms between which to retain pyramids.

        Args:
            selection (str | Sequence[str] | PyramidsDescriptor | Sequence[PyramidsDescriptor] | Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): The criterion or criteria by which to retain pyramids in the selection. Defaults to None meaning that all pyramids will be added back into the selection.
            inplace (bool, optional): Whether the selection should be updated in-place. Defaults to False.

        Returns:
            Self: The updated selection
        """
        new_selection = set()
        if isinstance(selection, str) or isinstance(selection, tuple):
            selection_list = [selection]
        elif selection is None:
            new_selection = self.pyramids.copy()
            selection_list = []
        else:
            selection_list = selection

        for entry in selection_list:
            if isinstance(entry, str):
                # Found smarts match:
                matches = self.__match_pattern(self.mol, entry)
                new_selection.update(self._new_pyramid_selection_from_atoms(matches))
            elif isinstance(entry, AtomDescriptor):
                # We have an atom selection list or list thereof. Consume it and stop iteration.
                new_selection.update(
                    self._new_pyramid_selection_from_atoms(selection_list)
                )
                break
            elif isinstance(entry, tuple):
                new_selection.add(entry)
            else:
                # We have a sequence of sequences of atoms:
                try:
                    if isinstance(entry[0], AtomDescriptor):
                        new_selection.update(
                            self._new_pyramid_selection_from_atoms(entry)
                        )
                except:
                    pass

        return self.copy_or_update(pyramids_selected=new_selection, inplace=inplace)

    def _new_pyramid_selection_from_atoms(
        self,
        atoms: Sequence[AtomDescriptor]
        | Sequence[Sequence[AtomDescriptor]]
        | None = None,
        consider_all: bool = False,
    ) -> set[PyramidsDescriptor]:
        """Internal helper to get pyramids selection instead of directly updating the selection

        Args:
            atoms (Sequence[AtomDescriptor] | Sequence[Sequence[AtomDescriptor]] | None, optional): Either a single set of atoms to keep pyramids between or multiple sets within which the pyramids should be kept. Defaults to None.
            consider_all (bool, optional): Whether to use the entire set of features in the whole molecule as basis or just the selected set. Defaults to using only the currently selected set.
        Returns:
            set[Pyramid]: The set of dihedral descriptors in current selection fully covered by these atoms
        """
        basis_set = self.pyramids if consider_all else self.pyramids_selected
        new_selection: set[PyramidsDescriptor]
        if atoms is None:
            new_selection = self.pyramids.copy()
        else:
            new_selection = set()

            for entry in atoms:
                filter_set: Sequence[AtomDescriptor]
                # Flag to allow for breaking out of the loop if the atoms array should have been used as filter instead.
                break_after = False
                if isinstance(entry, tuple):
                    # atoms is a sequence of atoms to select from
                    filter_set = atoms
                    break_after = True
                else:
                    filter_set = entry

                for pyramid in basis_set:
                    x, (sides) = pyramid
                    if x in filter_set and all(x in filter_set for x in sides):
                        new_selection.add(pyramid)

                if break_after:
                    break
        return new_selection

    def select_bats(
        self,
        smarts: str | Sequence[str] | None = None,
        idxs: FeatureDescriptor | Sequence[FeatureDescriptor] | None = None,
        mode: Literal['intersect', 'ext', 'sub'] = 'intersect',
        inplace: bool = False,
    ) -> Self:
        """Update entire selection on this molecule to a subset of available atoms, bonds, angles or dihedrals.

        Updates can be requested by providing smarts strings or by providing specific ids of features, where the feature type will be
        determined based on the length of the tuple.


        Args:
            smarts (str | Sequence[str] | None, optional): One or more smarts to identify subsets of the molecule and the features therein. Defaults to None.
            idxs (FeatureDescriptor | Sequence[FeatureDescriptor] | None, optional): Either a single tuple or a sequence of tuples to use for the update. Will be assigned based on the length of the tuple. Defaults to None.
            mode (Literal['intersect', 'ext', 'sub'], optional): The mode for the update. The new selection can either be the intersection of the current selection and the features covered by the new update set, it can be extended to contain the new update set ('ext') or the new update set can be removed from the current selection (`sub`). Defaults to 'intersect'.
            inplace (bool, optional): Whether the selection should be updated in-place. Defaults to False.

        Returns:
            Self: the updated selection
        """
        new_atoms_selection = set()
        new_bonds_selection = set()
        new_angles_selection = set()
        new_dihedrals_selection = set()
        new_pyramids_selection = set()

        consider_all_flag = mode == 'ext'

        if smarts is None and idxs is None:
            logging.warning("No selection criteria provided.")
            return self

        if smarts is not None:
            if isinstance(smarts, str):
                smarts = [smarts]
            for smarts_string in smarts:
                matches = self.__match_pattern(self.mol, smarts_string)

                # TODO: FIXME: Propagate the extension mode further down.
                # currently, we do not add bonds back if they are not in the selection.
                # TODO: FIXME: We should add a flag to iterate over the entire set of features instead of just the selected set.

                new_atoms_selection.update(self._flatten(matches))
                new_bonds_selection.update(
                    self._new_bond_selection_from_atoms(
                        matches, consider_all=consider_all_flag
                    )
                )
                new_angles_selection.update(
                    self._new_angle_selection_from_atoms(
                        matches, consider_all=consider_all_flag
                    )
                )
                new_dihedrals_selection.update(
                    self._new_dihedral_selection_from_atoms(
                        matches, consider_all=consider_all_flag
                    )
                )
                new_pyramids_selection.update(
                    self._new_pyramid_selection_from_atoms(
                        matches, consider_all=consider_all_flag
                    )
                )

        if idxs is not None:
            if isinstance(idxs, AtomDescriptor) or isinstance(idxs, tuple):
                idxs = [idxs]

            for idx in idxs:
                if isinstance(idx, AtomDescriptor):
                    new_atoms_selection.add(idx)
                else:
                    tuple_len = len(idx)
                    if tuple_len == 2:
                        if isinstance(idx[1], tuple):
                            new_pyramids_selection.add(idx)
                        else:
                            new_bonds_selection.add(idx)
                    elif tuple_len == 3:
                        new_angles_selection.add(idx)
                    elif tuple_len == 4:
                        new_dihedrals_selection.add(idx)

        if mode == 'intersect':
            new_atoms_selection = new_atoms_selection.intersection(self.atoms_selected)
            new_bonds_selection = new_bonds_selection.intersection(self.bonds_selected)
            new_angles_selection = new_angles_selection.intersection(
                self.angles_selected
            )
            new_dihedrals_selection = new_dihedrals_selection.intersection(
                self.dihedrals_selected
            )
            new_pyramids_selection = new_pyramids_selection.intersection(
                self.pyramids_selected
            )
        elif mode == 'ext':
            new_atoms_selection.update(self.atoms_selected)
            new_bonds_selection.update(self.bonds_selected)
            new_angles_selection.update(self.angles_selected)
            new_dihedrals_selection.update(self.dihedrals_selected)
            new_pyramids_selection.difference_update(self.pyramids_selected)
        elif mode == 'sub':
            new_atoms_selection.difference_update(self.atoms_selected)
            new_bonds_selection.difference_update(self.bonds_selected)
            new_angles_selection.difference_update(self.angles_selected)
            new_dihedrals_selection.difference_update(self.dihedrals_selected)
            new_pyramids_selection.difference_update(self.pyramids_selected)

        return self.copy_or_update(
            atoms_selected=new_atoms_selection,
            bonds_selected=new_bonds_selection,
            angles_selected=new_angles_selection,
            dihedrals_selected=new_dihedrals_selection,
            pyramids_selected=new_pyramids_selection,
            inplace=inplace,
        )

    @staticmethod
    def __match_pattern(
        mol: Mol | None, smarts: str
    ) -> list[Sequence[AtomDescriptor]] | None:
        """
        Find all substructure matches of a SMARTS pattern in a molecule.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol|None
            RDKit molecule object.
        smarts : str
            SMARTS pattern to search for.

        Returns
        -------
        list of Sequence of AtomDescriptor
            Each Entry contains all atom indices corresponding to one match of the SMARTS pattern.
            Returns an empty list if no match is found.
        None
            if the provided SMARTS string was invalid.

        Raises:
        -------
        ValueError: If no `mol` object is provided. Cannot match if not molecule object provided.
        """
        if mol is None:
            raise ValueError("No Molecule set in selection. Cannot match SMARTS.")

        pattern = rdkit.Chem.MolFromSmarts(smarts)

        if pattern is None:
            # TODO: FIXME: Raise ValueError instead?
            logging.info(
                f"Invalid SMARTS '{smarts}'. Falling back to full reference molecule."
            )
            matches = None
        else:
            matches = list(mol.GetSubstructMatches(pattern))

        return matches

    @staticmethod
    def _flatten(obj) -> Iterator:
        """Helper functiont to flatten nested lists

        Args:
            obj (list|Any): A potentially nested set of lists.

        Yields:
            Iterator[Any]: The iterator to iterate over all entries in the flattened list.
        """
        if isinstance(obj, list):
            for item in obj:
                yield from StructureSelection._flatten(item)
        else:
            yield obj

    def draw(
        self,
        flag_level: FeatureLevelOptions = 'bonds',
        highlight_color: tuple[float, float, float] | str = st_yellow,
        width=300,
        height=300,
    ) -> SVG:
        """Helper function to allow visualization of the structure represented in this selection.

        Args:
            flag_level (FeatureLevelOptions, optional): Currently unused. Defaults to 'bonds'.
            highlight_color (tuple[float, float, float] | str, optional): Color to use for highlights of the active parts. Defaults to st_yellow.
            width (int, optional): Width of the figure. Defaults to 300.
            height (int, optional): Height of the figure. Defaults to 300.

        Returns:
            SVG: _description_
        """
        from rdkit.Chem.Draw import rdMolDraw2D

        if isinstance(highlight_color, str):
            highlight_color = tuple(hex2rgb(highlight_color))

        # draw molecule with highlights
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.drawOptions().fillHighlights = True
        drawer.drawOptions().addAtomIndices = True
        drawer.drawOptions().setHighlightColour(highlight_color)
        drawer.drawOptions().clearBackground = False

        active_bonds = self.__get_active_bonds(flag_level=flag_level)
        active_atoms = self.__get_active_atoms(flag_level=flag_level)

        if len(active_bonds) == 0:
            active_bonds = None
        else:
            active_bonds = self.__bond_descriptor_to_mol_index(active_bonds)

        # print(f"{drawer=}")
        # print(f"{self.mol=}")
        # print(f"{active_atoms=}")
        # print(f"{active_bonds=}")

        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            self.mol,
            highlightAtoms=active_atoms,
            highlightBonds=active_bonds,
        )
        # drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        img_text = drawer.GetDrawingText()
        img = SVG(img_text)

        return img

    def __get_active_atoms(
        self, flag_level: FeatureLevelOptions = 'atoms'
    ) -> list[AtomDescriptor]:
        """Helper function to get the list of atoms involved in the selection at feature level `flag_level`.

        Available levels are:
        - 1 | 'atoms': atoms selected
        - 2 | 'bonds': bonds selected
        - 3 | 'angles': angles selected
        - 4 | 'dihedrals': dihedrals selected
        - 5 | 'pyramids': pyramids selected

        Args:
            flag_level (Literal[1, 2, 3, 4], optional): The level of selection that should be used for finding all involved bonds. Defaults to 1.

        Returns:
            list[AtomDescriptor]: The list of atom indices involved in the selection at the desired feature level.s
        """
        if flag_level == 1 or flag_level == 'atoms':
            return list(self.atoms_selected)
        elif flag_level == 2 or flag_level == 'bonds':
            return list(set([at for bon in self.bonds_selected for at in bon]))
        elif flag_level == 3 or flag_level == 'angles':
            return list(set([at for an in self.angles_selected for at in an]))
        elif flag_level == 4 or flag_level == 'dihedrals':
            return list(set([at for dih in self.dihedrals_selected for at in dih]))
        elif flag_level == 5 or flag_level == 'pyramids':
            return list(
                set(
                    [at for (x, sides) in self.pyramids_selected for at in sides]
                    + [x for (x, sides) in self.pyramids_selected]
                )
            )

        return []

    def __get_active_bonds(
        self, flag_level: FeatureLevelOptions = 'bonds'
    ) -> list[BondDescriptor]:
        """Get the list of active bonds at a certain level of selection, i.e. in bonds, in angles, in dihedrals, or in pyramids.

        Args:
            flag_level (FeatureLevelOptions, optional): The level of selection that should be used for finding all involved bonds. Needs to be at least 2 (`bonds`) to yield any bonds. Defaults to 2.

        Returns:
            list[BondDescriptor]: The list of involved bond descriptors at that feature level.
        """
        if flag_level == 1 or flag_level == 'atoms':
            return []
        elif flag_level == 2 or flag_level == 'bonds':
            return list(self.bonds_selected)

        if flag_level == 3 or flag_level == 'angles':
            base_set = self.angles_selected
        elif flag_level == 4 or flag_level == 'dihedrals':
            base_set = self.dihedrals_selected
        else:
            bond_set = set()
            for x, sides in self.pyramids_selected:
                bond_set.add((x, sides[0]))
                bond_set.add((x, sides[1]))
                bond_set.add((x, sides[2]))
            return list(bond_set)

        bond_set = set()

        for tup in base_set:
            # for i in range(len(tup) - 1):
            #     a = tup[i]
            #     b = tup[i + 1]
            #     bond_set.add((int(a), int(b)))
            for a in tup:
                for b in tup:
                    if a != b and (a, b) in self.bonds:
                        bond_set.add((int(a), int(b)))

        return list(bond_set)

    def __bond_descriptor_to_mol_index(
        self, bond_descriptors: list[BondDescriptor]
    ) -> list[int]:
        """Helper function to translate a list of Bond descriptors into RDKit bond indices.

        Args:
            bond_descriptors (list[BondDescriptor]): The list of BondDescriptor tuples that we want to translate into RDKit mol internal bond indices.

        Returns:
            list[int]: The mapped list of RDKit self.mol internal bond indices.

        Raise:
            AssertionError: if self.mol is None, no mapping can be performed.
        """
        res: list[int] = []
        assert (
            self.mol is not None
        ), 'No molecule set for this selection. Cannot resolve bond ids.'
        for entry in bond_descriptors:
            bond = self.mol.GetBondBetweenAtoms(entry[0], entry[1])
            res.append(bond.GetIdx())

        return res

    @staticmethod
    def _to_feature_level_str(ft: FeatureLevelOptions) -> FeatureLevelType:
        if isinstance(ft, str):
            assert (
                ft in FEATURE_LEVELS
            ), f"Unknown feature level: {ft} supported are only {FEATURE_LEVELS}"
            return ft
        elif isinstance(ft, int):
            ft_offset = max(ft - 1, 0)
            assert ft_offset < len(
                FEATURE_LEVELS
            ), f"Unknown feature level: {ft} supported are only 1-{len(FEATURE_LEVELS)}"
            return FEATURE_LEVELS[ft_offset]
        else:
            raise ValueError(
                f"Unknown feature level {ft}. Type was {type(ft)} supported are only int or str"
            )
