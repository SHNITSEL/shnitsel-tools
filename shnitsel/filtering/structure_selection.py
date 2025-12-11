from dataclasses import dataclass
import logging
from typing import Literal, Self, Sequence, TypeAlias
import rdkit
import xarray as xr

from rdkit.Chem.rdchem import Mol

from shnitsel.bridges import to_mol


AtomDescriptor: TypeAlias = int
BondDescriptor: TypeAlias = tuple[int, int]
AngleDescriptor: TypeAlias = tuple[int, int, int]
DihedralDescriptor: TypeAlias = tuple[int, int, int, int]

FeatureDescriptor: TypeAlias = (
    AtomDescriptor | BondDescriptor | AngleDescriptor | DihedralDescriptor
)

FeatureList: TypeAlias = list[FeatureDescriptor]

ActiveFlag: TypeAlias = bool


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
            )

    @classmethod
    def init_from_dataset(
        cls: type[Self],
        dataset: xr.Dataset,
        default_selection: list[Literal['atoms', 'bonds', 'angles', 'dihedrals']] = [
            'atoms',
            'bonds',
        ],
    ) -> Self:
        """Alternative constructor that creates an initial StructureSelection object from a dataset using the entire structural information in it.

        Args:
            cls (type[StructureSelection]): The type of this StructureSelection so that we can create instances of it.
            dataset (xr.Dataset): The dataset to extract the structure information out of.
                Must have at least `atXYZ` variable and a `atom` dimension.
                Ideally, an `atom` coordinate for feature selection is also provided.
                Should only represent a single frame of data.
            default_selection (list[Literal['atoms', 'bonds', 'angles', 'dihedrals']], optional): List of features to activate as selected by default. Defaults to [ 'atoms', 'bonds', ].

        Raises:
            ValueError: If no structural information could be extracted from the dataset

        Returns:
            StructureSelection: A structure selection object initially covering all atoms and structural features.
        """
        filtered_dataset = dataset.squeeze()

        if 'frame' in filtered_dataset or 'time' in filtered_dataset:
            raise ValueError(
                "The dataset should not contain frame or data but represent a single frame of data."
            )

        # TODO: FIXME: Consider the charges needing to be set from the dataset settings.s
        mol = to_mol(
            filtered_dataset.atXYZ,
            to2D=False,
        )
        # Create an initial state selection
        return cls.init_from_mol(mol, default_selection=default_selection)

    @classmethod
    def init_from_mol(
        cls: type[Self],
        mol: Mol,
        default_selection: list[Literal['atoms', 'bonds', 'angles', 'dihedrals']] = [
            'atoms',
            'bonds',
        ],
    ) -> Self:
        """Alternative constructor that creates an initial StructureSelection object from an RDKit Mol object

        Args:
            cls (type[StructureSelection]): The type of this StructureSelection so that we can create instances of it.
            mol (rdkit.rdchem.Mol): The RDKit Mol object to extract all initial structural information out of
            default_selection (list[Literal['atoms', 'bonds', 'angles', 'dihedrals']], optional): List of features to activate as selected by default. Defaults to [ 'atoms', 'bonds', ].

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

        are_atoms_selected = 'atom' in default_selection
        are_bonds_selected = 'bonds' in default_selection
        are_angles_selected = 'angles' in default_selection
        are_dihedrals_selected = 'dihedrals' in default_selection

        for atom in mol.GetAtoms():
            atomid = (atom.GetIdx(),)
            atom_type = atom.GetSymbol()
            atoms.add(atomid)
            atoms_types[atomid] = atom_type

        if are_atoms_selected:
            atoms_selected.update(atoms)

        for bond in mol.GetBonds():
            beginIdx = (bond.GetBeginAtomIdx(),)
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
        )

    def select_atoms(
        self,
        selection: str
        | Sequence[str]
        | AtomDescriptor
        | Sequence[AtomDescriptor]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        return self.copy_or_update(inplace=inplace)

    def select_bonds(
        self,
        selection: str
        | Sequence[str]
        | BondDescriptor
        | Sequence[BondDescriptor]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        return self.copy_or_update(inplace=inplace)

    def select_angles(
        self,
        selection: str
        | Sequence[str]
        | AngleDescriptor
        | Sequence[AngleDescriptor]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        return self.copy_or_update(inplace=inplace)

    def select_dihedrals(
        self,
        selection: str
        | Sequence[str]
        | DihedralDescriptor
        | Sequence[DihedralDescriptor]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        return self.copy_or_update(inplace=inplace)

    def select_bats(
        self,
        selection: str
        | Sequence[str]
        | DihedralDescriptor
        | Sequence[DihedralDescriptor]
        | None = None,
        inplace: bool = False,
    ) -> Self:
        return self.copy_or_update(inplace=inplace)

    @staticmethod
    def __match_pattern(mol: Mol, smarts: str) -> list[Sequence[AtomDescriptor]] | None:
        """
        Find all substructure matches of a SMARTS pattern in a molecule.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
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
        """
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
