from dataclasses import dataclass
from typing import Self, TypeAlias
import xarray as xr


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

    atoms: list[tuple[ActiveFlag, AtomDescriptor]]
    bonds: list[tuple[ActiveFlag, BondDescriptor]]
    angles: list[tuple[ActiveFlag, AngleDescriptor]]
    dihedrals: list[tuple[ActiveFlag, DihedralDescriptor]]

    def copy_or_update(
        self,
        atoms: list[tuple[ActiveFlag, AtomDescriptor]] | None = None,
        bonds: list[tuple[ActiveFlag, BondDescriptor]] | None = None,
        angles: list[tuple[ActiveFlag, AngleDescriptor]] | None = None,
        dihedrals: list[tuple[ActiveFlag, DihedralDescriptor]] | None = None,
        inplace: bool = False,
    ) -> Self:
        """Function to create a copy with replaced member values.

        Meant as a helper for the `Frozen` logic of the selection, i.e. method calls return a new instance
        instead of updating the existing instance.

        Args:
            atoms (list[tuple[ActiveFlag, AtomDescriptor]], optional): The list of new atom flags. Defaults to None.
            bonds (list[tuple[ActiveFlag, BondDescriptor]], optional): The list of new Bond flags. Defaults to None.
            angles (list[tuple[ActiveFlag, AngleDescriptor]] | None, optional): List of angles and a flag whether they are currently part of the selection. Defaults to None.
            dihedrals (list[tuple[ActiveFlag, DihedralDescriptor]] | None, optional): List of dihedral index tuples and a flag whether they are currently a part of the selection. Defaults to None.
            inplace (bool, optional): Flag to allow for in-place updates instead of returning a new cop. Defaults to False.

        Returns:
            StructureSelection: The selection update with the new members set. Can either be a copy if `inplace=False` or the old instance with updated members otherwise.
        """
        if inplace:
            # Update and create
            if atoms is not None:
                self.atoms = atoms
            if bonds is not None:
                self.bonds = bonds
            if angles is not None:
                self.angles = angles
            if dihedrals is not None:
                self.dihedrals = dihedrals

            return self
        else:
            if atoms is None:
                atoms = self.atoms
            if bonds is None:
                bonds = self.bonds
            if angles is None:
                angles = self.angles
            if dihedrals is None:
                dihedrals = self.dihedrals

            return type(self)(
                atoms=atoms,
                bonds=bonds,
                angles=angles,
                dihedrals=dihedrals,
            )

    @classmethod
    def init_from_dataset(cls: type[Self], dataset: xr.Dataset) -> Self:
        """Alternative constructor that creates an initial StructureSelection object from a dataset using the entire structural information in it.

        Args:
            cls (type[StructureSelection]): The type of this StructureSelection so that we can create instances of it.
            dataset (xr.Dataset): The dataset to extract the structure information out of. Must have at least `atXYZ` variable and a `atom` dimension. Ideally, an `atom` coordinate for feature selection is also provided.

        Raises:
            ValueError: If no structural information could be extracted from the dataset

        Returns:
            StructureSelection: A structure selection object initially covering all atoms and structural features.
        """
        # TODO: FIXME: Implement actual feature selection with geomatch

        atoms = list()
        bonds = list()
        angles = list()
        dihedrals = list()

        # Create an initial state selection
        return cls(
            atoms=atoms,
            bonds=bonds,
            angles=angles,
            dihedrals=dihedrals,
        )
