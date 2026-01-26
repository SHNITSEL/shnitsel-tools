from typing import Hashable, Literal, NamedTuple, TypeAlias, TypeVar
import xarray as xr

# For general analysis and data passing
AtXYZ: TypeAlias = xr.DataArray
DimName: TypeAlias = str
Frames: TypeAlias = xr.Dataset
InterState: TypeAlias = xr.Dataset
PerState: TypeAlias = xr.Dataset

# For spectra calculation
SpectraDictType: TypeAlias = dict[tuple[float, tuple[int, int]], xr.DataArray]

# Types for stacked and unstacked trajectories
Stacked: TypeAlias = xr.Dataset | xr.DataArray
Unstacked: TypeAlias = xr.Dataset | xr.DataArray

DatasetOrArray = TypeVar("DatasetOrArray", bound=xr.Dataset | xr.DataArray)

DataArrayOrVar = TypeVar("DataArrayOrVar", bound=xr.Variable | xr.DataArray)

# Types for selecting and managing states and state combinations
StateId: TypeAlias = int
StateCombination: TypeAlias = tuple[StateId, StateId]

MultiplicityLabel: TypeAlias = Literal[
    's', 'S', 'singlet', 'd', 'D', 'doublet', 't', 'T', 'triplet'
]

MultiplicityLabelValues: set[MultiplicityLabel] = {
    's',
    'S',
    'singlet',
    'd',
    'D',
    'doublet',
    't',
    'T',
    'triplet',
}


class StateInfo(NamedTuple):
    id: StateId
    name: str
    multiplicity: int | None
    charge: int | None


class StateCombInfo(NamedTuple):
    ids: StateCombination
    name: str


StateTypeSpecifier = (
    Literal['s', 'S', 'd', 'D', 't', 'T', 'singlet', 'doublet', 'triplet'] | int
)

ErrorOptions = Literal["raise", "ignore"]
ErrorOptionsWithWarn = Literal["raise", "warn", "ignore"]
