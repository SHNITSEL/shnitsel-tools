from typing import Hashable, TypeAlias, TypeVar
import xarray as xr

# For general analysis and data passing
AtXYZ: TypeAlias = xr.DataArray
DimName: TypeAlias = Hashable
Frames: TypeAlias = xr.Dataset
InterState: TypeAlias = xr.Dataset
PerState: TypeAlias = xr.Dataset

# For spectra calculation
SpectraDictType: TypeAlias = dict[tuple[float, tuple[int, int]], xr.DataArray]

# Types for stacked and unstacked trajectories
Stacked: TypeAlias = xr.Dataset | xr.DataArray
Unstacked: TypeAlias = xr.Dataset | xr.DataArray

DatasetOrArray = TypeVar("DatasetOrArray", bound=xr.Dataset | xr.DataArray)
