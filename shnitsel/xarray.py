import xarray as xr

from . import _accessors
from ._generated_accessors import DataArrayAccessor, DatasetAccessor
from . import _state

from ._accessors import DataTreeAccessor

xr.register_dataarray_accessor(_accessors.DATAARRAY_ACCESSOR_NAME)(DataArrayAccessor)
_state.DATAARRAY_ACCESSOR_NAME = _accessors.DATAARRAY_ACCESSOR_NAME
_state.DATAARRAY_ACCESSOR_REGISTERED = True

xr.register_dataset_accessor(_accessors.DATASET_ACCESSOR_NAME)(DatasetAccessor)
_state.DATASET_ACCESSOR_NAME = _accessors.DATASET_ACCESSOR_NAME
_state.DATASET_ACCESSOR_REGISTERED = True

xr.register_datatree_accessor(_accessors.DATATREE_ACCESSOR_NAME)(DataTreeAccessor)
_state.DATATREE_ACCESSOR_NAME = _accessors.DATATREE_ACCESSOR_NAME
_state.DATATREE_ACCESSOR_REGISTERED = True