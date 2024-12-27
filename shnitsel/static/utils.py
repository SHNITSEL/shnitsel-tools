import xarray as xr
import numpy as np
from shnitsel.dynamic.postprocess import Converter

def frobenius_norm(prop_arr: xr.DataArray) -> xr.DataArray:
    return xr.apply_ufunc(
        np.linalg.norm,
        prop_arr,
        input_core_dims=[prop_arr.dims[-2:]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[prop_arr.dtype],
        keep_attrs=True,
        kwargs={'ord': 'fro', 'axis': (-2, -1)}
    )
    
distance_converter = Converter('positions', dict(
    bohr=1.0,
    angstrom=0.529177249)
)



