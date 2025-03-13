import xarray as xr
import numpy as np
from shnitsel.dynamic.postprocess import Converter
import rdkit


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

def vector_norm(prop_arr: xr.DataArray) -> xr.DataArray:
    return xr.apply_ufunc(
        np.linalg.norm,
        prop_arr,
        input_core_dims=[(prop_arr.dims[-1],)],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[prop_arr.dtype],
        keep_attrs=True,
        kwargs={'axis': -1}
    )
    
convert_distance = Converter('positions', dict(
    bohr=1.0,
    angstrom=0.529177249)
)

convert_force = Converter('forces', {
    'hartree/bohr': 1.0,
    'eV/angstrom': 51.42206749896734
})

def array_to_xyz(positions, symbols):
    try:
        positions = convert_distance(positions, 'angstrom')
    except:
        print("Warning: Could not convert positions to angstrom.")
    
    num_atoms = len(symbols)
    header = f"{num_atoms}\n \n"
    body_lines = []
    for symbol, pos in zip(symbols, positions):
        line = f"{symbol}    {pos[0]:.8f}    {pos[1]:.8f}    {pos[2]:.8f}"
        body_lines.append(line)
        
    body = "\n".join(body_lines)
    return header + body








