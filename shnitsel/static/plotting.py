import matplotlib.pyplot as plt
import xarray as xr
from shnitsel.dynamic.postprocess import convert_energy, norm
from .utils import frobenius_norm
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from typing import List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def plot_energy_hist(data: xr.DataArray, num_bins: int=100, limit: tuple=None, unit: str='eV', cbar_shrink: float=0.9):
    energies = data.energy if 'energy' in data else data
    energies = convert_energy(energies, unit)
    energies -= energies.min().values
    energies = energies.transpose()
    titles = data.state.values
    axis_label = f'Energy relative to minimum S0 [{unit}]'

    fig, axs = plt.subplots(energies.shape[0], figsize=(7.0, energies.shape[0]))

    max_value = float(energies.max().values)
    x_range = (0, max_value) if limit is None else limit

    counts_array = np.zeros((energies.shape[0], num_bins), dtype=int)

    for idx, energy_slice in enumerate(energies):
        counts, _ = np.histogram(energy_slice, bins=num_bins, range=x_range, density=False)
        counts_array[idx] = counts

    vmin, vmax = counts_array.min(), counts_array.max()
    cmap = plt.get_cmap('viridis')
    norm = Normalize(vmin=vmin, vmax=vmax)

    for idx, energy_slice in enumerate(energies):
        counts = counts_array[idx]
        _, bins = np.histogram(energy_slice, bins=num_bins, range=x_range, density=False)
        
        for i in range(len(bins) - 1):
            axs[idx].fill_between(
                [bins[i], bins[i+1]],
                1, 0,
                color=cmap(norm(counts[i]))
            )
        axs[idx].get_yaxis().set_visible(False)
        axs[idx].set_title(titles[idx])

    axs[-1].set_xlabel(axis_label)
    plt.tight_layout()

    cbar = plt.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        ax=axs,
        label='Frequency',
        shrink=cbar_shrink
    )
    
    return fig, axs

def plot_hist(
    prop_arr: xr.DataArray,
    num_bins: int = 100,
    limit: Optional[Tuple[float, float]] = None,
    scale_type: str = 'linear',
    combine_states: bool = False,
    states: Optional[List[str]] = None
) -> Optional[Tuple[plt.Figure, np.ndarray]]:
    """
    Plot histograms for a specified property in an xarray DataArray.

    Parameters
    ----------
    prop_arr : xr.DataArray
        The input data containing the property to plot.
    num_bins : int, optional
        Number of bins for the histogram (default is 100).
    limit : tuple of float, optional
        The lower and upper range of the bins. If not provided, defaults to the data range.
    scale_type : str, optional
        Scale type for the y-axis ('linear', 'log', etc.) (default is 'linear').
    combine_states : bool, optional
        If True, plot all states together (default is False).
    states : list of str, optional
        Specific states to include in the plot. Must be present in the DataArray (default is None).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plots.
    axs : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects for each histogram.
    """
    # Input validation
    if not isinstance(prop_arr, xr.DataArray):
        raise TypeError(f"'data' must be an xarray.DataArray, got {type(prop_arr)}.")

    if not isinstance(num_bins, int) or num_bins <= 0:
        raise ValueError("'num_bins' must be a positive integer.")

    if limit is not None:
        if (not isinstance(limit, tuple) or len(limit) != 2 or
            not all(isinstance(x, (int, float)) for x in limit)):
            raise ValueError("'limit' must be a tuple of two numbers (min, max).")

    valid_scales = ['linear', 'log', 'symlog', 'logit']
    if scale_type not in valid_scales:
        raise ValueError(f"'scale_type' must be one of {valid_scales}.")

    
    # Select specific states if provided   
    if states is not None:
        if prop_arr.dims[1] not in ['state', 'statecomb']:
            raise ValueError(f'DataArray must have a state or statecomb dimension if states are provided. Provided dimensions: {prop_arr.dims}')
        available_states = prop_arr[prop_arr.dims[1]].values
        
        if not all([state in available_states for state in states]):
            raise ValueError(f'States to select must be in {available_states}')
        prop_arr = prop_arr.sel({prop_arr.dims[1]: states})
        logger.info(f'Selected states: {states}')

        if not all(state in available_states for state in states):
            raise ValueError(f"States to select must be in {available_states}")
        
    dims = prop_arr.dims
    print(dims)
    name = prop_arr.name    
    unit = prop_arr.attrs.get('unit', 'unknown') 
    labels = prop_arr.coords[prop_arr.dims[1]].values.tolist()
    print(prop_arr.shape)   	
    if len(dims) == 4: # apply frobenius norm along dimensions direction and atom
        prop_arr = frobenius_norm(prop_arr)
    elif len(dims) == 3: # apply norm along dimension direction
        prop_arr = norm(prop_arr)
            
    prop_arr = prop_arr.transpose()
    print(prop_arr.shape)
    
    if name == 'socs': # flatten data for soc dimension
        combine_states = True
    
    if not combine_states:
        fig, axs = plt.subplots(prop_arr.shape[0], figsize=(7.0, prop_arr.shape[0] * 1.5), squeeze=False)
        axs = axs.flatten()
        
        for idx, (ax, prop_slice, label) in enumerate(zip(axs, prop_arr, labels)):
            ax.hist(prop_slice, bins=num_bins, range=limit, density=True)
            ax.set_title(label)
            ax.set_yscale(scale_type)

        axs[-1].set_xlabel(f'{name} [{unit}]')
        plt.tight_layout()
        return fig, axs     
    else:
        prop_arr = prop_arr.values.flatten()
        fig, ax = plt.subplots(figsize=(7.0, 3))
        ax.hist(prop_arr, bins=num_bins, range=limit, density=True)
        ax.set_title(name)
        ax.set_yscale(scale_type)
        ax.set_xlabel(f'{name} [{unit}]')
        plt.tight_layout()
        return fig, ax
    
    
def dihedral_plot(positions: xr.DataArray, dihedral_index: List[int], prop: xr.DataArray):
    """
    Plot dihedral angles.

    Parameters
    ----------
    positions : xr.DataArray
        The input positions.
    dihedral_index : list of int
        The indices of the dihedral angles to plot.
    color_prop : xr.DataArray, optional
        Optional data array for coloring.
    """
    dih_pos = positions[:, dihedral_index, :]
    p1, p2, p3, p4 = dih_pos[:, 0, :], dih_pos[:, 1, :], dih_pos[:, 2, :], dih_pos[:, 3, :]
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    angle = np.arctan2(
        np.linalg.norm(np.cross(n1, n2), axis=1),
        np.einsum('ij,ij->i', n1, n2)
    )
    angle = np.degrees(angle)
    angle[angle < 0] += 360
    angle[angle > 180] = 360 - angle[angle > 180]

    plt.scatter(angle, prop, cmap='viridis')
    plt.xlabel('Dihedral Angle (degrees)')
    plt.ylabel(f'{prop.name} relative to minimum [{prop.attrs.get("unit", "unknown")}]')
    
    
    
    
    


