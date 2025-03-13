import matplotlib.pyplot as plt
import xarray as xr
from shnitsel.dynamic.postprocess import convert_energy, norm
from .utils import frobenius_norm
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from typing import List, Tuple, Optional, Union
import logging
import rdkit
from .utils import convert_distance, array_to_xyz

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def plot_color_hist(
    data: xr.DataArray,
    num_bins: int = 100,
    limit: Optional[Tuple[float, float]] = None,
    cbar_shrink: float = 0.9,
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Plot energy histograms for each state in an xarray DataArray using colored patches.
    
    For each state, the function computes histogram counts and then uses 
    `fill_between` to draw a colored patch for each bin (spanning y=0 to y=1) 
    where the color is mapped to the frequency count.
    
    Parameters
    ----------
    data : xr.DataArray
        Input energy data. If the DataArray contains an 'energy' attribute,
        it is used; otherwise, the data itself is assumed to be the energies.
    num_bins : int, optional
        Number of bins for the histogram (default is 100).
    limit : tuple of float, optional
        The (min, max) range for the histogram bins. If None, the range is determined 
        from the data (0 to maximum energy).
    cbar_shrink : float, optional
        Shrink factor for the colorbar (default is 0.9).
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the subplots.
    axs : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects for each histogram.
    """
    
    dims = data.dims
    if dims[0] == 'frame':
        raise ValueError("DataArray can't have 'frame' as the first dimension, please reorder the dimensions.")

    # Get titles from the 'state' coordinate if available, else use default labels
    titles = data.coords[data.dims[0]].values
        
    axis_label = f'{data.name} [{data.attrs.get("units", "unknown")}]'

    # Create subplots for each state
    fig, axs = plt.subplots(data.shape[0], figsize=(7.0, data.shape[0] * 2))
    if data.shape[0] == 1:
        axs = np.array([axs])  # ensure axs is always an array

    # Determine histogram range and compute bin edges
    max_energy = float(data.max().values)
    x_range = (0, max_energy) if limit is None else limit
    bins = np.linspace(x_range[0], x_range[1], num_bins + 1)

    # Compute histogram counts for each state
    counts_array = np.zeros((data.shape[0], num_bins), dtype=int)
    for idx, energy_slice in enumerate(data):
        counts, _ = np.histogram(energy_slice, bins=bins, density=False)
        counts_array[idx] = counts

    # Create a normalization object for mapping counts to colors
    vmin, vmax = counts_array.min(), counts_array.max()
    cmap = plt.get_cmap('viridis')
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Plot using fill_between to create colored patches per bin
    for idx in range(data.shape[0]):
        counts = counts_array[idx]
        for i in range(len(bins) - 1):
            axs[idx].fill_between(
                [bins[i], bins[i+1]],
                [0, 0],
                [1, 1],
                color=cmap(norm(counts[i]))
            )
        axs[idx].set_title(titles[idx])
        #axs[idx].ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        axs[idx].get_yaxis().set_visible(False)
    
    axs[-1].set_xlabel(axis_label)
    plt.tight_layout()

    # Add a colorbar representing frequency count
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    plt.colorbar(sm, ax=axs, label='Frequency', shrink=cbar_shrink)
    
    return fig, axs


    
def plot_hist(
    prop_arr: xr.DataArray,
    num_bins: int = 100,
    limit: Optional[Tuple[float, float]] = None,
    scale_type: str = 'linear',
    combine_states: bool = False,
    states: Optional[List[str]] = None,
    color: Optional[str] = 'blue',
    font_size: int = 18,
    hist_kwargs: Optional[dict] = None,
    ax: Optional[plt.Axes] = None,
    state_dim: Optional[str] = None,
) -> Tuple[plt.Figure, Union[np.ndarray, plt.Axes]]:
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
    color : str, optional
        Color for the histogram (default is 'blue').
    font_size : int, optional
        Font size for the plot (default is 18).
    hist_kwargs : dict, optional
        Additional keyword arguments to pass to plt.hist.
    ax : matplotlib.axes.Axes, optional
        If provided (and combine_states is True), the histogram will be plotted on this axis.
    state_dim : str, optional
        Name of the dimension representing states. If None, the function will search for common state dimensions (e.g. 'state', 'statecomb').

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot(s).
    axs : matplotlib.axes.Axes or numpy.ndarray of matplotlib.axes.Axes
        The axis or axes objects for the histogram(s).
    """
    import matplotlib.pyplot as plt
    import logging
    logger = logging.getLogger(__name__)
    
    # Input validation
    if not isinstance(prop_arr, xr.DataArray):
        raise TypeError(f"'prop_arr' must be an xarray.DataArray, got {type(prop_arr)}.")
    if not isinstance(num_bins, int) or num_bins <= 0:
        raise ValueError("'num_bins' must be a positive integer.")
    if limit is not None:
        if (not isinstance(limit, tuple) or len(limit) != 2 or
            not all(isinstance(x, (int, float)) for x in limit)):
            raise ValueError("'limit' must be a tuple of two numbers (min, max).")
    valid_scales = ['linear', 'log', 'symlog', 'logit']
    if scale_type not in valid_scales:
        raise ValueError(f"'scale_type' must be one of {valid_scales}.")

    if hist_kwargs is None:
        hist_kwargs = {}

    # Select specific states if provided using a flexible state dimension
    if states is not None:
        possible_dims = [state_dim] if state_dim is not None else ['state', 'statecomb']
        found_dim = None
        for dim in prop_arr.dims:
            if dim in possible_dims:
                found_dim = dim
                break
        if found_dim is None:
            raise ValueError(f"DataArray must have one of the following state dimensions: {possible_dims}. Available dimensions: {prop_arr.dims}")
        available_states = prop_arr.coords[found_dim].values.tolist()
        if not all(state in available_states for state in states):
            raise ValueError(f"States to select must be among {available_states}")
        prop_arr = prop_arr.sel({found_dim: states})
        logger.info(f"Selected states: {states}")

    dims = prop_arr.dims
    name = prop_arr.name or 'property'
    unit = prop_arr.attrs.get('units', 'unknown')
    frobenius = False
    # Determine labels based on available state coordinates if possible
    if len(dims) > 1:
        # Use the provided state_dim if available or search for common names
        label_dim = state_dim if state_dim and state_dim in dims else next((d for d in dims if d in ['state', 'statecomb']), dims[1])
        labels = prop_arr.coords[label_dim].values.tolist()
    else:
        labels = [name]

    # Apply norm transformations if applicable
    if len(dims) == 4:  # apply Frobenius norm along last 2 dimensions (e.g., direction and atom)
        prop_arr = frobenius_norm(prop_arr)
        frobenius = True
    elif len(dims) == 3: # apply vector norm along last dimension 
        prop_arr = norm(prop_arr)

    prop_arr = prop_arr.transpose()

    # Auto-set combine_states for specific property names if needed
    if name == 'socs':
        combine_states = True

    # Update plot font size (note: this updates the global rcParams)
    plt.rcParams.update({'font.size': font_size})

    if not combine_states:
        # Create new subplots if multiple histograms are required.
        num_plots = prop_arr.shape[0]
        fig, axs = plt.subplots(
            num_plots, figsize=(7.0, num_plots * 3), squeeze=False
        )
        axs = axs.flatten()

        # Plot individual histograms
        for current_ax, prop_slice, label in zip(axs, prop_arr, labels):
            current_ax.hist(prop_slice, bins=num_bins, range=limit, density=False, color=color, **hist_kwargs)
            current_ax.set_title(label)
            current_ax.set_yscale(scale_type)
            current_ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
            current_ax.set_ylabel('Frequency')
        xlabel = f'Frobenius Norm of {name} \n [{unit}]' if frobenius else f'{name} [{unit}]'
        axs[-1].set_xlabel(xlabel)
        plt.tight_layout()
        return fig, axs
    else:
        # Combine states: flatten the data for a single histogram.
        flat_data = prop_arr.values.flatten()
        if ax is None:
            fig, current_ax = plt.subplots(figsize=(7.0, 4))
        else:
            fig = ax.figure
            current_ax = ax
        current_ax.hist(flat_data, bins=num_bins, range=limit, density=False, color=color, **hist_kwargs)
        current_ax.set_title(name)
        current_ax.set_yscale(scale_type)
        current_ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        xlabel = f'Frobenius Norm of {name} [{unit}]' if frobenius else f'{name} [{unit}]'
        current_ax.set_xlabel(xlabel)
        current_ax.set_ylabel('Frequency')
        plt.tight_layout()
        return fig, current_ax

    

def dihedral_plot(positions: xr.DataArray,
                  dihedral_indices: List[int],
                  prop: xr.DataArray,
                  color_prop: Optional[xr.DataArray] = None,
                  ax: Optional[plt.Axes] = None,
                 ) -> Tuple[plt.Figure, plt.Axes, object]:

    # Compute the dihedral angle
    dih_pos = positions[:, dihedral_indices, :]
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

    # Create or use provided axis
    if ax is None:
        fig, ax = plt.subplots()
        plot_cbar = True
    else:
        fig = ax.figure
        plot_cbar = False

    # Plot and optionally add a colorbar
    if color_prop is not None:
        scatter = ax.scatter(angle, prop, c=color_prop, s=10, cmap='viridis')
        ax.set_xlabel('Dihedral Angle [°]')
        ax.set_ylabel(f'{prop.name} [{prop.attrs.get("units", "unknown")}]')
        if plot_cbar:
            cbar = plt.colorbar(ax.collections[0], ax=ax)
            cbar.set_label(f'{color_prop.name} [{color_prop.attrs.get("units", "unknown")}]')
    else:
        scatter = ax.scatter(angle, prop, s=10)
        ax.set_xlabel('Dihedral Angle [°]')
        ax.set_ylabel(f'{prop.name}[{prop.attrs.get("units", "unknown")}]')
    
    return fig, ax, scatter



def display_atom_indices(positions, symbols, charge=0):
    xyz = array_to_xyz(positions, symbols)
    mol = rdkit.Chem.rdmolfiles.MolFromXYZBlock(xyz)
    rdkit.Chem.rdDetermineBonds.DetermineBonds(mol, charge=charge)
    rdkit.Chem.AllChem.Compute2DCoords(mol)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("atomLabel", str(atom.GetIdx())+":"+atom.GetSymbol())
        
    img = rdkit.Chem.AllChem.Draw.MolToImage(mol, size=(500, 500))
    return img
    
    
    
    


