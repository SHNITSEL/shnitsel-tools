import xarray as xr
import numpy as np
from typing import List, Tuple, Optional
from shnitsel.static.utils import convert_distance


class SOAP_creator():
    def __init__(self, species: List[str], rcut: float=5.0, nmax: int=8, lmax: int=8, sigma: float=0.1, rbf: str='gto', average: str='off'):
        self.rcut = rcut
        self.nmax = nmax
        self.lmax = lmax
        self.sigma = sigma
        self.species = species
        self.soap = None
        self.rbf = rbf
        if average not in ['off', 'inner', 'outer']:
            raise ValueError('average must be one of "off", "inner", "outer"')
        self.average = average
        
    def create_descriptors(self, data: xr.Dataset, store_results: bool=False, n_jobs: int=-1) -> List[np.ndarray]:
        from dscribe.descriptors import SOAP
        from ase import Atoms
        from ase.build import molecule
        from tqdm import tqdm
        import multiprocessing

        # Create the SOAP descriptor
        self.soap = SOAP(
            r_cut=self.rcut,
            n_max=self.nmax,
            l_max=self.lmax,
            sigma=self.sigma,
            rbf=self.rbf,
            periodic=False,
            species=self.species,
            sparse=False,
            average=self.average
        )
        
        # check for correct shapes
        if data.positions.dims != ('frame', 'atom', 'direction'):
            raise ValueError('DataArray must have dimensions (frame, atom, direction), got {}'.format(data.dims))
        
        # convert coordinates to dscribe default Angstrom
        positions = convert_distance(data.positions, 'angstrom')
        symbols = data.symbols.values
        
        n_jobs = n_jobs if n_jobs != -1 else multiprocessing.cpu_count()
        print(f'Creating SOAP descriptors using {n_jobs} cores')
        
        # convert to ase for compatibility with dscribe
        systems = []
        for i in range(len(positions)):
            atoms = Atoms(
                positions=positions[i].values,
                symbols=symbols
            )
            systems.append((atoms,))
            
        soaps = self.soap.create_parallel(systems, n_jobs=-1, func=self.soap.create)
        return np.stack(soaps, axis=0)
    
def pca_plot(pca_res:np.array, explained_variance: List[float]=None, color_prop: xr.DataArray=None, fontsize: int=15, state: str=None):
    
    import matplotlib.pyplot as plt

    n_components = pca_res.shape[1]
    
    x_min, x_max = np.min(pca_res[:, 0]), np.max(pca_res[:, 0])
    y_min, y_max = np.min(pca_res[:, 1]), np.max(pca_res[:, 1])
    
    params = {
        'legend.fontsize': fontsize,
        'axes.labelsize': fontsize-2,
        'axes.titlesize': fontsize,
        'axes.titlepad': 20,
        'xtick.labelsize': fontsize-5,
        'ytick.labelsize': fontsize-5,
        }
    plt.rcParams.update(params)
    
    try:
        key = list(color_prop.coords.keys())[0]
        val = color_prop.coords[key].values
    except:
        key = None
        val = None  

    if n_components == 3:
        z_min, z_max = np.min(pca_res[:, 2]), np.max(pca_res[:, 2])
        fig = plt.figure(figsize=(9,5))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(pca_res[:,0], pca_res[:,1], pca_res[:,2], c=color_prop, cmap='viridis', s=1) if color_prop is not None else ax.scatter(pca_res[:,0], pca_res[:,1], pca_res[:,2], s=5)
        ax.set_xlabel(f'PC1 ({explained_variance[0]*100:.2f}% var.)')
        ax.set_ylabel(f'PC2 ({explained_variance[1]*100:.2f}% var)')
        ax.set_zlabel(f'PC3 ({explained_variance[2]*100:.2f}% var)')
        plt.colorbar(scatter, label = f'{val} {color_prop.name} [{color_prop.attrs.get("unit", "unknown")}]') if color_prop is not None else None
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        plt.subplots_adjust(left=0.1, right=2.5, top=0.9, bottom=0.1)
    elif n_components == 2:
        fig, ax = plt.subplots(figsize=(5,5))
        scatter = ax.scatter(pca_res[:,0], pca_res[:,1], c=color_prop, cmap='viridis', s=5) if color_prop is not None else ax.scatter(pca_res[:,0], pca_res[:,1], s=5)
        ax.set_xlabel(f'PC1 ({explained_variance[0]*100:.2f}%)')
        ax.set_ylabel(f'PC2 ({explained_variance[1]*100:.2f}%)')
        plt.colorbar(scatter, label = f'{val} {color_prop.name} [{color_prop.attrs.get("unit", "unknown")}]') if color_prop is not None else None
    else:
        raise ValueError('n_components must be 2 or 3')

    plt.tight_layout()
    return pca_res, explained_variance, fig
    
    
    
    
    
    
        
        
    
