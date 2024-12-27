import xarray as xr
import numpy as np
from typing import List, Tuple, Optional
from shnitsel.static.utils import distance_converter


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
        positions = distance_converter(data.positions, 'angstrom')
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
    
def pca_and_plot(descriptors:np.array, n_components: int=2, color_prop: xr.DataArray=None):	
    
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import MinMaxScaler  
    
    #check shapes of descriptors and color_prop
    if color_prop is not None:
        if color_prop.shape[0] != descriptors.shape[0]:
            raise ValueError('color_prop and descriptors must have the same number of frames')
        if color_prop.ndim != 1:
            raise ValueError('color_prop must be a 1D array')
        prop_norm = color_prop - color_prop.min()
    if descriptors.ndim > 2:
        descriptors = descriptors.reshape(descriptors.shape[0], -1) # flatten the descriptors
        
    pca = PCA(n_components=n_components)
    scaler = MinMaxScaler()
    descriptors = scaler.fit_transform(descriptors)
    pca_res = pca.fit_transform(descriptors)
    explained_variance = pca.explained_variance_ratio_  
    print(f'Explained variance: {explained_variance}')  
    
    if n_components == 3:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(pca_res[:,0], pca_res[:,1], pca_res[:,2], c=prop_norm, cmap='viridis', s=5)
        ax.set_xlabel(f'PC1 ({explained_variance[0]*100:.2f}%)')
        ax.set_ylabel(f'PC2 ({explained_variance[1]*100:.2f}%)')
        ax.set_zlabel(f'PC3 ({explained_variance[2]*100:.2f}%)')
        plt.colorbar(scatter, label = f'{color_prop.name} relative to minimum [{color_prop.attrs.get("unit", "unknown")}]')
        plt.tight_layout()
        plt.show()
        
    elif n_components == 2:
        fig, ax = plt.subplots(figsize=(7,8))
        scatter = ax.scatter(pca_res[:,0], pca_res[:,1], c=prop_norm, cmap='viridis', s=5)
        ax.set_xlabel(f'PC1 ({explained_variance[0]*100:.2f}%)')
        ax.set_ylabel(f'PC2 ({explained_variance[1]*100:.2f}%)')
        plt.colorbar(scatter, label = f'{color_prop.name} relative to minimum [{color_prop.attrs.get("unit", "unknown")}]')
        plt.tight_layout()
        plt.show()
        
    else:
        raise ValueError('n_components must be 2 or 3')
    
    
    return pca_res, explained_variance, fig
    
    
    
    
    
    
        
        
    
