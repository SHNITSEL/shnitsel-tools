# Tutorials

## 1 Data parsing and Conversion

`shnitsel.dynamic.parse` enables the construction of SHNITSEL datasets from the following sources:

- [SHARC](https://sharc-md.org/) output files: initial conditions (`ICOND_*` folders) or trajectories (`TRAJ_*` folders)
- [NewtonX](https://newtonx.org/) output files: trajectories 
- [ASE](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) databases

The respective tutorial notebooks are: [01_sharc2hdf5.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/01_sharc2hdf5.ipynb), 
[02_nx2hdf5.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/02_nx2hdf5.ipynb), and 
[03_ase2hdf5.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/03_ase2hdf5.ipynb)

## 2 Data Processing and Visualization

### 2.1 Static

The `shnitsel.static` module provides tools for analyzing and visualizing the static datasets. 
The tutorial notebook is: [15_stattic_module.ipynb] (https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/15_static_module.ipynb)

### 2.2 Dynamic 

- [ ] Data Sheet Summary

optional:
- [ ] Histograms (E, F, Dipoles, ...)
- [ ] E vs. t
- [ ] Absorption spectra at chosen (delay) times
- [ ] Population Plots
- [ ] PCA plot (dynamic)
