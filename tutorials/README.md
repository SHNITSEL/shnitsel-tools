# Tutorials

## Data parsing and conversion

`shnitsel-tools` can read the following:

- [SHARC](https://sharc-md.org/) output files: initial conditions (`ICOND_*` folders) or trajectories (`TRAJ_*` folders)
- [NewtonX](https://newtonx.org/) output files: trajectories 
- [PyRAI2MD](https://github.com/lopez-lab/PyRAI2MD) output files: trajectories 
- [ASE](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) databases

The respective tutorial notebooks are: [01_sharc2hdf5.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_1_sharc2hdf5.ipynb), 
[0_2_nx2hdf5.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_2_nx2hdf5.ipynb), 
[0_3_pyrai2md2hdf5.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_3_pyrai2md2hdf5.ipynb) and
[0_4_ase2hdf5.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_4_ase2hdf5.ipynb)

## Geometry space
Structures can be displayed via RDKit. 
Evolution of nuclear positions in the course of a trajectory can be visualized by applying PCA to pairwise distances. 
Geometric features such as bond lengths, angles and dihedrals can be calculated and e.g. used to colour the plot.
Relevant tutorials: [1_1_GS_PCA.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/1_1_GS_PCA.ipynb), [1_2_GS_cluster_dihedrals.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/1_2_GS_cluster_dihedrals.ipynb), [1_3_GS_multi_compound_biplot.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/1_3_GS_multi_compound_biplot.ipynb), [1_4_GS_diheral_PCA.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/1_4_GS_diheral_PCA.ipynb).

## Property space
Relevant tutorials: [2_1_PS_filtration.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_1_PS_filtration.ipynb), [2_2_PS_explore.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_2_PS_explore.ipynb), [2_3_PS_energy_forces.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_3_PS_energy_forces.ipynb), [2_4_PS_spectra.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_4_PS_spectra.ipynb)

## Datasheet
Plots to give an overview of an ensemble, including
- Histograms (E, F, Dipoles, ...)
- $\Delta E$ vs. $t$
- Absorption spectra at chosen (delay) times
- Population Plots
- PCA plot

Relevant tutorial: [3_1_datasheet.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/3_1_datasheet.ipynb)
<!--
## 2 Data Processing and Visualization

- [ ] Data Sheet Summary

optional:
- [ ] Histograms (E, F, Dipoles, ...)
- [ ] E vs. t
- [ ] Absorption spectra at chosen (delay) times
- [ ] Population Plots
- [ ] PCA plot (dynamic)
-->
