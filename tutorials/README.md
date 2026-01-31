# Tutorials

These tutorials can be downloaded from this repo, or viewed in HTML form [here](https://shnitsel.github.io/tools/docs/_build/tutorials.html).

## Full tutorial 

A good place to start, this tutorial showcases a wide array of `shnitsel-tools` functionality.

Relevant tutorial: [0_full_tutorial.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_full_tutorial.ipynb)

## General IO with `shnitsel-tools`

`shnitsel-tools` can read the following:

- [SHARC](https://sharc-md.org/) output files: initial conditions (`ICOND_*` folders) or trajectories (`TRAJ_*` folders)
- [NewtonX](https://newtonx.org/) output files: trajectories 
- [PyRAI2MD](https://github.com/lopez-lab/PyRAI2MD) output files: trajectories 
- [ASE](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) databases

Relevant tutorial: [2_2_general_io.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_2_general_io.ipynb)

## Selection

Learn how to flexibly and expressively select geometrical features and electronic states for focussed analysis.

Relevant tutorials: [2_4a_filtering_structure.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_4a_filtering_structure.ipynb), [2_4b_filtering_states.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_4b_filtering_states.ipynb)

## Datasheet
`shnitsel-tools` provides preset plots to give an overview of ensembles, including
- Histograms (Energy, Forces, Dipoles, ...)
- $\Delta E$ vs. $t$
- Absorption spectra at chosen (delay) times
- Population Plots
- PCA of geometrical features

Relevant tutorial: [2_6_b_datasheet.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_6_b_datasheet.ipynb)

## Walkthrough: across compounds
`shnitsel-tools` can accommodate ensembles of multiple compounds in the same data structure.
In this case, we compare retinal and methyleneimmonium.

Relevant tutorial: [3_across_compounds.ipynb](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/3_across_compounds.ipynb)

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
