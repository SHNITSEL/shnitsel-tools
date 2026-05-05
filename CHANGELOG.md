# Changelog

Shnitsel tools is a project under active development. 
to keep track of the updates happening in various version, we list the 
changes here


## shnitsel-tools v2026.4.1

Bug fixes compared to previous version and addition of new features:

### Breaking Changes

* To keep the structure of shnitsel clearer, we moved the PCA, LDA, PLS, TICA and all dimensionality reduction related methods into the `shnitsel.analyze.dimred` submodule whereas they were previously located directly in `shnitsel.analyze`


### Features

* State-based and label-based coloring suppport in dimensionality reduction plots via `biplot_kde()`
* Tentative support for Linear Discriminant analysis (LDA) on datasets
* Initial support for Time-lagged independent component analysis (TICA) to analyze time-correlation between features. Not yet fully integrated into the default dimensionality reduction framework in `shnitsel.analyze`
* Added accessors for creating timeplots and performing pca analysis to `shnitsel.xarray` support.
* Improved VMD visualization support
* Add generic debug information and version information support with `shnitsel.__version__` and `shnitsel.show_debug_info()` to help with error debugging support
* Implement `focus_hops()` and `hops()` functionality to obtain more detailed analysis tools related to hop-relative statistics within trajectories
* Add reading of `state_coeffs_diag` from sharc trajectories, i.e. the coefficients associated with the diagonal basis states.
* Add reading of `prob_hop_diag` (hopping probability in diagonal basis) and `u_matrix` (Basis transform matrix from diagonal to MCH state coefficients) from SHARC data.

### Fixes

* Icons and CSS files for tree visualization were not included in the initial release package version on pip
* Colorbar label was not set in `biplot_kde()` for `time` coloring property
* Imports of SHARC formats did not correctly read and set kinetic energy in imported dataset
* Fix broken imports in `shnitsel.vis.support` module
* Fix broken logic in `hop_mask_from_active_state()` function
* Fix construction of `default_mol()` if mol information is stored in `atXYZ` variable of a dataset for legacy support.
* Fix issues arising from too strict constraints on different ShnitselDataset container types leading to issues when filtering data.
* Fix state selection and other state-related operations failing if the state dimension was reduced to a scalar coordinate during filtering.
* Fix NewtonX import failing due to bad type specifier for atom type names
* Fix `sanity_check()` failing on stacked trajectories due to bad invocation of `.cumprod()`
* Fix plots failing due to change in naming conventions.
* Fix issues with constructing default molecular structures from layered trajectories
* Fix filtering for physical sanity not working on layered trajectories
* Fix multi-index methods failing on layered trajectories. Will now simply not modify the trajectory if they cannot be applied and print out a warning.


### Known issues

* Selection of `pyramids` via `SMARTS` patterns where `StructureSelection` syntax is allowed is currently not working as intended
* For some collections of trajectories, `as_stacked()` or other combiner methods may break if there are empty trajectories, e.g. after `sanity_check()`


## shnitsel-tools v2026.1.2

Initial version of the published toolkit.


### Features
* Import of `SHARC`, `PyrAI2md` and `NewtonX` output formats
* Export to ASE databses and Netcdf 4.0/HDF5 files 
* General `Datasheet` support for generating a quick overview for datasets
* Tree manipulation support of imported data
* Extensive Metadata retention in tree-based xarray-compatible formats
* Filtering of data for physically reasonable properties via `sanity_check()`
* Filtering for hopping points SH trajectories
* Selection of specific states and state transitions for analysis with `StateSelection` framework
* Selection of specific geometrical features (positions, distances, angles, dihedrals/torsions, pyramidalization) via the `StructureSelection` framework
* Postprocessing analysis for state population, oscillator strength, dimensionality reduction, geometrical feature analysis, etc.
* Visualization support of dimensionality reduction, main contributing features, time series data, 