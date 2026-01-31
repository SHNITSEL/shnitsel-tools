<div align="center">
  <h1>shnitsel-tools</h1>
  <img src="https://raw.githubusercontent.com/SHNITSEL/shnitsel-tools/main/logo_shnitsel_tools.png" alt="SHNITSEL-TOOLS Logo" width="200px">
  <h3>Surface Hopping Nested Instances Training Set for Excited-state Learning Tools</h3>
  <br>
  <a href="https://shnitsel.github.io/"><img src="https://img.shields.io/badge/Website-shnitsel.github.io-yellow.svg" alt="DOI"></a>
  &nbsp;&nbsp;&nbsp;
  <a href="https://shnitsel.github.io/tools/docs/_build/index.html"><img src="https://img.shields.io/badge/Docs-shnitsel.github.io-yellow.svg" alt="DOI"></a>
</div>

--------------------

## About

`shnitsel-tools` is designed to to support the entire data lifecycle of surface hopping (SH) trajectory data upon simulation: data managment, storage, processing, visualization and interpretation. 
The tool is compatible with surface hopping data generated using the software packages [SHARC 3/4](https://sharc-md.org/), [Newton-X](https://newtonx.org/), and [PyRAI2MD](https://github.com/lopez-lab/PyRAI2MD).
The package leverages [Xarray](https://xarray.dev/) to benefit from efficient multidimensional data handling, improved metadata management, and a structure that aligns naturally with the needs of quantum chemical datasets.

## Installation

`shnitsel-tools` is normally used interactively via Jupyter Notebook on a local machine.
However, some users might find it convenient to convert trajectories to NetCDF
on-cluster, as the NetCDF file will likely download faster than the raw text files.
Either way the following should work as usual, ideally in a fresh virtual (e.g. `conda`) environment:

  ```bash
  pip install shnitsel-tools[vis]
  ```

For more detailed installation instructions, see [here](#detailed-installation-instructions)

## Usage
<!-- TODO: Make documentation more visible -->
For documentation including an API reference, please see https://shnitsel.github.io/tools/docs/_build/index.html.

`shnitsel-tools` mostly exposes data as Xarray (`xr`) objects, so familiarity with that library is beneficial.
Xarray is somewhat like Pandas for higher-dimensional data, or like Numpy with labels and other metadata.
- [Overview of data structures](https://tutorial.xarray.dev/intermediate/datastructures-intermediate.html)
- [Official quick overview](https://docs.xarray.dev/en/stable/getting-started-guide/quick-overview.html)
- [Xarray in 45 minutes](https://tutorial.xarray.dev/overview/xarray-in-45-min.html) for a more detailed introduction

### Tutorials
For a quick start, see the [tutorials](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials) directory,
which contains Jupyter Notebooks showing the workflow for parsing, writing and loading SHNITSEL databases as well as how to postprocess and visualize the respective data.

<!--
TODO: Adapt to new tutorials!
#### Collection & storage
- [parsing trajcetory and initial condition data obtained by SHARC](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_1_sharc2hdf5.ipynb)
- [parsing trajectory data produced with Newton-X](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_2_nx2hdf5.ipynb)
- [convert ASE databases](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/0_4_ase2hdf5.ipynb)
#### Management
- [exploration of electronic properties](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_2_PS_explore.ipynb)
#### Postprocessing & visualization of data
- [datasheet for trajectory data](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/3_1_datasheet.ipynb)
- [principal component analysis and trajectory classification](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/1_1_GS_PCA.ipynb)

### Workflow walkthrough
Four [notebooks](https://github.com/SHNITSEL/shnitsel-tools/tree/main/tutorials/walkthrough) demonstrate a workflow for the comparative
analysis of homologous/isoelectronic molecules, from filtration via dimensional reduction and clustering to kinetics.
-->

## Tests
Running the test-suite can help confirm that `shnitsel-tools` is correctly installed in your environment.
To run the tests, please install `shnitsel-tools` by cloning the repo as described [here](#for-tutorials-tests-or-development-conda) (conda users)
or [here](#for-tutorials-tests-or-development-uv) (uv users).
Then run the `pytest` program without arguments while inside the repository directory, with the virtual environment activated.
More information on how to invoke `pytest` is available [here](https://docs.pytest.org/en/stable/how-to/usage.html#usage).

If you are a contributor, please consider using the `tox` command to run tests across the different supported Python versions.

## Overview of repository contents
```console
$ tree -I '_*' shnitsel/
shnitsel
├── analyze
│   ├── generic.py
│   ├── hops.py
│   ├── lda.py
│   ├── pca.py
│   ├── pls.py
│   ├── populations.py
│   ├── spectra.py
│   └── stats.py
├── bridges.py
├── clean
│   ├── common.py
│   ├── dispatch_plots.py
│   ├── filter_energy.py
│   ├── filter_geo.py
│   └── filtration_class.py
├── cli
│   ├── convert_to_shnitsel.py
│   ├── generate_datasheet.py
│   └── merge_shnitsel_files.py
├── core
│   ├── feature_detection.py
│   └── typedefs.py
├── data
│   ├── atom_helpers.py
│   ├── charge_helpers.py
│   ├── dataset_containers
│   │   ├── data_series.py
│   │   ├── dataset_vis.py
│   │   ├── frames.py
│   │   ├── inter_state.py
│   │   ├── multi_layered.py
│   │   ├── multi_series.py
│   │   ├── multi_stacked.py
│   │   ├── per_state.py
│   │   ├── shared.py
│   │   ├── trajectory_collection.py
│   │   ├── trajectory.py
│   │   └── xr_conversion.py
│   ├── helpers.py
│   ├── multi_indices.py
│   ├── proxy_class.py
│   ├── shnitsel_db
│   │   └── db_function_decorator.py
│   ├── shnitsel_db_helpers.py
│   ├── state_helpers.py
│   ├── traj_combiner_methods.py
│   ├── trajectory_grouping_params.py
│   ├── trajectory_variable_description.py
│   └── tree
│       ├── child_support_functions.py
│       ├── compound.py
│       ├── data_group.py
│       ├── data_leaf.py
│       ├── datatree_level.py
│       ├── node.py
│       ├── selection.py
│       ├── support_functions.py
│       ├── tree_completion.py
│       ├── tree.py
│       ├── tree_vis.py
│       └── xr_conversion.py
├── filtering
│   ├── filter.md
│   ├── helpers.py
│   ├── state_selection.py
│   └── structure_selection.py
├── geo
│   ├── alignment.py
│   ├── analogs.py
│   ├── geocalc_
│   │   ├── algebra.py
│   │   ├── angles.py
│   │   ├── bla_chromophor.py
│   │   ├── dihedrals.py
│   │   ├── distances.py
│   │   ├── helpers.py
│   │   ├── positions.py
│   │   └── pyramids.py
│   └── geocalc.py
├── io
│   ├── ase
│   │   ├── format_reader.py
│   │   ├── parse.py
│   │   └── write.py
│   ├── format_reader_base.py
│   ├── format_registry.py
│   ├── molcas
│   │   └── molcas_opt.py
│   ├── newtonx
│   │   ├── format_reader.py
│   │   └── parse.py
│   ├── pyrai2md
│   │   ├── format_reader.py
│   │   └── parse.py
│   ├── read.py
│   ├── sharc
│   │   ├── format_reader.py
│   │   ├── parse_initial_conditions.py
│   │   ├── parse_trajectory.py
│   │   └── qm_helpers.py
│   ├── shared
│   │   ├── helpers.py
│   │   ├── messages.py
│   │   ├── trajectory_finalization.py
│   │   ├── trajectory_setup.py
│   │   └── variable_flagging.py
│   ├── shnitsel
│   │   ├── format_reader.py
│   │   ├── parse.py
│   │   └── write.py
│   ├── xr_io_compatibility.py
│   └── xyz
│       └── parse.py
├── rd.py
├── test_support
│   └── trajectory_verification.py
├── units
│   ├── conversion.py
│   ├── defaults.py
│   └── definitions.py
├── vis
│   ├── colormaps.py
│   ├── datasheet
│   │   ├── datasheet_page.py
│   │   ├── datasheet.py
│   │   └── figures
│   │       ├── common.py
│   │       ├── dip_trans_hist.py
│   │       ├── energy_bands.py
│   │       ├── hist.py
│   │       ├── nacs_hist.py
│   │       ├── per_state_hist.py
│   │       ├── socs_hist.py
│   │       ├── soc_trans_hist.py
│   │       ├── structure.py
│   │       └── time.py
│   ├── plot
│   │   ├── common.py
│   │   ├── filtration.py
│   │   ├── kde.py
│   │   ├── p3mhelpers.py
│   │   ├── pca_biplot.py
│   │   ├── polychrom.py
│   │   ├── select.py
│   │   ├── spectra3d.py
│   │   └── time.py
│   ├── static
│   │   ├── css
│   │   │   └── style.css
│   │   └── html
│   │       └── icons-svg-inline.html
│   ├── support
│   │   ├── ipython_visualization.py
│   │   ├── multi_plot.py
│   │   └── visualizeable.py
│   └── vmd
│       └── script.tcl
└── xarray.py
```

## Detailed installation instructions

### Optional dependencies
In the following, the `[vis]` suffix causes optional plotting dependencies to be
installed. If you are using `shnitsel-tools` on an HPC, you can omit it.

If you would like to contribute to `shnitsel-tools`, you may find the development
dependencies useful. These can be obtained by adding `[dev]` at the end of the
package name.

To install all optional dependencies, please add `[vis,dev]` after the package name.

### Installation using conda

#### For tutorials, tests or development (conda)
If you would like to work through the tutorials, run the tests, or modify the code, please use the following commands:

```bash
git clone 'https://github.com/SHNITSEL/shnitsel-tools.git'
cd shnitsel-tools
conda create -n shnitsel python==3.12 pip
conda activate shnitsel
pip install -e .[vis,dev]
```

#### For ordinary use
If you would just like to use the package, it is unnecessary to clone the repository.
Instead, it should suffice to run the following command with the conda environment active:

```bash
conda create -n shnitsel python==3.12 pip
conda activate shnitsel
pip install shnitsel-tools[vis]
```

### Instalation using uv
This tool is typically faster and more light-weight than `pip` and `conda`.
Unlike `conda`, it creates traditional Python virtual environments, which are stored in the folder in which the command is run and activated by sourcing
a shell-script. It can be installed by following the instructions at https://docs.astral.sh/uv/.


#### For tutorials, tests or development (uv)
If you would like to work through the tutorials or modify the code, run the tests, please use the following commands:

```bash
git clone 'https://github.com/SHNITSEL/shnitsel-tools.git'
cd shnitsel-tools
uv venv --python 3.12  # create an environment under ./.venv
source .venv/bin/activate  # activate the new environment
uv pip install -e .[vis,dev]
```

#### For ordinary use
If you would just like to use the package, it is unnecessary to clone the repository.
Instead, it should suffice to run the following commands:

```bash
uv venv --python 3.12 shnitsel  # create a directory here named ./shnitsel
source shnitsel/bin/activate  # activate the new environment
uv pip install shnitsel-tools[vis]
```

#### If you only want to use the command-line
The following will ensure the command-line programs provided are always available,
without requiring environments to be activated first.  
This may be the easiest approach for use on HPC.
It is not expected to support interactive notebook-based workflows.
```bash
uv tool install shnitsel-tools
```

<!--
### For developers
  
  We recommend installation using the `uv` tool, available at https://docs.astral.sh/uv/.
  Please clone this repo and run the following in the `shnitsel-tools` directory:

  ```bash
  git clone 'https://github.com/SHNITSEL/shnitsel-tools.git'
  cd shnitsel-tools
  uv venv  # create an environment under ./.venv
  source .venv/bin/activate  # activate the new environment
  uv pip install -e .[dev]  # install shnitsel in editable mode
  ```

  In the above, the option `-e` installs in editable mode, meaning that Python will see changes you make
  to the source, while `[dev]` installs the optional development dependencies.  

  If you would like to contribute your changes,
  please [fork](https://github.com/SHNITSEL/shnitsel-tools/fork) this repo,
  and make a pull request.

-->

## Further Information

[![Website](https://img.shields.io/badge/Website-shnitsel.github.io-yellow.svg)](https://shnitsel.github.io/)


