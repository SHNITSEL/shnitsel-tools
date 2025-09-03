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
Either way the following should work as usual, ideally in a virtual (e.g. `conda`) environment:
  
  ```bash
  pip install shnitsel-tools
  ```

## Usage

The package is organized into top-level functions for reading data,
accessor methods available on `xr.Dataset` and `xr.DataArray` objects, plotting routines found in the `shnitsel.plot` namespace,
and functions taking an RDKit `Mol` object as their principal argument under `shnitsel.rd`.
The adventurous may find something useful under `shnitsel.core`, though this should be considered internal and therefore subject to change.

### Tutorials
For a quick start, see the [tutorials](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials) directory,
which contains Jupyter Notebooks showing the workflow for parsing, writing and loading SHNITSEL databases as well as how to postprocess and visualize the respective data.

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

## Tree

```bash
shnitsel
├── core
│   ├── ase.py
│   ├── datasheet
│   │   ├── colormaps.py
│   │   ├── common.py
│   │   ├── dip_trans_hist.py
│   │   ├── hist.py
│   │   ├── __init__.py
│   │   ├── nacs_hist.py
│   │   ├── oop.py
│   │   ├── per_state_hist.py
│   │   ├── structure.py
│   │   └── time.py
│   ├── filter_unphysical.py
│   ├── filtre.py
│   ├── indexes.py
│   ├── __init__.py
│   ├── parse
│   │   ├── common.py
│   │   ├── __init__.py
│   │   ├── nx.py
│   │   ├── pyrai2md.py
│   │   ├── sharc_icond.py
│   │   ├── sharc_traj.py
│   │   └── xyz.py
│   ├── plot
│   │   ├── __init__.py
│   │   ├── kde.py
│   │   ├── p3mhelpers.py
│   │   ├── pca_biplot.py
│   │   ├── polychrom.py
│   │   ├── select.py
│   │   └── spectra3d.py
│   ├── postprocess.py
│   ├── spectra.py
│   └── xrhelpers.py
├── __init__.py
├── plot
│   └── __init__.py
├── rd.py
└── xarray.py
```

## Development
  
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

## Further Information

[![Website](https://img.shields.io/badge/Website-shnitsel.github.io-yellow.svg)](https://shnitsel.github.io/)


