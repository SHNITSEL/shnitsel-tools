<div align="center">
  <h1>shnitsel-tools</h1>
  <img src="https://raw.githubusercontent.com/SHNITSEL/shnitsel-tools/main/logo_shnitsel_tools.png" alt="SHNITSEL-TOOLS Logo" width="200px">
  <h3>Surface Hopping Nested Instances Training Set for Excited-state Learning Tools</h3>
  <br>
  <a href="https://shnitsel.github.io/">
    <img src="https://img.shields.io/badge/Website-shnitsel.github.io-yellow.svg" alt="DOI">
  </a>
</div>

--------------------

## About

`shnitsel-tools` is designed to to support the entire data lifecycle of surface hopping (SH) trajectory data upon simulation: data managment, storage, processing, visualization and interpretation. 
The tool is compatible with surface hopping data generated using the software packages [SHARC 3/4](https://sharc-md.org/), [Newton-X](https://newtonx.org/), and [PyRAI2MD](https://github.com/lopez-lab/PyRAI2MD).
The package leverages [Xarray](https://xarray.dev/) to benefit from efficient multidimensional data handling, improved metadata management, and a structure that aligns naturally with the needs of quantum chemical datasets.

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
[exploration of electronic properties](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/2_2_PS_explore.ipynb)
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

## Installation

You can create the environment with a custom path using one of the following methods:

<details open>
  <summary><strong>Option 1: Using pip</strong></summary>
  The following should work as usual:
  
  ```bash
  pip install shnitsel-tools
  ```
</details>

<details>
  <summary><strong>Option 2: Using venvs with `uv`</strong></summary>
  If you would like to work on the source code,
  we recommend installing using the `uv` tool, available at https://docs.astral.sh/uv/.  
  Run the following in the `shnitsel-tools` directory:

  ```bash
  uv venv  # create an environment under ./.venv
  . .venv/bin/activate  # activate the new environment
  uv pip install -e .  # install shnitsel in editable mode
  ```

  To install the optional development dependencies run

  ```bash
  uv pip install -e '.[dev]'
  ```
  
</details>

<details>
  <summary><strong>Option 3: Using the `--prefix` Flag</strong></summary>
  
  You can create the environment and specify the desired path by using the `conda env create` command with the `--prefix` flag:
  
  ```bash
  conda env create --prefix /home/user/anaconda3/envs/shnitsel -f shnitsel-tools.yml
  ```
</details>

<details>
  <summary><strong>Option 4: Adding the Path to the .yml File</strong></summary>
  
  Alternatively, you can manually add the desired path to the shnitsel-tools.yml file and create the environment directly:
    
  1) Open the shnitsel-tools.yml file for editing:
  
  ```bash
  vi shnitsel-tools.yml
  ```
  
  2) Add the following line to the file:
  
  
  ```
  prefix: /home/user/anaconda3/envs/shnitsel
  ```
  
  3) Create the environment with a custom path. 
  
  ```bash
  conda env create -f shnitsel-rdkit.yml
  ```
</details>

## Further Information

[![Website](https://img.shields.io/badge/Website-shnitsel.github.io-yellow.svg)](https://shnitsel.github.io/)


