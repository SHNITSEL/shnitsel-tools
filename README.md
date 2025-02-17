<div align="center">
  <h1>SHNITSEL</h1>
  <img src="https://github.com/SHNITSEL/shnitsel-tools/blob/main/shnitsel_logo.png" alt="SHNITSEL Logo" width="200px">
  <h3>Surface Hopping Newly Invented Training Set for Excited-state Learning</h3>
</div>

--------------------

## About

Code and tools for parsing, analyzing, and visualizing data from static and dynamic calculations in SHARC and Newton-X format.


> [!TIP]
> ### Tutorials
> #### Generation of SHNITSEL databases
> - [parsing trajcetory and initial condition data obtained by SHARC](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/01_sharc2hdf5.ipynb)
> - [parsing trajectory data produced with Newton-X](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/02_nx2hdf5.ipynb)
> - [convert ASE databases](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/03_ase2hdf5.ipynb)
>
> #### Postprocessing and visualization of data
> - static data
> - [processing of dynamic data](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/10_usage.ipynb)
> - [datasheet for trajectory data](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/11_datasheet.ipynb)
> - [principal component analysis and trajectory classification](https://github.com/SHNITSEL/shnitsel-tools/blob/main/tutorials/12_butene_pca.ipynb)
>

## Installation

You can create the environment with a custom path using one of the following methods:

<details open>
  <summary><strong>Option 1: Using `uv`</strong></summary>
  I recommend using the `uv` tool, available at https://docs.astral.sh/uv/.  
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

<details open>
  <summary><strong>Option 2: Using the `--prefix` Flag</strong></summary>
  
  You can create the environment and specify the desired path by using the `conda env create` command with the `--prefix` flag:
  
  ```bash
  conda env create --prefix /home/user/anaconda3/envs/shnitsel -f shnitsel-tools.yml
  ```
</details>

<details>
  <summary><strong>Option 3: Adding the Path to the .yml File</strong></summary>
  
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

## Citation

