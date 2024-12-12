<div align="center">
  <h1>SHNITSEL</h1>
  <img src="https://github.com/SHNITSEL/shnitsel-tools/blob/main/shnitsel_logo.png" alt="SHNITSEL Logo" width="200px">
  <h3>Surface Hopping Newly Invented Training Set for Excited-state Learning</h3>
</div>

--------------------

## About

Code and tools for parsing, analyzing, and visualizing data from static and dynamic calculations in SHARC and Newton-X format.

## Installation

You can create the environment with a custom path using one of the following methods:

### Option 1: Using the `--prefix` Flag

You can create the environment and specify the desired path by using the `conda env create` command with the `--prefix` flag:

```bash
conda env create --prefix /home/user/anaconda3/envs/shnitsel -f shnitsel-tools.yml
```

### Option 2: Adding the Path to the .yml File

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

## Citation

