# The `shnitsel` code

The `shnitsel` Python package is designed to read, process, and visualize data generated from workflows using [SHARC](https://sharc-md.org/), [Newton-X](https://newtonx.org/), and [ASE databases](https://wiki.fysik.dtu.dk/ase/ase/db/db.html).
The `shnitsel` code leverages Xarray to benefit from efficient multidimensional data handling, improved metadata management, and a structure that aligns naturally with the needs of quantum chemical datasets.

The package is organized into two main modules, reflecting the division between static and dynamic computational photochemistry methods:

- `shnitsel.static`: For processing time-independent (static) data.
- `shnitsel.dynamic`: For processing single or multiple trajectories (dynamic) data.
    
