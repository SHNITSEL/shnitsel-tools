# The `shnitsel` code

The `shnitsel` Python package is designed to read, process, and visualize data generated from workflows using SHARC, Newton-X, and ASE databases.
The `shnitsel` code leverages Xarray to benefit from efficient multidimensional data handling, improved metadata management, and a structure that aligns naturally with the needs of quantum chemical datasets.

The package is organized into two main modules, reflecting the division between static and dynamic computational photochemistry methods:

- `shnitsel.static`: For processing time-independent (static) data.
- `shnitsel.dynamic`: For processing single or multiple trajectories (dynamic) data.
    
