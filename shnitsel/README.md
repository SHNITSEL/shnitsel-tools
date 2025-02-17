# The `shnitsel` code

The `shnitsel` Python package is designed to read, process, and visualize data generated from workflows using [SHARC](https://sharc-md.org/), [Newton-X](https://newtonx.org/), and [ASE databases](https://wiki.fysik.dtu.dk/ase/ase/db/db.html).
The `shnitsel` code leverages Xarray to benefit from efficient multidimensional data handling, improved metadata management, and a structure that aligns naturally with the needs of quantum chemical datasets.

The package is organized into two main modules, reflecting the division between static and dynamic computational photochemistry methods:

- `shnitsel.static`: For processing time-independent (static) data.
- `shnitsel.dynamic`: For processing single or multiple trajectories (dynamic) data.

```bash
shnitsel
├── dynamic
│   ├── datasheet
│   │   ├── calc
│   │   │   ├── __init__.py
│   │   │   └── spectra.py
│   │   ├── __init__.py
│   │   ├── oop.py
│   │   ├── plot
│   │   │   ├── colormaps.py
│   │   │   ├── common.py
│   │   │   ├── dip_trans_hist.py
│   │   │   ├── hist.py
│   │   │   ├── __init__.py
│   │   │   ├── nacs_hist.py
│   │   │   ├── per_state_hist.py
│   │   │   ├── structure.py
│   │   │   └── time.py
│   ├── filter_unphysical.py
│   ├── indexes.py
│   ├── __init__.py
│   ├── parse
│   │   ├── common.py
│   │   ├── __init__.py
│   │   ├── nx.py
│   │   ├── sharc_icond.py
│   │   ├── sharc_traj.py
│   │   └── xyz.py
│   ├── pca_biplot.py
│   ├── plot
│   │   ├── dihedral_kde.py
│   │   ├── __init__.py
│   │   ├── p3mhelpers.py
│   ├── plotting.py
│   ├── postprocess.py
│   └── xrhelpers.py
├── __init__.py
└── static
    ├── __init__.py
    ├── pca.py
    ├── plotting.py
    └── utils.py
```
