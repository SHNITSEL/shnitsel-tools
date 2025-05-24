import xarray as xr
import pandas as pd
import numpy as np


def read_traj(traj_path):
    # with open(os.path.join(traj_path, 'RESULTS', 'nx.log')) as f:
    #     single_traj = parse_nx_log(f)

    # nsteps = single_traj.sizes['ts']

    nstates = 3

    with open(os.path.join(traj_path, 'RESULTS', 'dyn.xyz')) as f:
        atNames, atXYZ = xyz.parse_xyz(f)

    # single_traj.attrs['atNames'] = atNames
    dims = ['ts', 'atom', 'direction']
    single_traj = single_traj.assign_coords({'atNames': ('atom', atNames)})

    if (
        not single_traj.attrs['completed']
        and atXYZ.shape[0] == single_traj.sizes['ts'] + 1
    ):
        logging.info("Geometry file contains ts after error. Truncating.")
        atXYZ = atXYZ[:-1]

    single_traj['atXYZ'] = xr.DataArray(
        atXYZ, coords={k: single_traj.coords[k] for k in dims}, dims=dims
    )

    return single_traj