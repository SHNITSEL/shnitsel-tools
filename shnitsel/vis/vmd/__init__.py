import os
import subprocess
import tempfile

import numpy as np

from shnitsel.bridges import traj_to_xyz

_tcl_script_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'script.tcl'
)


def traj_vmd(atXYZ, groupby='atrajectory'):
    """Open geometries in the VMD viewer, if installed

    Parameters
    ----------
    atXYZ
        The geometries to transmit
    groupby, optional
        A set of frames will be grouped into a VMD molecule if
        they have the same value in this coordinate, by default this is 'trajid'
        so that each trajectory
    """
    # See git history of this file for an attempt to communicate
    # settings to VMD via a generated file
    with tempfile.TemporaryDirectory() as d:
        paths = []
        # TODO: Why not use `.groupby` and then `.squeeze`?
        trajids = np.unique(atXYZ.coords[groupby].values)
        for trajid in trajids:
            traj = atXYZ.loc[{groupby: trajid}]
            path = os.path.join(d, f"{trajid}.xyz")
            with open(path, 'w') as f:
                print(traj_to_xyz(traj), file=f)
            paths.append(path)
        subprocess.call(['vmd', '-e', _tcl_script_path, '-m'] + paths)