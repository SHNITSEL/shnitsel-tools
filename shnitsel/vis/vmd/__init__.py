import os
import subprocess
import tempfile

import numpy as np

from shnitsel.bridges import traj_to_xyz

_tcl_script_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'script.tcl'
)


def traj_vmd(atXYZ, groupby='trajid', scale=0.5):
    with tempfile.TemporaryDirectory() as d:
        # settings_path = os.path.join(d, "settings.tcl")
        # with open(settings_path, 'w') as f:
        #     print("running settings file", file=f)
        #     print(f"set shnitsel_scale {scale}", file=f)
        #     print("echo shnitsel_scale=$shnitsel_scale", file=f)

        paths = []
        trajids = np.unique(atXYZ.coords[groupby].values)
        for trajid in trajids:
            traj = atXYZ.sel(trajid=trajid)
            path = os.path.join(d, f"{trajid}.xyz")
            with open(path, 'w') as f:
                print(traj_to_xyz(traj), file=f)
            paths.append(path)
        subprocess.call(
            # ['vmd', '-e', settings_path, '-e', _tcl_script_path, '-m'] + paths
            ['vmd', '-e', _tcl_script_path, '-m'] + paths
            # ['vmd', '-e', settings_path, '-m'] + paths
        )