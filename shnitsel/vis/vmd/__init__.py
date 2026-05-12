import logging
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any

import numpy as np

from shnitsel.bridges import traj_to_xyz
from shnitsel.data.dataset_containers.shared import ShnitselDataset
from shnitsel.data.tree.node import TreeNode

import xarray as xr

_tcl_script_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'script.tcl'
)


def traj_vmd(
    atXYZ: xr.Dataset
    | xr.DataArray
    | ShnitselDataset
    | TreeNode[Any, xr.Dataset | xr.DataArray | ShnitselDataset],
    groupby='atrajectory',
    path: os.PathLike | None = None,
):
    """Open geometries in the VMD viewer, if installed

    Parameters
    ----------
    atXYZ : xr.Dataset| xr.DataArray | ShnitselDataset | TreeNode[Any, xr.Dataset| xr.DataArray | ShnitselDataset]
        The geometries to transmit. Or datasets holding the geometries
    groupby : Hashable, default='atrajectory'
        A set of frames will be grouped into a VMD molecule if
        they have the same value in this coordinate, by default this is 'trajid'
        so that each trajectory has its own entity in vmd
    path : os.PathLike, optional
        Optionally a path to a directory that the vmd files should be written to.

    """
    # See git history of this file for an attempt to communicate
    # settings to VMD via a generated file

    positional_data = atXYZ
    if not isinstance(atXYZ, xr.DataArray):
        if isinstance(atXYZ, TreeNode):
            atXYZ = atXYZ.as_stacked

        positional_data = atXYZ.atXYZ

    assert isinstance(positional_data, xr.DataArray), (
        "Could not extract positional information from provided `atXYZ` parameter."
    )

    # If no path is provided, we create a temporary directory
    if path is None:
        tmpdir = tempfile.TemporaryDirectory()
        directory = Path(tmpdir.name)
    else:
        tmpdir = None
        directory = Path(path)
        directory.mkdir(exist_ok=True)

    paths = []
    # TODO: Why not use `.groupby` and then `.squeeze`?
    trajids = np.unique(positional_data.coords[groupby].values)
    for trajid in trajids:
        _tmp_traj = positional_data.loc[{groupby: trajid}]
        _path = os.path.join(directory, f"{trajid}.xyz")

        with open(_path, 'w') as f:
            print(traj_to_xyz(_tmp_traj), file=f)
        paths.append(_path)


    if path is not None:
        # Copy over tcl script
        with open(_tcl_script_path) as f_in:
            with open(directory / "script.tcl", "w") as f_out:
                print(f_in.readlines(), file=f_out)
        print(f"VMD files written to `{path}`")
        cmd_paths = [f"\"{_path}\"" for _path in paths]
        print(
            f"Use `vmd -e {directory / 'script.tcl'} -m {' '.join(cmd_paths)}` to show vmd results."
        )

    try:
        print("Sending structure to vmd...")
        res = subprocess.call(['vmd', '-e', _tcl_script_path, '-m'] + paths)
        if res:
            logging.error("Error trying to run `vmd`")
    except Exception as e:
        logging.error(f"Failed to run `vmd`: {e}")

    # Trigger cleanup
    if tmpdir is not None:
        tmpdir.cleanup()
