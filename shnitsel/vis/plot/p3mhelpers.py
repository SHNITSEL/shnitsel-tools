from math import sqrt, ceil

import numpy as np
import xarray as xr

from shnitsel.bridges import to_xyz, traj_to_xyz
from shnitsel._contracts import needs
from shnitsel._state import HAS_IPYTHON

if HAS_IPYTHON:
    import py3Dmol

@needs(coords_or_vars={'atNames'}, dims={'atom', 'direction'}, not_dims={'frame'})
def frame3D(atXYZ_frame: str | xr.DataArray):
    """Display a single geometry using py3Dmol

    Parameters
    ----------
    atXYZ_frame
        The geometry to display

    Returns
    -------
        The View object created
    """
    if isinstance(atXYZ_frame, xr.DataArray):
        atXYZ_frame = to_xyz(atXYZ_frame)
    view = py3Dmol.view()
    view.addModel(atXYZ_frame)

    view.setStyle({'stick': {'showNonBonded': True}})
    view.zoomTo()
    return view


@needs(groupable={'frame'}, dims={'atom', 'direction'}, coords_or_vars={'atNames'})
def frames3Dgrid(atXYZ: xr.DataArray):
    """Display several geometries stacked along a ``frame`` dimension using py3Dmol

    Parameters
    ----------
    atXYZ
        The geometry to display

    Returns
    -------
        The View object created
    """
    n = ceil(sqrt(atXYZ.sizes['frame']))
    view = py3Dmol.view(viewergrid=(n, n), width=1000, height=800, linked=True)

    for i, (label, frameXYZ) in enumerate(atXYZ.groupby('frame')):
        if 'frame' in frameXYZ.dims:
            frameXYZ = frameXYZ.squeeze('frame')
        data = frameXYZ.pipe(to_xyz)
        viewer = (i // n, i % n)
        view.addModel(data, viewer=viewer)
        view.addLabel(
            label,
            {
                'useScreen': True,
                'screenOffset': {"x": 25, "y": -50},
            },
            viewer=viewer,
        )

    view.setStyle({'stick': {'showNonBonded': True}})
    view.zoomTo()
    return view


@needs(groupable={'time'}, dims={'atom', 'direction'}, coords_or_vars={'atNames'})
def traj3D(traj: str | xr.DataArray):
    """Display a trajectory using py3Dmol

    Parameters
    ----------
    traj
        The trajectory geometries to display

    Returns
    -------
        The View object created
    """
    if isinstance(traj, xr.DataArray):
        traj = traj_to_xyz(traj)
    view = py3Dmol.view()
    view.addModelsAsFrames(traj)

    view.setStyle({'stick': {'showNonBonded': True}})
    view.zoomTo()
    view.animate({'loop': "forward"})
    return view


@needs(
    groupable={'time'},
    dims={'atom', 'direction'},
    coords={'atrajectory'},
    coords_or_vars={'atNames'},
)
def trajs3Dgrid(
    atXYZ: xr.DataArray, trajids: list[int | str] | None = None, loop: str = 'forward'
):
    """Display a trajectory using py3Dmol

    Parameters
    ----------
    traj
        The trajectory geometries to display
    trajids
        If given, only show these trajectories
    loop
        Passed to ``py3Dmol``'s ``view.animate``. Accepted values include
        'forward', 'backward', 'backAndForth'.

    Returns
    -------
        The View object created
    """
    if trajids is None:
        trajids = np.unique(atXYZ.coords['atrajectory'].values)

    n = ceil(sqrt(len(trajids)))
    view = py3Dmol.view(viewergrid=(n, n), width=1000, height=800, linked=True)

    for i, trajid in enumerate(trajids):
        data = atXYZ.sel(atrajectory=trajid).pipe(traj_to_xyz)
        viewer = (i // n, i % n)
        view.addModelsAsFrames(data, viewer=viewer)
        view.addLabel(
            trajid,
            {
                'useScreen': True,
                'screenOffset': {"x": 25, "y": -50},
            },
            viewer=viewer,
        )

    view.setStyle({'stick': {'showNonBonded': True}})
    view.zoomTo()
    view.animate({'loop': loop})
    return view