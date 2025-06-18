import py3Dmol
import xarray as xr
from math import sqrt, ceil
from .. import postprocess as P


def frame3D(atXYZ_frame: str | xr.DataArray):
    if isinstance(atXYZ_frame, xr.DataArray):
        atXYZ_frame = P.to_xyz(atXYZ_frame)
    view = py3Dmol.view()
    view.addModel(atXYZ_frame)

    view.setStyle({'stick': {'showNonBonded': True}})
    view.zoomTo()
    view.show()
    return view


def frames3Dgrid(atXYZ: xr.DataArray):
    n = ceil(sqrt(atXYZ.sizes['frame']))
    view = py3Dmol.view(viewergrid=(n, n), width=1000, height=800, linked=True)

    for i, (label, frameXYZ) in enumerate(atXYZ.groupby('frame')):
        data = frameXYZ.pipe(P.to_xyz)
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
    view.show()
    return view


def traj3D(traj: str | xr.DataArray):
    if isinstance(traj, xr.DataArray):
        traj = P.traj_to_xyz(traj)
    view = py3Dmol.view()
    view.addModelsAsFrames(traj)

    view.setStyle({'stick': {'showNonBonded': True}})
    view.zoomTo()
    view.animate({'loop': "forward"})
    view.show()
    return view


def trajs3Dgrid(atXYZ: xr.DataArray, trajids: list[int | str]):
    n = ceil(sqrt(len(trajids)))
    view = py3Dmol.view(viewergrid=(n, n), width=1000, height=800, linked=True)

    for i, trajid in enumerate(trajids):
        data = atXYZ.sel(trajid=trajid).pipe(P.traj_to_xyz)
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
    view.animate({'loop': "forward"})
    view.show()
    return view