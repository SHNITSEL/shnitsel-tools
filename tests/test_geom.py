from shnitsel.io import read
from shnitsel.data.tree import tree_to_frames
from shnitsel.geo.geocalc import get_bats

from pytest import fixture


@fixture(
    params=[
        ('tutorials/tut_data/traj_I02.nc', -1),
    ]
)
def frames(request):
    path, charge = request.param
    db = read(path)
    res = tree_to_frames(db)
    res['atXYZ'].attrs['charge'] = charge
    return res


def test_get_bats(frames):
    res = get_bats(frames['atXYZ'])
    assert 'descriptor' in res.dims
