from shnitsel.io import read
from shnitsel.geo.geocalc import get_bats

from pytest import fixture


@fixture(
    params=[
        ('tutorials/tut_data/traj_I02.nc', -1),
    ]
)
def frames(request):
    path, charge = request.param
    db = read(path).set_charge(charge)
    return db.as_frames()


def test_get_bats(frames):
    res = get_bats(frames['atXYZ'])
    assert 'descriptor' in res.dims
