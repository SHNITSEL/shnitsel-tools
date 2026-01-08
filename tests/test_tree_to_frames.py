from pytest import fixture

from shnitsel.io import read
from shnitsel.data.tree import tree_to_frames


@fixture(
    params=[
        'tutorials/tut_data/traj_I02.nc',
        'tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs',
        'tutorials/test_data/newtonx/test_I01_v2.6',
    ]
)
def tree(request):
    res = read(request.param)
    return res


def test_tree_to_frames(tree):
    frames = tree_to_frames(tree)
