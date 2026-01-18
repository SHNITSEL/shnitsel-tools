from pytest import fixture

from shnitsel.io import read
# from shnitsel.data.tree import tree_to_frames


@fixture(
    params=[
        'tutorials/test_data/shnitsel/traj_I02.nc',
        'tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs',
        'tutorials/test_data/newtonx/test_I01_v2.6',
    ]
)
def tree(request):
    res = read(request.param)
    return res


def test_tree_to_stacked(tree):
    frames = tree.as_stacked
    assert ('atrajectory' in frames.dims)
    assert ('trajectory' in frames.dims)


def test_tree_to_layered(tree):
    layered = tree.as_layered
    assert ('atrajectory' not in layered.dims)
    assert ('trajectory' in layered.dims)
