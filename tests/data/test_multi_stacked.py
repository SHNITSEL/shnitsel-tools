from pytest import fixture

from shnitsel.data.dataset_containers.multi_layered import MultiSeriesLayered
from shnitsel.data.dataset_containers.multi_stacked import MultiSeriesStacked
from shnitsel.data.tree.node import TreeNode
from shnitsel.io import read


@fixture(
    params=[
        ('tutorials/test_data/shnitsel/traj_I02.nc', 1),
        ('tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs', 1),
        ('tutorials/test_data/newtonx/test_pyrazene_v2.6', 0),
    ]
)
def tree(request) -> TreeNode:
    path, charge = request.param
    res = read(path).set_charge(charge)
    return res


@fixture()
def stacked(tree) -> MultiSeriesStacked:
    return tree.as_stacked


def test_stacked_type(stacked):
    assert isinstance(stacked, MultiSeriesStacked)

    assert 'atrajectory' in stacked.coords
    assert 'frame' in stacked.dims
    assert 'trajectory' in stacked.dims


def test_stacked_conversion_idempotent(stacked):
    assert stacked.as_stacked is stacked


def test_stacked_conversion_to_layered(stacked):
    stacked_to_layered = stacked.as_layered
    assert isinstance(stacked_to_layered, MultiSeriesLayered)
    assert 'atrajectory' not in stacked_to_layered.dims
    assert 'atrajectory' not in stacked_to_layered.coords
    assert 'time' in stacked_to_layered.dims
    assert 'frame' not in stacked_to_layered.dims
    assert 'trajectory' in stacked_to_layered.dims
    assert stacked_to_layered.sizes['trajectory'] == stacked.sizes['trajectory']
