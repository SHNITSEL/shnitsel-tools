from pytest import fixture
import pytest

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
def layered(tree) -> MultiSeriesLayered:
    return tree.as_layered


def test_layered_type(layered):
    # print(layered)
    assert isinstance(layered, MultiSeriesLayered)
    assert 'atrajectory' not in layered.dims
    assert 'atrajectory' not in layered.coords
    assert 'time' in layered.dims
    assert 'frame' not in layered.dims
    assert 'trajectory' in layered.dims


def test_layered_conversion_idempotent(layered):
    assert layered.as_layered is layered

@pytest.mark.xfail
def test_layered_conversion_to_stacked(layered):
    layered_to_stacked = layered.as_stacked
    assert isinstance(layered_to_stacked, MultiSeriesStacked)
    assert 'atrajectory' in layered_to_stacked.coords
    assert 'frame' in layered_to_stacked.dims
    assert 'trajectory' in layered_to_stacked.dims
    assert layered_to_stacked.sizes['trajectory'] == layered.sizes['trajectory']
