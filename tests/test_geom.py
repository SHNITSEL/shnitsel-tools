from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.io import read
from shnitsel.geo.geocalc import (
    get_bats,
    get_angles,
    get_dihedrals,
    get_distances,
    get_pyramidalization,
    get_max_chromophor_BLA,
)
from shnitsel.filtering.structure_selection import StructureSelection

from pytest import fixture


@fixture(
    params=[
        ('./tutorials/test_data/shnitsel/traj_I02.nc', 1),
    ]
)
def db(request):
    path, charge = request.param
    db = read(path, expect_dtype=ShnitselDB[Trajectory]).set_charge(charge)
    return db


@fixture()
def selection(db):
    return StructureSelection.init_from_dataset(
        list(db.collect_data())[0], ['bonds', 'angles', 'dihedrals', 'pyramids']
    )


def test_get_bats(db):
    res = get_bats(db)
    assert all('descriptor' in x.dims for x in res.collect_data())


def test_get_bats_selection(db, selection):
    res = get_bats(db, selection)
    assert all('descriptor' in x.dims for x in res.collect_data())


def test_get_bats_all(db):
    res = get_bats(db, default_features=['bonds', 'angles', 'dihedrals', 'pyramids'])
    assert all('descriptor' in x.dims for x in res.collect_data())


def test_get_distances(db):
    res = get_distances(db)
    assert all('descriptor' in x.dims for x in res.collect_data())


def test_get_angles(db):
    res = get_angles(db)
    assert all('descriptor' in x.dims for x in res.collect_data())


def test_get_dihedrals(db):
    res = get_dihedrals(db)
    assert all('descriptor' in x.dims for x in res.collect_data())


def test_get_pyramidalization(db):
    res = get_pyramidalization(db)
    assert all('descriptor' in x.dims for x in res.collect_data())


def test_get_max_chromophor_BLA(db):
    res = get_max_chromophor_BLA(db)
    assert all('descriptor' in x.dims for x in res.collect_data())
