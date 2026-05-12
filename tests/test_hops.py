from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.tree.tree import ShnitselDB
from shnitsel.io import read
from shnitsel.analyze.hops import (
    assign_hop_time,
    hops_mask_from_active_state,
    hops,
    filter_data_at_hops,
    focus_hops,
)
from shnitsel.filtering.state_selection import StateSelection

from pytest import fixture


@fixture(
    params=[
        ('./tutorials/test_data/shnitsel/traj_I02.nc', 1),
    ]
)
def db(request) -> ShnitselDB[Trajectory]:
    path, charge = request.param
    db = read(path, expect_dtype=ShnitselDB[Trajectory]).set_charge(charge)
    return db


@fixture()
def full_state_selection(db: ShnitselDB[Trajectory]):
    return StateSelection.init_from_dataset(list(db.collect_data())[0])


labels = ["db", "stacked", "layered"]
static_selections = ["0, 1, 1->0", "S, S<>S"]


@fixture()
def frames(db: ShnitselDB[Trajectory]):
    return db.as_stacked


@fixture()
def layered(db: ShnitselDB[Trajectory]):
    return db.as_layered


# Testing whether hops_mask_from_active_state() runs without crashing
def test_hops_mask_from_active_state(
    subtests, db, frames, layered, full_state_selection
):
    for i, (name, object) in enumerate(zip(labels, [db, frames, layered])):
        with subtests.test(msg=f"Input format: {name}", i=i):
            hops_mask_from_active_state(object)
            hops_mask_from_active_state(object, full_state_selection)
            for static_sel in static_selections:
                hops_mask_from_active_state(object, static_sel)


# Testing whether assign_hop_time() runs without crashing
def test_assign_hop_time(subtests, db, frames, layered, full_state_selection):
    for i, (name, object) in enumerate(zip(labels, [db, frames, layered])):
        with subtests.test(msg=f"Input format: {name}", i=i):
            assign_hop_time(object)
            assign_hop_time(object, full_state_selection)
            for static_sel in static_selections:
                assign_hop_time(object, static_sel)


# Testing whether hops() runs without crashing
def test_hops(subtests, db, frames, layered, full_state_selection):
    for i, (name, object) in enumerate(zip(labels, [db, frames, layered])):
        with subtests.test(msg=f"Input format: {name}", i=i):
            hops(object)
            hops(object, full_state_selection)
            for static_sel in static_selections:
                hops(object, static_sel)


# Testing whether filter_data_at_hops runs without crashing
def test_filter_data_at_hops(subtests, db, frames, layered, full_state_selection):
    for i, (name, object) in enumerate(zip(labels, [db, frames, layered])):
        with subtests.test(msg=f"Input format: {name}", i=i):
            filter_data_at_hops(object)
            filter_data_at_hops(object, full_state_selection)
            for static_sel in static_selections:
                filter_data_at_hops(object, static_sel)


# Testing whether focus_hops runs without crashing
def test_focus_hops(subtests, db, frames, layered, full_state_selection):
    for i, (name, object) in enumerate(zip(labels, [db, frames, layered])):
        with subtests.test(msg=f"Input format: {name}", i=i):
            focus_hops(object)
            focus_hops(object, full_state_selection)
            for static_sel in static_selections:
                focus_hops(object, static_sel)
