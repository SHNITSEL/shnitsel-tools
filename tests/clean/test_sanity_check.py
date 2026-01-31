import pytest

import shnitsel.clean as F
from shnitsel.clean.common import true_upto, _filter_mask_from_filtranda
from shnitsel.clean.filter_energy import filter_by_energy, calculate_energy_filtranda
from shnitsel.clean.filter_geo import filter_by_length, calculate_bond_length_filtranda
from shnitsel.data.dataset_containers.frames import Frames
from shnitsel.io import read


@pytest.fixture(
    params=[
        ('tutorials/test_data/shnitsel/traj_I02.nc', -1),
    ]
)
def frames(request) -> Frames:
    path, charge = request.param
    db = read(path)
    res = db.set_charge(charge)
    return res.as_stacked


def test_filtranda(frames):
    calculate_energy_filtranda(frames)


@pytest.fixture
def energy_filtranda(frames):
    return calculate_energy_filtranda(frames)


@pytest.fixture
def length_filtranda(frames):
    return calculate_energy_filtranda(frames)


@pytest.fixture
def ds_filtranda(frames, energy_filtranda):
    return frames.assign(filtranda=energy_filtranda)


def test_cutoffs_from_filtranda(energy_filtranda):
    true_upto(_filter_mask_from_filtranda(energy_filtranda), 'criterion')


def test_mask_from_filtranda(energy_filtranda):
    _filter_mask_from_filtranda(energy_filtranda)


def test_filter_energy(frames):
    filter_by_energy(frames)


def test_filter_length(frames):
    filter_by_length(frames)
