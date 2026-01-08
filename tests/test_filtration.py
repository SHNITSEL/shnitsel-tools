import pytest

import shnitsel.clean as F
from shnitsel.clean.common import cutoffs_from_filtranda
from shnitsel.data.tree import tree_to_frames
from shnitsel.io import read


@pytest.fixture(
    params=[
        ('tutorials/tut_data/traj_I02.nc', -1),
    ]
)
def frames(request):
    path, charge = request.param
    db = read(path)
    res = tree_to_frames(db)
    res['atXYZ'].attrs['charge'] = charge
    res.attrs['charge'] = charge
    return res


def test_filtranda(frames):
    F.energy_filtranda(frames)


@pytest.fixture
def filtranda(frames):
    return F.energy_filtranda(frames)


@pytest.fixture
def ds_filtranda(frames, filtranda):
    return frames.assign(filtranda=filtranda)


def test_cutoffs_from_filtranda(filtranda):
    cutoffs_from_filtranda(filtranda)


def test_cum_mask_from_filtranda(filtranda):
    F.cum_mask_from_filtranda(filtranda)
