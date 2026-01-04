import pytest

import xarray as xr

import shnitsel.xarray
import shnitsel.clean as F
from shnitsel.clean.common import cutoffs_from_filtranda


def load_frames(path):
    frames = (
        xr.load_dataset(path).set_xindex(['from', 'to']).set_xindex(['trajid', 'time'])
    )
    for attr in list(frames.attrs):
        if attr.startswith('_'):
            del frames.attrs[attr]

    return frames


@pytest.fixture
def frames():
    frames = load_frames('/nc/st-refactor/h24-oldstyle.nc')
    frames['energy'] = frames['energy'].st.convert_energy('eV')
    frames['e_kin'] = frames['e_kin'].st.convert_energy('eV')
    return frames


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
