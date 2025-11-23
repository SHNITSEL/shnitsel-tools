import pytest

import xarray as xr

import shnitsel.xarray
import shnitsel.core.filtration2 as F

def load_frames(path):
    frames = xr.load_dataset(path)\
        .set_xindex(['from', 'to'])\
        .set_xindex(['trajid', 'time'])
    for attr in list(frames.attrs):
        if attr.startswith('_'):
            del frames.attrs[attr]
        
    return frames

@pytest.fixture
def frames():
    frames = load_frames('/c/st-refactor/h24-oldstyle.nc')
    frames['energy'] = frames['energy'].st.convert_energy('eV')
    frames['e_kin'] = frames['e_kin'].st.convert_energy('eV')
    return frames

def test_filtranda(frames):
    return F.energy_filtranda(frames)

@pytest.fixture
def filtranda(frames):
    return F.energy_filtranda(frames)

@pytest.fixture
def ds_filtranda(frames, filtranda):
    return frames.assign(filtranda=filtranda)

def test_cutoffs_from_filtranda(filtranda):
    return F.cutoffs_from_filtranda(filtranda)

def test_cum_mask_from_filtranda(filtranda):
    return F.cum_mask_from_filtranda(filtranda)

