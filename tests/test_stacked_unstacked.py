import pytest
from xarray.testing import assert_equal, assert_allclose


import shnitsel as st
import shnitsel.xarray

PER_TRAJ_COORDS = [
    'trajectory_input_path',
    't_max',
    'est_level',
    'theory_basis_set',
    'completed',
    'max_ts',
    'has_forces',
    'input_format',
]


@pytest.fixture
def stacked():
    path = './tutorials/test_data/shnitsel/fixtures/butene_dynamic/data_new.nc'

    wrapped = st.io.read(path, input_state_names=['S0', 'S1', 'S2'])
    frames = wrapped.set_charge(0)
    # frames.attrs['charge'] = 0
    # frames.attrs['mol'] = frames.st.default_mol()
    # frames.atXYZ.attrs['mol'] = frames.attrs['mol']
    return frames


@pytest.fixture
def unstacked(stacked):
    return stacked.st.unstack_trajs()


@pytest.mark.parametrize(
    "func,data_var,kws",
    [
        (st.analyze.generic.norm, None, {}),
        (st.analyze.generic.center, 'atXYZ', dict(dim='atom')),
        (st.analyze.generic.subtract_combinations, 'atXYZ', dict(dim='atom')),
        # skipping replace_total -- hard to test and metadata-agnostic anyway
        (st.analyze.generic.relativize, 'energy', {}),
        # (st.analyze.generic.pwdists, None, {}), # FIXME: ds.dataset
        # (st.analyze.hops.hops, None, {}),
        # (st.analyze.hops.focus_hops, None, {}),
        # (st.analyze.hops.assign_hop_time, None, {}),
        # (st.analyze.populations.calc_classical_populations, None, {}),
        # (st.analyze.spectra.get_spectra, None, {}),
        # (st.analyze.stats.time_grouped_confidence_interval, 'energy', {}),
    ],
)
def test_stacked_equal_unstacked(stacked, unstacked, func, data_var, kws):
    from shnitsel.data.multi_indices import unstack_trajs

    # args = args or []
    kws = kws or {}
    if data_var is not None:
        input_stacked = stacked[data_var]
        input_unstacked = unstacked[data_var].drop_vars(
            PER_TRAJ_COORDS, errors='ignore'
        )
    else:
        input_stacked = stacked
        input_unstacked = unstacked
    res_stacked = func(input_stacked, **kws)
    res_unstacked = func(input_unstacked, **kws)
    assert_equal(unstack_trajs(res_stacked), res_unstacked)
