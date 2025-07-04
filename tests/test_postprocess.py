import pytest
from hypothesis import assume, given
import hypothesis.strategies as st
from xarray.testing import assert_equal
import xarray.testing.strategies as xrst

import numpy as np
import xarray as xr

import shnitsel
from shnitsel import postprocess
from shnitsel.core.postprocess import (
    norm,
    subtract_combinations,
    pca,
    _sudi_groupby,
    sudi,
)


@pytest.fixture
def traj_butene():
    frames = shnitsel.read_trajs('tutorials/test_data/sharc/traj_butene', kind='sharc')
    return postprocess.ts_to_time(frames)


@given(
    xrst.variables(
        dims=st.just({'test1': 2, 'direction': 3, 'test2': 5}),
        dtype=st.just(float),  # type: ignore
    ),
)
def test_norm(da):
    res = norm(da)
    if not np.isnan(da).any():
        assert (res >= 0).all()
    assert len(res.dims) == len(da.dims) - 1


@given(
    xrst.variables(
        dims=st.just({'test1': 2, 'target': 3, 'test2': 5}),
        dtype=st.just(float),  # type: ignore
    ),
)
def test_subtract_combinations(da):
    assume((da != np.inf).all())
    assume((da != -np.inf).all())
    assume((~np.isnan(da)).all())  # no NaNs allowed
    da = xr.DataArray(da)
    res = subtract_combinations(da, 'target')
    for c, i, j in [(0, 1, 0), (1, 2, 0), (2, 2, 1)]:
        da_diff = da.isel(target=i) - da.isel(target=j)
        to_check = res.isel(targetcomb=c)
        assert_equal(da_diff, to_check)

@given(
    xrst.variables(
        dims=st.just({'test': 2, 'target': 4}),
        dtype=st.just(float),  # type: ignore
    ),
)
def test_pca(da):
    assume(not np.isinf(da).any())  # no +/-inf allowed
    assume(not np.isnan(da).any())  # no NaNs allowed
    res = pca(da, dim='target')
    assert isinstance(res, xr.DataArray) or isinstance(res, xr.Variable)
    assert 'PC' in res.dims


def test_pairwise_dists_pca():
    pass

def test_sudi2(traj_butene):
    da = traj_butene.energy
    left = _sudi_groupby(da)
    right = sudi(da)
    assert_equal(left, right)