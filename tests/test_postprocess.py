import pytest
from hypothesis import assume, given
import hypothesis.strategies as st
from hypothesis.strategies import just, composite
import hypothesis.extra.numpy as hnp
from xarray.testing import assert_equal
import xarray.testing.strategies as xrst

import numpy as np
import xarray as xr

from shnitsel.analyze.generic import (
    norm,
    subtract_combinations,
    keep_norming,
    replace_total,
    relativize,
    pwdists,
)
from shnitsel.analyze.hops import (
    hops,
    focus_hops,
    assign_hop_time,
)
from shnitsel.analyze.populations import calc_classical_populations
from shnitsel.analyze.pca import pca
# from shnitsel.data.helpers import ts_to_time
# from shnitsel.io import read


def dim_name_supersets_of(names):
    required = st.just(names)
    extra = xrst.dimension_names()
    together = st.builds(lambda r, e: r + e, required, extra)

    return together.flatmap(st.permutations)


def arrays_no_nan(dtype, shape):
    if np.issubdtype(dtype, np.floating):
        elements = st.floats(allow_nan=False)
    elif np.issubdtype(dtype, np.complexfloating):
        elements = st.complex_numbers(allow_nan=False)
    else:
        elements = None  # Allow any valid element

    return hnp.arrays(dtype=st.just(dtype), shape=st.just(shape), elements=elements)

class TestProcessing:
    """Class to test all functions of the shnitsel tools related to postprocessing"""

    # TODO: Fixture currently unused. Remove?
    # @pytest.fixture
    # def traj_butene(self):
    #     # TODO: FIXME: This does not work with default settings of the new read() function
    #     frames = read('tutorials/test_data/sharc/traj_butene', kind='sharc')
    #     return ts_to_time(frames)

    @given(
        xrst.variables(
            dims=st.just({'test1': 2, 'direction': 3, 'test2': 5}),
            dtype=st.just(float),  # type: ignore
        ),
    )
    def test_norm(self, da):
        da = xr.DataArray(da)
        res = norm(da)
        if not np.isnan(da).any():
            assert (res >= 0).all()
        assert len(res.dims) == len(da.dims) - 1

    @composite
    def inputs_for_subtract_combinations(draw):
        dims = draw(xrst.dimension_names(min_dims=1))

        dvar = xrst.variables(
            dims=st.just(dims),
            dtype=st.just(float),  # type: ignore
        )
        da = xr.DataArray(draw(dvar))

        axis = draw(st.integers(min_value=0, max_value=len(dims) - 1))
        target_dim = dims[axis]

        return da, target_dim

    @given(inputs_for_subtract_combinations())
    def test_subtract_combinations(self, inputs):
        from itertools import combinations

        da, target_dim = inputs
        assume((da != np.inf).all())
        assume((da != -np.inf).all())
        assume((~np.isnan(da)).all())  # no NaNs allowed
        da = xr.DataArray(da)
        res = subtract_combinations(da, target_dim)
        combs = combinations(range(da.sizes[target_dim]), 2)
        for c, (i, j) in enumerate(combs):
            da_diff = da[{target_dim: j}] - da.isel({target_dim: i})
            to_check = res[{target_dim + "comb": c}]
            assert_equal(da_diff, to_check)

    # TODO: test center

    @composite
    def inputs_for_keep_norming(draw):
        dims = draw(xrst.dimension_names(min_dims=1))
        dims_in_var = draw(xrst.unique_subset_of(dims))
        dims_to_exclude = draw(xrst.unique_subset_of(dims))
        da = draw(xrst.variables(dims=st.just(dims_in_var), dtype=st.just(float)))
        return da, dims_to_exclude

    @given(inputs_for_keep_norming())
    def test_keep_norming(self, inputs):
        # NB. The output of keep_norming isn't necessarily >= 0
        da, exclude = inputs
        da = xr.DataArray(da)
        res = keep_norming(da, exclude=exclude)
        assert set(da.dims) - set(res.dims) == set(da.dims) - set(exclude)

    @composite
    def inputs_for_replace_total(draw):
        def make_array_strategy(dtype, shape):
            if np.issubdtype(dtype, np.floating):
                elements = st.floats(allow_nan=False)
            elif np.issubdtype(dtype, np.complexfloating):
                elements = st.complex_numbers(allow_nan=False)
            else:
                elements = None  # Allow any valid element

            return hnp.arrays(
                dtype=st.just(dtype), shape=st.just(shape), elements=elements
            )

        da = draw(xrst.variables(array_strategy_fn=make_array_strategy))
        to_replace = np.unique(da)
        nitems = len(to_replace)
        # TODO: Make replace_total() robust against mixtures of numbers and strings
        kws = dict(min_size=nitems, max_size=nitems, unique=True)
        value = draw(
            st.lists(st.integers(), **kws),
            st.lists(st.floats(), **kws),
            st.lists(st.text(), **kws),
        )
        return da, to_replace, value

    @given(inputs_for_replace_total())
    def test_replace_total(self, inputs):
        """Test round trip -- replace and replace again"""
        # FIXME: the function being tested expects inputs that
        # support `=` and `<`; this exculdes nan and mixtures of
        # string and int
        da, to_replace, value = inputs
        da = xr.DataArray(da)
        da2 = replace_total(da, to_replace, value)
        da3 = replace_total(da2, to_replace=value, value=to_replace)
        assert_equal(da, da3)

    @given(xrst.variables(dtype=st.one_of(st.just(float))))
    def test_relativize(self, da):
        # TODO: test **sel parameter
        da = xr.DataArray(da)
        res = relativize(da)
        # FIXME: Fails for extremely large numbers on some platforms. Change test? Fix function? Document?
        if (~(np.isnan(res) | np.isinf(res))).any():
            assert res.min() == 0

    @pytest.mark.xfail
    @given(xrst.variables(dims=dim_name_supersets_of(['atom', 'direction'])))
    def test_pwdists(self, da):
        da = xr.DataArray(da)
        res = pwdists(da)

        assert set(da.dims) - set(res.dims) == {'atom', 'direction'}
        assert set(res.dims) - set(da.dims) == {'atomcomb'}
        if not np.isnan(da).any() and not np.isinf(da).any():
            assert (res >= 0).all()

    ##############
    # analyze.hops
    @composite
    def frames_for_hops(draw):
        nframes = draw(hnp.array_shapes(min_dims=1, max_dims=1))[0]

        def get_var(dtype):
            return draw(
                xrst.variables(dims=st.just({'frame': nframes}), dtype=st.just(dtype))
            )

        astate = get_var(int)
        time = get_var(float)
        trajid = get_var(int)

        da = xr.DataArray(
            astate, coords={'time': ('frame', time), 'trajid': ('frame', trajid)}
        ).set_xindex(['trajid', 'time'])
        return xr.Dataset({'astate': da})
    
    @pytest.mark.xfail
    @given(frames_for_hops())
    def test_hops(self, frames):
        res = hops(frames)
        assert 'tidx' in res
        assert 'hop_from' in res
        assert 'hop_to' in res

    @pytest.mark.xfail
    @given(frames_for_hops())
    def test_focus_hops(self, frames):
        res = focus_hops(frames)
        # Check hop-independent coord dimensions
        assert res['hop_time'].dims == ('hop_time',)
        assert res['hop_tidx'].dims == ('hop_time',)

        # Check per-hop 1D coord dimensions
        assert res['hop_from'].dims == ('hop',)
        assert res['hop_to'].dims == ('hop',)
        assert res['trajid'].dims == ('hop',)

        # Check per-hop 2D coord dimensions
        assert set(res['time'].dims) == {'hop', 'hop_time'}
        assert set(res['tidx'].dims) == {'hop', 'hop_time'}

    @pytest.mark.xfail
    @given(frames_for_hops(), st.booleans())
    def test_assign_hop_time(self, frames, choose_first):
        which = 'first' if choose_first else 'last'
        res = assign_hop_time(frames, which=which)
        assert 'hop_time' in res.coords
        assert res.coords['hop_time'].dims == ('frame',)

    #############
    # Populations

    @composite
    def frames_for_populations(draw):
        nframes = draw(hnp.array_shapes(min_dims=1, max_dims=1))[0]

        # calc_classical_populations expects state IDs to be integers starting at 1 without gaps.
        # Does it make sense to lift that assumption?
        # Previous approach:
        # states = draw(st.lists(st.integers(min_value=1), min_size=1, unique=True))
        states = range(1, draw(st.integers(min_value=2, max_value=20)))

        def make_array_strategy(dtype, shape):
            assert dtype is int
            return hnp.arrays(
                dtype=st.just(dtype),
                shape=st.just(shape),
                elements=st.sampled_from(states),
            )

        astate = draw(
            xrst.variables(
                array_strategy_fn=make_array_strategy,
                dims=st.just({'frame': nframes}),
                dtype=st.just(int),
            )
        )

        def get_var(dtype):
            return draw(
                xrst.variables(dims=st.just({'frame': nframes}), dtype=st.just(dtype))
            )

        time = get_var(float)
        trajid = get_var(int)

        da = xr.DataArray(
            astate, coords={'time': ('frame', time), 'trajid': ('frame', trajid)}
        ).set_xindex(['trajid', 'time'])
        return xr.Dataset({'astate': da}, coords={'state': states})

    @pytest.mark.xfail
    @given(frames_for_populations())
    def test_calc_classical_populations(self, frames):
        res = calc_classical_populations(frames)
        assert 'state' in res.dims
        assert ((0 <= res) & (res <= 1)).all()

    #################################
    # Dimensional reduction functions
    @pytest.mark.xfail
    @given(
        xrst.variables(
            dims=st.just({'test': 2, 'target': 4}),
            dtype=st.just(float),  # type: ignore
        ),
    )
    def test_pca(self, da):
        assume(not np.isinf(da).any())  # no +/-inf allowed
        assume(not np.isnan(da).any())  # no NaNs allowed
        da = xr.DataArray(da)

        res = pca(da, dim='target')
        assert isinstance(res, xr.DataArray) or isinstance(res, xr.Variable)
        assert 'PC' in res.dims

    def test_pairwise_dists_pca(self):
        pass
