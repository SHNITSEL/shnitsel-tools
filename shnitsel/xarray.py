from collections import namedtuple
from itertools import chain

import xarray as xr
from .dynamic import postprocess as P, xrhelpers

M = namedtuple(
    'M',
    [
        'func',
        'required_dims',
        'required_coords',
        'required_name',
        'required_attrs',
        'required_attrs_values',
        'incompatible_dims',
    ],
    defaults=[None, None, None, None, None, None, None],
)

DA_METHODS: dict[str, M] = {
    ## From postprocess:
    'norm': M(P.norm),
    'subtract_combinations': M(P.subtract_combinations),
    'pca': M(P.pca),
    'pairwise_dists_pca': M(P.pairwise_dists_pca, required_dims={'atom'}),
    'sudi': M(P.sudi, required_dims={'frame'}),
    'hop_indices': M(P.hop_indices, required_dims={'frame'}, required_name='astate'),
    'relativize': M(P.relativize),
    'convert_energy': M(P.convert_energy, required_attrs={'units'}),
    'convert_dipoles': M(P.convert_dipoles, required_attrs={'units'}),
    'convert_length': M(P.convert_length, required_attrs={'units'}),
    'ts_to_time': M(P.ts_to_time, required_coords={'ts'}),
    'keep_norming': M(P.keep_norming),
    'calc_ci': M(P.xr_calc_ci),  # name differs!
    'to_xyz': M(
        P.to_xyz,
        required_coords={'atNames'},
        incompatible_dims={'frames'},
    ),
    'traj_to_xyz': M(
        P.traj_to_xyz, required_dims={'frame'}, required_coords={'atNames'}
    ),
    'dihedral': M(P.dihedral, required_dims={'atom'}),
    'angle': M(P.angle, required_dims={'atom'}),
    'distance': M(P.distance, required_dims={'atom'}),
    'trajs_with_hops': M(P.trajs_with_hops, required_name='astate'),
    'get_hop_types': M(P.get_hop_types, required_name='astate'),
    # Do not include legacy functions pca_for_plot() or changes(). NB. changes expects astate to be a coord
    # Do not include functions that expect multiple DataArrays: broaden_gauss(), dcross(), ddot(), dnorm()
    # Do not include functions that expect an ndarrayNOT calc_ci NOR ci_agg_last_dim: EXCLUDE because takes ndarray
    #
    ## From xrhelpers:
    'sel_trajs': M(xrhelpers.sel_trajs, {'frame'}),
}

M2 = namedtuple(
    'M2',
    [
        'func',
        'required_vars',
        'required_dims',
        'required_coords',
        'required_name',
        'required_attrs',
        'required_attrs_values',
        'incompatible_dims',
    ],
    defaults=[None, None, None, None, None, None, None],
)

DS_METHODS: dict[str, M2] = {
    'pca_and_hops': M2(P.pca_and_hops, required_vars={'atXYZ', 'astate'}),
    'validate': M2(P.validate),
    'ts_to_time': M2(P.ts_to_time, required_coords={'ts'}),
    'assign_fosc': M2(P.assign_fosc, required_vars={'energy', 'dip_trans'}),
    'broaden_gauss': M2(
        P.ds_broaden_gauss, required_vars={'energy', 'fosc'}
    ),  # name differs!
    'get_per_state': M2(P.get_per_state, required_dims={'state'}),
    'get_inter_state': M2(P.get_inter_state, required_dims={'statecomb'}),
    # calc_pops: should really be a DataArray method; uses only the 'astate' variable, requires 'state' and 'frame' dims
    'time_grouped_ci': M2(
        P.time_grouped_ci, required_coords={'time'}, required_dims={'frame'}
    ),
    'find_hops': M2(P.find_hops, required_coords={'trajid'}, required_vars={'astate'}),
    ## From xrhelpers:
    'save_frames': M2(xrhelpers.save_frames),
    'sel_trajs': M2(xrhelpers.sel_trajs, required_dims={'frame'}),
}


class ShnitselAccessor:
    _obj: xr.DataArray | xr.Dataset
    _methods: list

    def _make_method(self, func):
        def method(*args, **kwargs):
            return func(self._obj, *args, **kwargs)

        try:
            method.__name__ = func.__name__
        except AttributeError:
            pass
        method.__doc__ = func.__doc__
        return method

    def __dir__(self):
        return chain(super().__dir__(), self._methods)


@xr.register_dataarray_accessor('sh')
class DAShnitselAccessor(ShnitselAccessor):
    def __init__(self, da):
        self._obj = da
        self._methods = []

        dims = set(da.dims)
        coords = set(da.coords)
        atkeys = set(da.attrs)
        for name, (
            func,
            rdims,
            rcoords,
            rname,
            ratks,
            ratd,
            xdims,
        ) in DA_METHODS.items():
            if (
                (rdims is None or rdims <= dims)
                and (rcoords is None or rcoords <= coords)
                and (rname is None or rname == da.name)
                and (ratks is None or ratks <= atkeys)
                and (ratd is None or all(ratd[k] == da.attrs[k] for k in ratd))
                and (xdims is None or xdims.isdisjoint(dims))
            ):
                setattr(self, name, self._make_method(func))
                self._methods.append(name)


@xr.register_dataset_accessor('sh')
class DSShnitselAccessor(ShnitselAccessor):
    def __init__(self, ds):
        self._obj = ds
        self._methods = []

        vars = set(ds.data_vars)
        dims = set(ds.dims)
        coords = set(ds.coords)
        atkeys = set(ds.attrs)
        for name, (
            func,
            rvars,
            rdims,
            rcoords,
            rname,
            ratks,
            ratd,
            xdims,
        ) in DS_METHODS.items():
            if (
                (rvars is None or rvars <= vars)
                and (rdims is None or rdims <= dims)
                and (rcoords is None or rcoords <= coords)
                and (rname is None or rname == ds.name)
                and (ratks is None or ratks <= atkeys)
                and (ratd is None or all(ratd[k] == ds.attrs[k] for k in ratd))
                and (xdims is None or (xdims & dims == {}))
            ):
                setattr(self, name, self._make_method(func))
                self._methods.append(name)