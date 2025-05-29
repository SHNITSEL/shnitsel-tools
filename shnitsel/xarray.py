from collections import namedtuple

import xarray as xr
from .dynamic import postprocess, xrhelpers

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
    'norm': M(postprocess.norm),
    'subtract_combinations': M(postprocess.subtract_combinations),
    'pca': M(postprocess.pca),
    'pairwise_dists_pca': M(postprocess.pairwise_dists_pca, required_dims={'atom'}),
    'sudi': M(postprocess.sudi, required_dims={'frame'}),
    'hop_indices': M(
        postprocess.hop_indices, required_dims={'frame'}, required_name='astate'
    ),
    'relativize': M(postprocess.relativize),
    'convert_energy': M(postprocess.convert_energy, required_attrs={'units'}),
    'convert_dipoles': M(postprocess.convert_dipoles, required_attrs={'units'}),
    'convert_length': M(postprocess.convert_length, required_attrs={'units'}),
    'ts_to_time': M(postprocess.ts_to_time, required_coords={'ts'}),
    'keep_norming': M(postprocess.keep_norming),
    'xr_calc_ci': M(postprocess.xr_calc_ci),
    'to_xyz': M(
        postprocess.to_xyz,
        required_coords={'atNames'},
        incompatible_dims={'frames'},
    ),
    'traj_to_xyz': M(
        postprocess.traj_to_xyz, required_dims={'frame'}, required_coords={'atNames'}
    ),
    'dnorm': M(postprocess.dnorm, required_dims={'direction'}),
    'dcross': M(postprocess.dcross, required_dims={'direction'}),
    'ddot': M(postprocess.ddot, required_dims={'direction'}),
    'dihedral': M(postprocess.dihedral, required_dims={'atom'}),
    'angle': M(postprocess.angle, required_dims={'atom'}),
    'distance': M(postprocess.distance, required_dims={'atom'}),
    'trajs_with_hops': M(postprocess.trajs_with_hops, required_name='astate'),
    'get_hop_types': M(postprocess.get_hop_types, required_name='astate'),
    # NOT pca_for_plot: da, legacy: leave out
    # NOT changes -- this is a legacy equivalent to hop_indices; notably, it assumes that 'astates' is a coord
    # NOT broaden_gauss: EXCLUDE because requires two DataArrays as input -- maybe let it get data from a coord  instead if one's 1D?
    # NOT calc_ci NOR ci_agg_last_dim: EXCLUDE because takes ndarray
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
    'pca_and_hops': M2(postprocess.pca_and_hops, required_vars={'atXYZ', 'astate'}),
    'validate': M2(postprocess.validate),
    'ts_to_time': M2(postprocess.ts_to_time, required_coords={'ts'}),
    'assign_fosc': M2(postprocess.assign_fosc, required_vars={'energy', 'dip_trans'}),
    # ds_broaden_gauss: requires 'energy' and 'fosc' -- but should it be so specific? Aren't there other thing one might want to broaden? Also should 'ds' stay in the name?
    # get_per_state: only really makes sense with at least one of {'energy', 'forces', 'dip_perm'}
    # get_inter_state: oh, just include it unconditionally
    # calc_pops: should really be a DataArray method; uses only the 'astate' variable, requires 'state' and 'frame' dims
    'time_grouped_ci': M2(
        postprocess.time_grouped_ci, required_coords={'time'}, required_dims={'frame'}
    ),
    'find_hops': M2(
        postprocess.find_hops, required_coords={'trajid'}, required_vars={'astate'}
    ),
    ## From xrhelpers:
    'save_frames': M2(xrhelpers.save_frames),
}


class ShnitselAccessor:
    _obj: xr.DataArray | xr.Dataset

    def _make_method(self, func):
        def method(*args, **kwargs):
            return func(self._obj, *args, **kwargs)

        try:
            method.__name__ = func.__name__
        except AttributeError:
            pass
        method.__doc__ = func.__doc__
        return method

@xr.register_dataarray_accessor('sh')
class DAShnitselAccessor(ShnitselAccessor):
    def __init__(self, da):
        self._obj = da

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
                and (xdims is None or len(xdims & dims) == 0)
            ):
                setattr(self, name, self._make_method(func))


@xr.register_dataset_accessor('sh')
class DSShnitselAccessor(ShnitselAccessor):
    def __init__(self, ds):
        self._obj = ds

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