from typing import Callable
from functools import partial  # , partialmethod

import xarray as xr
from .dynamic import postprocess, xrhelpers

DA_METHODS: dict[str, tuple[Callable, set[str]]] = {
    'to_xyz': (postprocess.to_xyz, {'atNames'}),
    'dihedral': (postprocess.dihedral, {'atom'}),
    'sudi': (postprocess.sudi, set()),
    'sel_trajs': (xrhelpers.sel_trajs, {'frame'}),
}


@xr.register_dataarray_accessor('sh')
class DAShnitselAccessor:
    def __init__(self, da):
        self._da = da

        coords = set(da.coords)
        for name, (func, required_coords) in DA_METHODS.items():
            if required_coords <= coords:
                setattr(self, name, self._make_method(func))

    def _make_method(self, func):
        def method(*args, **kwargs):
            return func(self._da, *args, **kwargs)

        method.__name__ = func.__name__
        method.__doc__ = func.__doc__
        return method

    # def __getattr__(self, name):
    #     if name in self._methods:
    #         return self._methods[name]
    #     raise AttributeError(
    #         f"'{type(self).__name__}' object has no attribute '{name}'"
    #     )

    # def __dir__(self):
    #     return list(super().__dir__()) + list(self._methods.keys())

    # let's try one way:
    # @_requires('atNames', 'atom')
    def _to_xyz(self, comment='#'):
        return postprocess.to_xyz(self._da, comment)

    def subtract_combinations(self, dim, labels=False):
        return postprocess.subtract_combinations(self._da, dim, labels)