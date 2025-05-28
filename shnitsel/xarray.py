from typing import Callable
from functools import partial  # , partialmethod

import xarray as xr
from .dynamic import postprocess


@xr.register_dataarray_accessor('sh')
class AAccessor:
    _potential_methods: dict[str, tuple[set[str], Callable]] = {}

    # @classmethod
    # def _requires(cls, *required_coords):
    #     def decorator(func):
    #         name = func.__name__[1:]  # remove underscore
    #         cls._potential_methods[name] = (set(required_coords), func)
    #         return func

    #     return decorator

    def __init__(self, da):
        self._da = da
        self._methods = self._build_methods()

    def _build_methods(self):
        methods = {}
        coords = set(self._da.coords)
        for name, (required_coords, func) in self._potential_methods.items():
            if required_coords <= coords:
                methods[name] = func
        return methods

    def __getattr__(self, name):
        if name in self._methods:
            return self._methods[name]
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    def __dir__(self):
        return list(super().__dir__()) + list(self._methods.keys())

    # let's try one way:
    # @_requires('atNames', 'atom')
    def _to_xyz(self, comment='#'):
        return postprocess.to_xyz(self._da, comment)

    _potential_methods['to_xyz'] = ({'atNames'}, _to_xyz)

    # and another:
    _potential_methods['dihedral'] = ({'atom'}, partial(postprocess.dihedral))