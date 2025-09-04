from collections import namedtuple
from functools import wraps

import xarray as xr


from ._generated_accessors import GeneratedDAAcessor, GeneratedDSAcessor

from .core import (
    postprocess as P,
    xrhelpers,
    parse,
    plot,
    filter_unphysical,
    filtre,
    ase,
    geom,
    ml
)
from . import _state
from .core.plot import select, p3mhelpers


CONVERTERS: dict[str, P.Converter] = {
    'convert_energy': P.convert_energy,
    'convert_forces': P.convert_forces, 
    'convert_dipoles': P.convert_dipoles,
    'convert_length': P.convert_length,
}


class ShnitselAccessor:
    _obj: xr.DataArray | xr.Dataset
    _methods: list

    def __init__(self, obj):
        self._obj = obj

    def __dir__(self):
        return self._methods

    def __getattr__(self, key):
        if key in self._potential_methods:
            raise TypeError(
                "This method is unavailable, because the "
                + f"{type(self._obj).__name__!r} object "
                + self._reasons_unavailable(key)
            )
        raise AttributeError(f"{type(self).__name__!r} object has no attribute {key!r}")

    def _reasons_unavailable(self, met):
        reasons = []
        entry = self._potential_methods[met]
        dims = set(self._obj.dims)
        coords = set(self._obj.coords)
        atkeys = set(self._obj.attrs)
        if 'required_vars' in entry._fields:
            vars_ = set(self._obj.data_vars)
            if not vars_ >= (rvars := (entry.required_vars or set())):
                reasons.append(f"is missing required data_vars {rvars - vars_}")
        if not dims >= (rdims := (entry.required_dims or set())):
            reasons.append(f"is missing required dims {rdims - dims}")
        if not coords >= (rcoords := (entry.required_coords or set())):
            reasons.append(f"is missing required coords {rcoords - coords}")
        if entry.required_name is not None and entry.required_name != self._obj.name:
            reasons.append(f"is not named '{entry.required_name}'")
        if not atkeys >= (ratks := (entry.required_attrs or set())):
            reasons.append(f"is missing required attrs {ratks - atkeys}")
        if entry.required_attrs is not None:
            for k, v in entry.required_attrs.items():
                if (actual := self._obj.attrs[k]) != v:
                    reasons.append(
                        f"has attr {k!r} set to value {actual!r} "
                        f"rather than expected {actual!r}"
                    )
        if isect := dims.intersection(entry.incompatible_dims or set()):
            reasons.append(f"has incompatible dims {isect}")
        return "; ".join(reasons)

    def _repr_html_(self):
        available = [f'<li>{met}</li>' for met in self.__dir__()]
        unavailable = [
            # f'<dt>{met}</dt><dd>{self._potential_methods[met]!r}</dd>'
            f"""
                <td>{met}</td>
                <td style='text-align:left'>{self._reasons_unavailable(met)}</td>
            """
            for met in self._potential_methods
            if met not in self.__dir__()
        ]
        return f"""
<div style='display:flex;column-gap:20px;'> 
    <div>
        <b>Available methods:</b>
        <ul>{''.join(available)}</ul>
    </div>
    <div>
        <details>
            <summary><b>Unavailable methods:</b></summary>
            <table>
                <thead>
                    <tr>
                        <th>Method</th>
                        <th style='text-align:left'>Method unavailable because object</th></tr>
                </thead>
                <tbody>
                    <tr>{'</tr><tr>'.join(unavailable)}</tr>
                </tbody>
            </table>
        </details>
    </div>
</div>
        """


@xr.register_dataarray_accessor('sh')
class DAShnitselAccessor(ShnitselAccessor, GeneratedDAAcessor):
    pass


_state.DATAARRAY_ACCESSOR_NAME = 'sh'
_state.DATAARRAY_ACCESSOR_REGISTERED = True

# class DerivedProperties:
#     derivers: dict[str, M2]
#     properties: dict[str, xr.DataArray]
#     groups: dict  # TODO Later!

#     def __init__(self, obj):
#         self._obj = obj

#     def __getitem__(self, *keys):
#         return NotImplemented

#     def keys(self):
#         return set().union(self.derivers, self.properties)

# class DSDerivedProperties(DerivedProperties):
#     derivers = dict(
#         fosc=M2(lambda ds, *a, **k: P.get_fosc(ds.energy, ds.dip_trans, *a, **k)),
#         bond_lengths=M2(lambda ds, *a, **k: geom.get_bond_lengths(ds.atXYZ, *a, **k)),
#         bond_angles=M2(lambda ds, *a, **k: geom.get_bond_angles(ds.atXYZ, *a, **k)),
#         bond_torsions=M2(lambda ds, *a, **k: geom.get_bond_torsions(ds.atXYZ, *a, **k)),
#         bats=M2(lambda ds, *a, **k: geom.get_bats(ds.atXYZ, *a, **k)),
#         pwdists=M2(
#             lambda ds, *a, **k: P.norm(P.subtract_combinations(ds.atXYZ, 'atom'))
#         ),
#     )
#     properties = {}

#     def __getitem__(self, keys):
#         if not isinstance(keys, tuple):
#             keys = (keys,)

#         if len(keys) == 1 and isinstance(keys[0], list):
#             force_ds = True
#             keys = keys[0]
#         else:
#             force_ds = False
#             keys = list(keys)

#         if len(keys) == 1 and not force_ds:
#             k = keys[0]
#             if k in self.derivers:
#                 return self.derivers[k].func(self._obj)
#             elif k in self.properties:
#                 return self.properties[k]
#             else:
#                 return self._obj[k]

#         if Ellipsis in keys:
#             keys.remove(Ellipsis)
#             if Ellipsis in keys:
#                 raise ValueError("Ellipsis ('...') should only be provided once")
#             selection = list(self._obj.data_vars) + keys
#         else:
#             selection = keys

#         to_assign = {
#             k: self.derivers[k].func(self._obj) for k in keys if k in self.derivers
#         }

#         return self._obj.assign(to_assign)[selection]


@xr.register_dataset_accessor('sh')
class DSShnitselAccessor(ShnitselAccessor, GeneratedDSAcessor):
    pass


_state.DATASET_ACCESSOR_NAME = 'sh'
_state.DATASET_ACCESSOR_REGISTERED = True
