from typing import Sequence

import xarray as xr

from shnitsel._contracts import Needs
# from shnitsel.data.shnitsel_db_helpers import concat_subtree

from shnitsel.filtering.structure_selection import (
    StructureSelection,
    FeatureLevelOptions,
)

# CONVERTERS: dict[str, P.Converter] = {
#     'convert_energy': P.convert_energy,
#     'convert_forces': P.convert_forces,
#     'convert_dipoles': P.convert_dipoles,
#     'convert_length': P.convert_length,
# }

DATASET_ACCESSOR_NAME = 'st'
DATAARRAY_ACCESSOR_NAME = 'st'
DATATREE_ACCESSOR_NAME = 'st'


class ShnitselAccessor:
    _obj: xr.DataArray | xr.Dataset
    _methods: list

    def __init__(self, obj):
        self._obj = obj
        self._method_needs = {
            met: getattr(getattr(self, met), '_needs', Needs()) for met in self._methods
        }
        self.suitable = []
        self.unsuitable = {}
        for met in self._method_needs:
            reasons = self._reasons_unavailable(met)
            if reasons:
                self.unsuitable[met] = reasons
            else:
                self.suitable.append(met)

    def _reasons_unavailable(self, met):
        reasons = []
        entry = self._method_needs[met]
        dims = set(self._obj.dims)
        coords = set(self._obj.coords)
        # atkeys = set(self._obj.attrs)
        if 'data_vars' in entry._fields:
            vars_ = set(getattr(self._obj, 'data_vars', []))
            if not vars_ >= (rvars := (entry.data_vars or set())):
                reasons.append(f"is missing required data_vars {rvars - vars_}")
        if not dims >= (rdims := (entry.dims or set())):
            reasons.append(f"is missing required dims {rdims - dims}")
        if not coords >= (rcoords := (entry.coords or set())):
            reasons.append(f"is missing required coords {rcoords - coords}")
        if entry.name is not None and entry.name != self._obj.name:
            reasons.append(f"is not named '{entry.required_name}'")
        # if not atkeys >= (ratks := (entry.required_attrs or set())):
        #     reasons.append(f"is missing required attrs {ratks - atkeys}")
        # if entry.required_attrs is not None:
        #     for k, v in entry.required_attrs.items():
        #         if (actual := self._obj.attrs[k]) != v:
        #             reasons.append(
        #                 f"has attr {k!r} set to value {actual!r} "
        #                 f"rather than expected {actual!r}"
        #             )
        if isect := dims.intersection(entry.not_dims or set()):
            reasons.append(f"has incompatible dims {isect}")
        if reasons:
            return "; ".join(reasons)
        else:
            return ""

    def _repr_html_(self):
        unavailable = [
            f"""
                <td>{met}</td>
                <td style='text-align:left'>{reasons}</td>
            """
            for met, reasons in self.unsuitable.items()
        ]
        ustr = f"""
        <details>
            <summary><b>Unavailable methods:</b></summary>
            <table>
                <thead>
                <tr>
                    <th>Method</th>
                    <th style='text-align:left'>Method unavailable because object</th>
                </tr>
                </thead>
            <tbody>
                <tr>{'</tr><tr>'.join(unavailable)}</tr>
            </tbody>
            </table>
        </details>
        
        """
        available = [f'<li>{met}</li>' for met in self.suitable]
        astr = f"""
        <div>
            <b>Available methods:</b>
            <ul>{''.join(available)}</ul>
        </div>
        """
        return f"""
        <div style='display:flex;column-gap:20px;'> 
            <div>
                {astr}
            </div>
            <div>
                {ustr}
            </div>
        </div>
        """

    def _ipython_key_completions_(self) -> list[str]:
        """Provide method for the key-autocompletions in IPython.
        See https://ipython.readthedocs.io/en/stable/config/integrating.html#tab-completion
        For the details.
        """

        return list(self.suitable)


class DAManualAccessor(ShnitselAccessor):
    pass


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


class DSManualAccessor(ShnitselAccessor):
    def struc_sel(
        self,
        default_selection: Sequence[FeatureLevelOptions] | None = None,
        to2D: bool = True,
        frame_index: int = 0,
    ) -> StructureSelection:
        """Create an initial StructureSelection object from a dataset using the entire structural information in it.

        Parameters
        ----------
        frame (xr.Dataset)
            Should only represent a single frame of data (i.e. no 'frame' dimension or of zero-size).
        default_selection (Sequence[FeatureLevelOptions], optional)
            List of features to activate as selected by default. Defaults to [ 'atoms', 'bonds', ].
        to2D (bool, optional)
            Flag to control whether a mol representation is converted to a 2d projection before use for visualization.
        frame_index
            Used to select a single frame (using ``.isel()``) if an object containing multiple frames is provided;
            by default, the first frame is used.

        Raises
        ------
        ValueError
            If no structural information could be extracted from the dataset

        Returns
        -------
        StructureSelection
            A structure selection object initially covering all atoms and structural features.

        Notes
        -----
        The dataset (single frame) to extract the structure information out of must have
        at least an `atXYZ` variable and a `atom` dimension.
        Ideally, an `atom` coordinate for feature selection is also provided.
        """

        default_selection = default_selection or ['atoms', 'bonds']

        if 'frame' in self._obj.dims:
            frame = self._obj.isel(frame=frame_index)
        elif 'time' in self._obj.dims:
            frame = self._obj.isel(time=frame_index)
        else:
            frame = self.obj_

        return StructureSelection.init_from_dataset(
            frame,
            default_selection=default_selection,
            to2D=to2D,
        )
