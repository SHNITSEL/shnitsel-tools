import collections
import numpy
import numpy.typing as npt
import sklearn
import typing
import xarray
import xarray as xr
from ._accessors import DAManualAccessor, DSManualAccessor
from ._contracts import needs
from numpy import ndarray
from rdkit.Chem.rdchem import Mol
from shnitsel.core.ase import write_ase
from shnitsel.core.filtre import energy_filtranda, get_cutoffs, last_time_where, truncate
from shnitsel.core.geom import get_bats, get_bond_angles, get_bond_lengths, get_bond_torsions, kabsch
from shnitsel.core.ml import lda, pca, pls, pls_ds
from shnitsel.core.parse.sharc_icond import iconds_to_frames
from shnitsel.core.plot.p3mhelpers import frame3D, frames3Dgrid, traj3D, trajs3Dgrid
from shnitsel.core.plot.select import FrameSelector, TrajSelector
from shnitsel.core.plot.spectra3d import spectra_all_times
from shnitsel.core.postprocess import angle, assign_fosc, broaden_gauss, calc_ci, calc_pops, convert_dipoles, convert_energy, convert_forces, convert_length, default_mol, dihedral, distance, find_hops, get_hop_types, get_inter_state, get_per_state, hop_indices, keep_norming, norm, pairwise_dists_pca, pca_and_hops, relativize, setup_frames, smiles_map, subtract_combinations, sudi, time_grouped_ci, to_mol, to_xyz, traj_to_xyz, trajs_with_hops, ts_to_time, validate
from shnitsel.core.xrhelpers import assign_levels, expand_midx, flatten_levels, mgroupby, msel, save_frames, sel_trajids, sel_trajs, stack_trajs, unstack_trajs
from typing import Dict, Hashable, List, Literal, Optional, Sequence, Union
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.groupby import DataArrayGroupBy, DatasetGroupBy


class DataArrayAccessor(DAManualAccessor):
    _methods = [
        'norm',
        'subtract_combinations',
        'pairwise_dists_pca',
        'sudi',
        'hop_indices',
        'relativize',
        'ts_to_time',
        'keep_norming',
        'calc_ci',
        'time_grouped_ci',
        'to_xyz',
        'traj_to_xyz',
        'dihedral',
        'angle',
        'distance',
        'trajs_with_hops',
        'get_hop_types',
        'to_mol',
        'smiles_map',
        'default_mol',
        'convert_energy',
        'convert_forces',
        'convert_dipoles',
        'convert_length',
        'flatten_levels',
        'expand_midx',
        'assign_levels',
        'mgroupby',
        'msel',
        'sel_trajs',
        'sel_trajids',
        'last_time_where',
        'get_bond_lengths',
        'get_bond_angles',
        'get_bond_torsions',
        'get_bats',
        'kabsch',
        'FrameSelector',
        'TrajSelector',
        'frame3D',
        'frames3Dgrid',
        'traj3D',
        'trajs3Dgrid',
        'pca',
        'lda',
        'pls',
    ]

    def norm(self, dim: Hashable='direction', keep_attrs: bool | str | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.norm`."""
        return norm(self._obj, dim=dim, keep_attrs=keep_attrs)

    def subtract_combinations(self, dim: Hashable, labels: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.subtract_combinations`."""
        return subtract_combinations(self._obj, dim, labels=labels)

    @needs(dims={'atom'})
    def pairwise_dists_pca(self, **kwargs) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.pairwise_dists_pca`."""
        return pairwise_dists_pca(self._obj, **kwargs)

    @needs(dims={'frame'})
    def sudi(self) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.sudi`."""
        return sudi(self._obj)

    def hop_indices(self) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.hop_indices`."""
        return hop_indices(self._obj)

    def relativize(self, **sel) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.relativize`."""
        return relativize(self._obj, **sel)

    @needs(coords={'ts'})
    def ts_to_time(self, delta_t: float | None=None, old: Literal='drop') -> xarray.core.dataset.Dataset | xarray.core.dataarray.DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.ts_to_time`."""
        return ts_to_time(self._obj, delta_t=delta_t, old=old)

    def keep_norming(self, exclude: Optional=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.keep_norming`."""
        return keep_norming(self._obj, exclude=exclude)

    def calc_ci(self, confidence: float=0.95) -> ndarray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.calc_ci`."""
        return calc_ci(self._obj, confidence=confidence)

    @needs(dims={'frame'}, groupable={'time'})
    def time_grouped_ci(self, confidence: float=0.9) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.postprocess.time_grouped_ci`."""
        return time_grouped_ci(self._obj, confidence=confidence)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def to_xyz(self, comment='#') -> str:
        """Wrapper for :py:func:`shnitsel.core.postprocess.to_xyz`."""
        return to_xyz(self._obj, comment=comment)

    @needs(dims={'atom', 'direction'}, groupable={'time'}, coords_or_vars={'atNames'})
    def traj_to_xyz(self) -> str:
        """Wrapper for :py:func:`shnitsel.core.postprocess.traj_to_xyz`."""
        return traj_to_xyz(self._obj)

    @needs(dims={'atom'})
    def dihedral(self, i: int, j: int, k: int, l: int, deg: bool=False, full: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.dihedral`."""
        return dihedral(self._obj, i, j, k, l, deg=deg, full=full)

    @needs(dims={'atom'})
    def angle(self, i: int, j: int, k: int, deg: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.angle`."""
        return angle(self._obj, i, j, k, deg=deg)

    @needs(dims={'atom'})
    def distance(self, i: int, j: int) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.distance`."""
        return distance(self._obj, i, j)

    def trajs_with_hops(self) -> list:
        """Wrapper for :py:func:`shnitsel.core.postprocess.trajs_with_hops`."""
        return trajs_with_hops(self._obj)

    def get_hop_types(self) -> dict:
        """Wrapper for :py:func:`shnitsel.core.postprocess.get_hop_types`."""
        return get_hop_types(self._obj)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def to_mol(self, charge=None, covFactor=1.2, to2D=True, molAtomMapNumber=None, atomNote=None, atomLabel=None) -> Mol:
        """Wrapper for :py:func:`shnitsel.core.postprocess.to_mol`."""
        return to_mol(self._obj, charge=charge, covFactor=covFactor, to2D=to2D, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def smiles_map(self, charge=0, covFactor=1.5) -> str:
        """Wrapper for :py:func:`shnitsel.core.postprocess.smiles_map`."""
        return smiles_map(self._obj, charge=charge, covFactor=covFactor)

    def default_mol(self) -> Mol:
        """Wrapper for :py:func:`shnitsel.core.postprocess.default_mol`."""
        return default_mol(self._obj)

    def convert_energy(self, to: str):
        """Wrapper for :py:func:`shnitsel.core.postprocess.convert_energy`."""
        return convert_energy(self._obj, to)

    def convert_forces(self, to: str):
        """Wrapper for :py:func:`shnitsel.core.postprocess.convert_forces`."""
        return convert_forces(self._obj, to)

    def convert_dipoles(self, to: str):
        """Wrapper for :py:func:`shnitsel.core.postprocess.convert_dipoles`."""
        return convert_dipoles(self._obj, to)

    def convert_length(self, to: str):
        """Wrapper for :py:func:`shnitsel.core.postprocess.convert_length`."""
        return convert_length(self._obj, to)

    def flatten_levels(self, idx_name: str, levels: Sequence[str], new_name: str | None=None, position: int=0, renamer: typing.Callable | None=None) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.flatten_levels`."""
        return flatten_levels(self._obj, idx_name, levels, new_name=new_name, position=position, renamer=renamer)

    def expand_midx(self, midx_name, level_name, value) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.expand_midx`."""
        return expand_midx(self._obj, midx_name, level_name, value)

    def assign_levels(self, levels: dict[str, npt.ArrayLike] | None=None, **levels_kwargs: npt.ArrayLike) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.assign_levels`."""
        return assign_levels(self._obj, levels=levels, **levels_kwargs)

    def mgroupby(self, levels: Sequence[str]) -> DataArrayGroupBy | DatasetGroupBy:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.mgroupby`."""
        return mgroupby(self._obj, levels)

    def msel(self, **kwargs) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.msel`."""
        return msel(self._obj, **kwargs)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert=False) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajids(self, trajids: npt.ArrayLike, invert=False) -> xr.Dataset:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.sel_trajids`."""
        return sel_trajids(self._obj, trajids, invert=invert)

    @needs(dims={'frame'}, coords={'time', 'trajid'})
    def last_time_where(self):
        """Wrapper for :py:func:`shnitsel.core.filtre.last_time_where`."""
        return last_time_where(self._obj)

    @needs(dims={'atom'})
    def get_bond_lengths(self, bond_types=None, mol=None):
        """Wrapper for :py:func:`shnitsel.core.geom.get_bond_lengths`."""
        return get_bond_lengths(self._obj, bond_types=bond_types, mol=mol)

    @needs(dims={'atom'})
    def get_bond_angles(self, angle_types=None, mol=None, deg=False):
        """Wrapper for :py:func:`shnitsel.core.geom.get_bond_angles`."""
        return get_bond_angles(self._obj, angle_types=angle_types, mol=mol, deg=deg)

    @needs(dims={'atom'})
    def get_bond_torsions(self, quadruple_types=None, mol=None, signed=False, deg=False):
        """Wrapper for :py:func:`shnitsel.core.geom.get_bond_torsions`."""
        return get_bond_torsions(self._obj, quadruple_types=quadruple_types, mol=mol, signed=signed, deg=deg)

    @needs(dims={'atom'})
    def get_bats(self, mol=None, signed=False, deg=False):
        """Wrapper for :py:func:`shnitsel.core.geom.get_bats`."""
        return get_bats(self._obj, mol=mol, signed=signed, deg=deg)

    @needs(dims={'atom', 'direction'})
    def kabsch(self, reference_or_indexers: xarray.core.dataarray.DataArray | dict | None=None, **indexers_kwargs):
        """Wrapper for :py:func:`shnitsel.core.geom.kabsch`."""
        return kabsch(self._obj, reference_or_indexers=reference_or_indexers, **indexers_kwargs)

    def FrameSelector(self, xname=None, yname=None, title='', allowed_ws_origin=None, webgl=True):
        """Wrapper for :py:func:`shnitsel.core.plot.select.FrameSelector`."""
        return FrameSelector(self._obj, xname=xname, yname=yname, title=title, allowed_ws_origin=allowed_ws_origin, webgl=webgl)

    def TrajSelector(self, xname=None, yname=None, title='', allowed_ws_origin=None, webgl=True):
        """Wrapper for :py:func:`shnitsel.core.plot.select.TrajSelector`."""
        return TrajSelector(self._obj, xname=xname, yname=yname, title=title, allowed_ws_origin=allowed_ws_origin, webgl=webgl)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def frame3D(self):
        """Wrapper for :py:func:`shnitsel.core.plot.p3mhelpers.frame3D`."""
        return frame3D(self._obj)

    @needs(dims={'atom', 'direction'}, groupable={'frame'}, coords_or_vars={'atNames'})
    def frames3Dgrid(self):
        """Wrapper for :py:func:`shnitsel.core.plot.p3mhelpers.frames3Dgrid`."""
        return frames3Dgrid(self._obj)

    @needs(dims={'atom', 'direction'}, groupable={'time'}, coords_or_vars={'atNames'})
    def traj3D(self):
        """Wrapper for :py:func:`shnitsel.core.plot.p3mhelpers.traj3D`."""
        return traj3D(self._obj)

    @needs(dims={'atom', 'direction'}, coords={'trajid'}, groupable={'time'}, coords_or_vars={'atNames'})
    def trajs3Dgrid(self, trajids: list[int | str] | None=None, loop='forward'):
        """Wrapper for :py:func:`shnitsel.core.plot.p3mhelpers.trajs3Dgrid`."""
        return trajs3Dgrid(self._obj, trajids=trajids, loop=loop)

    def pca(self, dim: str, n_components: int=2, return_pca_object: bool=False) -> tuple[xarray.core.dataarray.DataArray, sklearn.decomposition._pca.PCA] | xarray.core.dataarray.DataArray:
        """Wrapper for :py:func:`shnitsel.core.ml.pca`."""
        return pca(self._obj, dim, n_components=n_components, return_pca_object=return_pca_object)

    def lda(self, dim, cats, n_components=2):
        """Wrapper for :py:func:`shnitsel.core.ml.lda`."""
        return lda(self._obj, dim, cats, n_components=n_components)

    def pls(self, yda, n_components=2, common_dim=None):
        """Wrapper for :py:func:`shnitsel.core.ml.pls`."""
        return pls(self._obj, yda, n_components=n_components, common_dim=common_dim)


class DatasetAccessor(DSManualAccessor):
    _methods = [
        'pca_and_hops',
        'validate',
        'ts_to_time',
        'setup_frames',
        'assign_fosc',
        'broaden_gauss',
        'get_per_state',
        'get_inter_state',
        'calc_pops',
        'find_hops',
        'default_mol',
        'flatten_levels',
        'expand_midx',
        'assign_levels',
        'mgroupby',
        'msel',
        'save_frames',
        'sel_trajs',
        'unstack_trajs',
        'stack_trajs',
        'iconds_to_frames',
        'spectra_all_times',
        'energy_filtranda',
        'get_cutoffs',
        'truncate',
        'write_ase',
        'pls_ds',
    ]

    @needs(coords_or_vars={'astate', 'atXYZ'})
    def pca_and_hops(self) -> tuple:
        """Wrapper for :py:func:`shnitsel.core.postprocess.pca_and_hops`."""
        return pca_and_hops(self._obj)

    def validate(self) -> ndarray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.validate`."""
        return validate(self._obj)

    @needs(coords={'ts'})
    def ts_to_time(self, delta_t: float | None=None, old: Literal='drop') -> xarray.core.dataset.Dataset | xarray.core.dataarray.DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.ts_to_time`."""
        return ts_to_time(self._obj, delta_t=delta_t, old=old)

    def setup_frames(self, to_time: bool | None=None, convert_to_eV: bool | None=None, convert_e_kin_to_eV: bool | None=None, relativize_energy: bool | None=None, relativize_selector=None) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.postprocess.setup_frames`."""
        return setup_frames(self._obj, to_time=to_time, convert_to_eV=convert_to_eV, convert_e_kin_to_eV=convert_e_kin_to_eV, relativize_energy=relativize_energy, relativize_selector=relativize_selector)

    @needs(data_vars={'dip_trans', 'energy'})
    def assign_fosc(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.postprocess.assign_fosc`."""
        return assign_fosc(self._obj)

    @needs(data_vars={'energy', 'fosc'})
    def broaden_gauss(self, fosc: DataArray, agg_dim: Hashable='frame', width: float=0.5, nsamples: int=1000, xmax: float | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.broaden_gauss`."""
        return broaden_gauss(self._obj, fosc, agg_dim=agg_dim, width=width, nsamples=nsamples, xmax=xmax)

    @needs(dims={'state'})
    def get_per_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.postprocess.get_per_state`."""
        return get_per_state(self._obj)

    def get_inter_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.postprocess.get_inter_state`."""
        return get_inter_state(self._obj)

    @needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
    def calc_pops(self) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.postprocess.calc_pops`."""
        return calc_pops(self._obj)

    @needs(coords={'trajid'}, data_vars={'astate'})
    def find_hops(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.postprocess.find_hops`."""
        return find_hops(self._obj)

    def default_mol(self) -> Mol:
        """Wrapper for :py:func:`shnitsel.core.postprocess.default_mol`."""
        return default_mol(self._obj)

    def flatten_levels(self, idx_name: str, levels: Sequence[str], new_name: str | None=None, position: int=0, renamer: typing.Callable | None=None) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.flatten_levels`."""
        return flatten_levels(self._obj, idx_name, levels, new_name=new_name, position=position, renamer=renamer)

    def expand_midx(self, midx_name, level_name, value) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.expand_midx`."""
        return expand_midx(self._obj, midx_name, level_name, value)

    def assign_levels(self, levels: dict[str, npt.ArrayLike] | None=None, **levels_kwargs: npt.ArrayLike) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.assign_levels`."""
        return assign_levels(self._obj, levels=levels, **levels_kwargs)

    def mgroupby(self, levels: Sequence[str]) -> DataArrayGroupBy | DatasetGroupBy:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.mgroupby`."""
        return mgroupby(self._obj, levels)

    def msel(self, **kwargs) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.msel`."""
        return msel(self._obj, **kwargs)

    def save_frames(self, path, complevel=9):
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.save_frames`."""
        return save_frames(self._obj, path, complevel=complevel)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert=False) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    def unstack_trajs(self) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.unstack_trajs`."""
        return unstack_trajs(self._obj)

    def stack_trajs(self) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.xrhelpers.stack_trajs`."""
        return stack_trajs(self._obj)

    @needs(dims={'icond'}, coords={'icond'}, not_dims={'time'})
    def iconds_to_frames(self):
        """Wrapper for :py:func:`shnitsel.core.parse.sharc_icond.iconds_to_frames`."""
        return iconds_to_frames(self._obj)

    @needs(coords={'frame', 'trajid'}, data_vars={'energy', 'fosc'})
    def spectra_all_times(self):
        """Wrapper for :py:func:`shnitsel.core.plot.spectra3d.spectra_all_times`."""
        return spectra_all_times(self._obj)

    @needs(data_vars={'e_kin', 'energy'})
    def energy_filtranda(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.filtre.energy_filtranda`."""
        return energy_filtranda(self._obj)

    @needs(dims={'frame'}, coords={'time', 'trajid'})
    def get_cutoffs(self):
        """Wrapper for :py:func:`shnitsel.core.filtre.get_cutoffs`."""
        return get_cutoffs(self._obj)

    @needs(dims={'frame'}, coords={'time', 'trajid'})
    def truncate(self, cutoffs):
        """Wrapper for :py:func:`shnitsel.core.filtre.truncate`."""
        return truncate(self._obj, cutoffs)

    @needs(dims={'frame'})
    def write_ase(self, db_path: str, kind: str | None, keys: Optional=None, preprocess: bool=True):
        """Wrapper for :py:func:`shnitsel.core.ase.write_ase`."""
        return write_ase(self._obj, db_path, kind, keys=keys, preprocess=preprocess)

    def pls_ds(self, xname, yname, n_components=2):
        """Wrapper for :py:func:`shnitsel.core.ml.pls_ds`."""
        return pls_ds(self._obj, xname, yname, n_components=n_components)

