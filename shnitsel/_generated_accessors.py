import collections
import numbers
import numpy
import numpy.typing as npt
import os
import pathlib
import rdkit
import sklearn
import typing
import xarray
import xarray as xr
from ._accessors import DAManualAccessor, DSManualAccessor
from ._contracts import needs
from numpy import nan, ndarray
from rdkit.Chem.rdchem import Mol
from shnitsel.bridges import default_mol, smiles_map, to_mol, to_xyz, traj_to_xyz
from shnitsel.clean.common import omit, transect, truncate
from shnitsel.clean.filter_energy import energy_filtranda, sanity_check
from shnitsel.clean.filter_geo import bond_length_filtranda, filter_by_length
from shnitsel.core.convenience import pairwise_dists_pca, pca_and_hops, validate
from shnitsel.core.generic import keep_norming, norm, subtract_combinations
from shnitsel.core.geom import angle, dihedral, distance, get_bats, get_bond_angles, get_bond_lengths, get_bond_torsions, get_pyramids, kabsch
from shnitsel.core.midx import assign_levels, expand_midx, flatten_levels, mdiff, mgroupby, msel, sel_trajids, sel_trajs, stack_trajs, unstack_trajs
from shnitsel.core.ml import lda, pca, pls, pls_ds
from shnitsel.core.plot.p3mhelpers import frame3D, frames3Dgrid, traj3D, trajs3Dgrid
from shnitsel.core.plot.select import FrameSelector, TrajSelector
from shnitsel.core.populations import classical
from shnitsel.core.spectra import assign_fosc, ds_broaden_gauss, spectra_all_times
from shnitsel.core.stats import calc_ci, get_inter_state, get_per_state, time_grouped_ci
from shnitsel.core.vmd import traj_vmd
from shnitsel.io.ase.write import write_ase_db
from shnitsel.io.sharc.parse_initial_conditions import iconds_to_frames
from shnitsel.io.shnitsel.write import write_shnitsel_file
from shnitsel.units.conversion import convert_dipole, convert_energy, convert_force, convert_length, convert_nacs, convert_time
from typing import Callable, Dict, Hashable, List, Literal, Optional, Sequence, Union
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.groupby import DataArrayGroupBy, DatasetGroupBy


class DataArrayAccessor(DAManualAccessor):
    _methods = [
        'norm',
        'subtract_combinations',
        'keep_norming',
        'calc_ci',
        'time_grouped_ci',
        'to_xyz',
        'traj_to_xyz',
        'to_mol',
        'smiles_map',
        'default_mol',
        'pairwise_dists_pca',
        'convert_energy',
        'convert_force',
        'convert_dipole',
        'convert_length',
        'convert_time',
        'convert_nacs',
        'convert_time',
        'mdiff',
        'flatten_levels',
        'expand_midx',
        'assign_levels',
        'mgroupby',
        'msel',
        'sel_trajs',
        'sel_trajids',
        'dihedral',
        'angle',
        'distance',
        'get_bond_lengths',
        'get_bond_angles',
        'get_bond_torsions',
        'get_pyramids',
        'get_bats',
        'kabsch',
        'FrameSelector',
        'TrajSelector',
        'frame3D',
        'frames3Dgrid',
        'traj3D',
        'trajs3Dgrid',
        'traj_vmd',
        'pca',
        'lda',
        'pls',
    ]

    def norm(self, dim: Hashable='direction', keep_attrs: bool | str | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.generic.norm`."""
        return norm(self._obj, dim=dim, keep_attrs=keep_attrs)

    def subtract_combinations(self, dim: Hashable, labels: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.generic.subtract_combinations`."""
        return subtract_combinations(self._obj, dim, labels=labels)

    def keep_norming(self, exclude: Optional=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.generic.keep_norming`."""
        return keep_norming(self._obj, exclude=exclude)

    def calc_ci(self, confidence: float=0.95) -> ndarray:
        """Wrapper for :py:func:`shnitsel.core.stats.calc_ci`."""
        return calc_ci(self._obj, confidence=confidence)

    @needs(dims={'frame'}, groupable={'time'})
    def time_grouped_ci(self, confidence: float=0.9) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.stats.time_grouped_ci`."""
        return time_grouped_ci(self._obj, confidence=confidence)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def to_xyz(self, comment='#') -> str:
        """Wrapper for :py:func:`shnitsel.bridges.to_xyz`."""
        return to_xyz(self._obj, comment=comment)

    @needs(dims={'atom', 'direction'}, groupable={'time'}, coords_or_vars={'atNames'})
    def traj_to_xyz(self) -> str:
        """Wrapper for :py:func:`shnitsel.bridges.traj_to_xyz`."""
        return traj_to_xyz(self._obj)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def to_mol(self, charge: int | None=None, covFactor: float=1.2, to2D: bool=True, molAtomMapNumber: Union=None, atomNote: Union=None, atomLabel: Union=None) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.to_mol`."""
        return to_mol(self._obj, charge=charge, covFactor=covFactor, to2D=to2D, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def smiles_map(self, charge=0, covFactor=1.5) -> str:
        """Wrapper for :py:func:`shnitsel.bridges.smiles_map`."""
        return smiles_map(self._obj, charge=charge, covFactor=covFactor)

    def default_mol(self) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.default_mol`."""
        return default_mol(self._obj)

    @needs(dims={'atom'})
    def pairwise_dists_pca(self, **kwargs) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.convenience.pairwise_dists_pca`."""
        return pairwise_dists_pca(self._obj, **kwargs)

    def convert_energy(self, to: str, convert_from: str | None=None):
        """Wrapper for :py:func:`shnitsel.units.conversion.convert_energy`."""
        return convert_energy(self._obj, to, convert_from=convert_from)

    def convert_force(self, to: str, convert_from: str | None=None):
        """Wrapper for :py:func:`shnitsel.units.conversion.convert_force`."""
        return convert_force(self._obj, to, convert_from=convert_from)

    def convert_dipole(self, to: str, convert_from: str | None=None):
        """Wrapper for :py:func:`shnitsel.units.conversion.convert_dipole`."""
        return convert_dipole(self._obj, to, convert_from=convert_from)

    def convert_length(self, to: str, convert_from: str | None=None):
        """Wrapper for :py:func:`shnitsel.units.conversion.convert_length`."""
        return convert_length(self._obj, to, convert_from=convert_from)

    def convert_time(self, to: str, convert_from: str | None=None):
        """Wrapper for :py:func:`shnitsel.units.conversion.convert_time`."""
        return convert_time(self._obj, to, convert_from=convert_from)

    def convert_nacs(self, to: str, convert_from: str | None=None):
        """Wrapper for :py:func:`shnitsel.units.conversion.convert_nacs`."""
        return convert_nacs(self._obj, to, convert_from=convert_from)

    def convert_time(self, to: str, convert_from: str | None=None):
        """Wrapper for :py:func:`shnitsel.units.conversion.convert_time`."""
        return convert_time(self._obj, to, convert_from=convert_from)

    @needs(dims={'frame'})
    def mdiff(self) -> xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.mdiff`."""
        return mdiff(self._obj)

    def flatten_levels(self, idx_name: str, levels: Sequence[str], new_name: str | None=None, position: int=0, renamer: Callable | None=None) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.flatten_levels`."""
        return flatten_levels(self._obj, idx_name, levels, new_name=new_name, position=position, renamer=renamer)

    def expand_midx(self, midx_name: str, level_name: str, value) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.expand_midx`."""
        return expand_midx(self._obj, midx_name, level_name, value)

    def assign_levels(self, levels: dict[str, npt.ArrayLike] | None=None, **levels_kwargs: npt.ArrayLike) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.assign_levels`."""
        return assign_levels(self._obj, levels=levels, **levels_kwargs)

    def mgroupby(self, levels: Sequence[str]) -> DataArrayGroupBy | DatasetGroupBy:
        """Wrapper for :py:func:`shnitsel.core.midx.mgroupby`."""
        return mgroupby(self._obj, levels)

    def msel(self, **kwargs) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.msel`."""
        return msel(self._obj, **kwargs)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert=False) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajids(self, trajids: npt.ArrayLike, invert=False) -> xr.Dataset:
        """Wrapper for :py:func:`shnitsel.core.midx.sel_trajids`."""
        return sel_trajids(self._obj, trajids, invert=invert)

    @needs(dims={'atom'})
    def dihedral(self, i: int, j: int, k: int, l: int, deg: bool=False, full: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.geom.dihedral`."""
        return dihedral(self._obj, i, j, k, l, deg=deg, full=full)

    @needs(dims={'atom'})
    def angle(self, i: int, j: int, k: int, deg: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.geom.angle`."""
        return angle(self._obj, i, j, k, deg=deg)

    @needs(dims={'atom'})
    def distance(self, i: int, j: int) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.geom.distance`."""
        return distance(self._obj, i, j)

    @needs(dims={'atom', 'direction'})
    def get_bond_lengths(self, matches_or_mol: dict | rdkit.Chem.rdchem.Mol | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.geom.get_bond_lengths`."""
        return get_bond_lengths(self._obj, matches_or_mol=matches_or_mol)

    @needs(dims={'atom', 'direction'})
    def get_bond_angles(self, matches_or_mol: dict | rdkit.Chem.rdchem.Mol | None=None, mol: rdkit.Chem.rdchem.Mol | None=None, ang: Literal=False):
        """Wrapper for :py:func:`shnitsel.core.geom.get_bond_angles`."""
        return get_bond_angles(self._obj, matches_or_mol=matches_or_mol, mol=mol, ang=ang)

    @needs(dims={'atom', 'direction'})
    def get_bond_torsions(self, matches_or_mol: dict | None=None, signed: bool | None=None, ang: Literal=False):
        """Wrapper for :py:func:`shnitsel.core.geom.get_bond_torsions`."""
        return get_bond_torsions(self._obj, matches_or_mol=matches_or_mol, signed=signed, ang=ang)

    def get_pyramids(self, pyramid_idxs: dict[int, list[int]] | None=None, mol: rdkit.Chem.rdchem.Mol | None=None, deg: bool=False, signed=True) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.geom.get_pyramids`."""
        return get_pyramids(self._obj, pyramid_idxs=pyramid_idxs, mol=mol, deg=deg, signed=signed)

    @needs(dims={'atom', 'direction'})
    def get_bats(self, matches_or_mol: dict | rdkit.Chem.rdchem.Mol | None=None, signed: bool | None=None, ang: Literal=False, pyr=False):
        """Wrapper for :py:func:`shnitsel.core.geom.get_bats`."""
        return get_bats(self._obj, matches_or_mol=matches_or_mol, signed=signed, ang=ang, pyr=pyr)

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

    def traj_vmd(self, groupby='trajid', scale=0.5):
        """Wrapper for :py:func:`shnitsel.core.vmd.traj_vmd`."""
        return traj_vmd(self._obj, groupby=groupby, scale=scale)

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
        'assign_fosc',
        'ds_broaden_gauss',
        'get_per_state',
        'get_inter_state',
        'classical',
        'default_mol',
        'flatten_levels',
        'expand_midx',
        'assign_levels',
        'mgroupby',
        'msel',
        'sel_trajs',
        'unstack_trajs',
        'stack_trajs',
        'write_shnitsel_file',
        'iconds_to_frames',
        'spectra_all_times',
        'energy_filtranda',
        'sanity_check',
        'bond_length_filtranda',
        'filter_by_length',
        'omit',
        'truncate',
        'transect',
        'write_ase_db',
        'pls_ds',
    ]

    @needs(coords_or_vars={'astate', 'atXYZ'})
    def pca_and_hops(self) -> tuple:
        """Wrapper for :py:func:`shnitsel.core.convenience.pca_and_hops`."""
        return pca_and_hops(self._obj)

    def validate(self) -> ndarray:
        """Wrapper for :py:func:`shnitsel.core.convenience.validate`."""
        return validate(self._obj)

    @needs(data_vars={'dip_trans', 'energy'})
    def assign_fosc(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.spectra.assign_fosc`."""
        return assign_fosc(self._obj)

    def ds_broaden_gauss(self, width: float=0.5, nsamples: int=1000, xmax: float | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.spectra.ds_broaden_gauss`."""
        return ds_broaden_gauss(self._obj, width=width, nsamples=nsamples, xmax=xmax)

    @needs(dims={'state'})
    def get_per_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.stats.get_per_state`."""
        return get_per_state(self._obj)

    def get_inter_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.core.stats.get_inter_state`."""
        return get_inter_state(self._obj)

    @needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
    def classical(self) -> DataArray:
        """Wrapper for :py:func:`shnitsel.core.populations.classical`."""
        return classical(self._obj)

    def default_mol(self) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.default_mol`."""
        return default_mol(self._obj)

    def flatten_levels(self, idx_name: str, levels: Sequence[str], new_name: str | None=None, position: int=0, renamer: Callable | None=None) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.flatten_levels`."""
        return flatten_levels(self._obj, idx_name, levels, new_name=new_name, position=position, renamer=renamer)

    def expand_midx(self, midx_name: str, level_name: str, value) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.expand_midx`."""
        return expand_midx(self._obj, midx_name, level_name, value)

    def assign_levels(self, levels: dict[str, npt.ArrayLike] | None=None, **levels_kwargs: npt.ArrayLike) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.assign_levels`."""
        return assign_levels(self._obj, levels=levels, **levels_kwargs)

    def mgroupby(self, levels: Sequence[str]) -> DataArrayGroupBy | DatasetGroupBy:
        """Wrapper for :py:func:`shnitsel.core.midx.mgroupby`."""
        return mgroupby(self._obj, levels)

    def msel(self, **kwargs) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.msel`."""
        return msel(self._obj, **kwargs)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert=False) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    def unstack_trajs(self) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.unstack_trajs`."""
        return unstack_trajs(self._obj)

    def stack_trajs(self) -> xr.Dataset | xr.DataArray:
        """Wrapper for :py:func:`shnitsel.core.midx.stack_trajs`."""
        return stack_trajs(self._obj)

    def write_shnitsel_file(self, savepath: str | os.PathLike | pathlib.Path, complevel: int=9):
        """Wrapper for :py:func:`shnitsel.io.shnitsel.write.write_shnitsel_file`."""
        return write_shnitsel_file(self._obj, savepath, complevel=complevel)

    @needs(dims={'icond'}, coords={'icond'}, not_dims={'time'})
    def iconds_to_frames(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.io.sharc.parse_initial_conditions.iconds_to_frames`."""
        return iconds_to_frames(self._obj)

    @needs(coords={'frame', 'trajid'}, data_vars={'energy', 'fosc'})
    def spectra_all_times(self):
        """Wrapper for :py:func:`shnitsel.core.spectra.spectra_all_times`."""
        return spectra_all_times(self._obj)

    def energy_filtranda(self, etot_drift: float | None=None, etot_step: float | None=None, epot_step: float | None=None, ekin_step: float | None=None, hop_epot_step: float | None=None, units='eV'):
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.energy_filtranda`."""
        return energy_filtranda(self._obj, etot_drift=etot_drift, etot_step=etot_step, epot_step=epot_step, ekin_step=ekin_step, hop_epot_step=hop_epot_step, units=units)

    def sanity_check(self, cut: Union='truncate', units='eV', etot_drift: float=nan, etot_step: float=nan, epot_step: float=nan, ekin_step: float=nan, hop_epot_step: float=nan):
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.sanity_check`."""
        return sanity_check(self._obj, cut=cut, units=units, etot_drift=etot_drift, etot_step=etot_step, epot_step=epot_step, ekin_step=ekin_step, hop_epot_step=hop_epot_step)

    def bond_length_filtranda(self, search_dict):
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.bond_length_filtranda`."""
        return bond_length_filtranda(self._obj, search_dict)

    def filter_by_length(self, cut: Union='truncate', search_dict: dict[str, numbers.Number] | None=None):
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.filter_by_length`."""
        return filter_by_length(self._obj, cut=cut, search_dict=search_dict)

    def omit(self):
        """Wrapper for :py:func:`shnitsel.clean.common.omit`."""
        return omit(self._obj)

    def truncate(self):
        """Wrapper for :py:func:`shnitsel.clean.common.truncate`."""
        return truncate(self._obj)

    def transect(self, cutoff: float):
        """Wrapper for :py:func:`shnitsel.clean.common.transect`."""
        return transect(self._obj, cutoff)

    @needs(data_vars={'atNames', 'atNums', 'atXYZ', 'energy'})
    def write_ase_db(self, db_path: str, db_format: Optional, keys_to_write: Optional=None, preprocess: bool=True):
        """Wrapper for :py:func:`shnitsel.io.ase.write.write_ase_db`."""
        return write_ase_db(self._obj, db_path, db_format, keys_to_write=keys_to_write, preprocess=preprocess)

    def pls_ds(self, xname, yname, n_components=2):
        """Wrapper for :py:func:`shnitsel.core.ml.pls_ds`."""
        return pls_ds(self._obj, xname, yname, n_components=n_components)

