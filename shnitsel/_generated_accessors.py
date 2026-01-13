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
from shnitsel.analyze.generic import keep_norming, norm, pwdists, subtract_combinations
from shnitsel.analyze.hops import assign_hop_time, focus_hops, hops
from shnitsel.analyze.lda import lda
from shnitsel.analyze.pca import pairwise_dists_pca, pca, pca_and_hops
from shnitsel.analyze.pls import pls, pls_ds
from shnitsel.analyze.populations import calc_classical_populations
from shnitsel.analyze.spectra import assign_fosc, ds_broaden_gauss, spectra_all_times
from shnitsel.analyze.stats import calc_confidence_interval, get_inter_state, get_per_state, time_grouped_confidence_interval
from shnitsel.bridges import default_mol, smiles_map, to_mol, to_xyz, traj_to_xyz
from shnitsel.clean.common import omit, transect, true_upto, truncate
from shnitsel.clean.filter_energy import energy_filtranda, sanity_check
from shnitsel.clean.filter_geo import bond_length_filtranda, filter_by_length
from shnitsel.core.typedefs import DatasetOrArray
from shnitsel.data.helpers import validate
from shnitsel.data.multi_indices import assign_levels, expand_midx, flatten_levels, mdiff, mgroupby, msel, sel_trajids, sel_trajs, stack_trajs, unstack_trajs
from shnitsel.geo.geocalc import angle, dihedral, distance, get_bats, get_bla_chromophor, get_bond_angles, get_bond_lengths, get_bond_torsions, get_pyramids, kabsch
from shnitsel.geo.geomatch_exact import get_bats_matching
from shnitsel.io.ase.write import write_ase_db
from shnitsel.io.shnitsel.write import write_shnitsel_file
from shnitsel.units.conversion import convert_dipole, convert_energy, convert_force, convert_length, convert_nacs, convert_time
from shnitsel.vis.plot.p3mhelpers import frame3D, frames3Dgrid, traj3D, trajs3Dgrid
from shnitsel.vis.plot.select import FrameSelector, TrajSelector
from shnitsel.vis.vmd import traj_vmd
from typing import Callable, Dict, Hashable, List, Literal, Optional, Sequence, Union
from xarray.core.dataarray import DataArray
from xarray.core.dataset import Dataset
from xarray.core.groupby import DataArrayGroupBy, DatasetGroupBy


class DataArrayAccessor(DAManualAccessor):
    _methods = [
        'norm',
        'subtract_combinations',
        'keep_norming',
        'pwdists',
        'calc_confidence_interval',
        'time_grouped_confidence_interval',
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
        'true_upto',
        'dihedral',
        'angle',
        'distance',
        'get_bond_lengths',
        'get_bond_angles',
        'get_bond_torsions',
        'get_pyramids',
        'get_bla_chromophor',
        'get_bats',
        'kabsch',
        'get_bats_matching',
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
        """Wrapper for :py:func:`shnitsel.analyze.generic.norm`."""
        return norm(self._obj, dim=dim, keep_attrs=keep_attrs)

    def subtract_combinations(self, dim: Hashable, labels: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.generic.subtract_combinations`."""
        return subtract_combinations(self._obj, dim, labels=labels)

    def keep_norming(self, exclude: Optional=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.generic.keep_norming`."""
        return keep_norming(self._obj, exclude=exclude)

    def pwdists(self, mean: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.generic.pwdists`."""
        return pwdists(self._obj, mean=mean)

    def calc_confidence_interval(self, confidence: float=0.95) -> ndarray:
        """Wrapper for :py:func:`shnitsel.analyze.stats.calc_confidence_interval`."""
        return calc_confidence_interval(self._obj, confidence=confidence)

    @needs(dims={'frame'}, groupable={'time'})
    def time_grouped_confidence_interval(self, confidence: float=0.9) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.stats.time_grouped_confidence_interval`."""
        return time_grouped_confidence_interval(self._obj, confidence=confidence)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def to_xyz(self, comment='#', units='angstrom') -> str:
        """Wrapper for :py:func:`shnitsel.bridges.to_xyz`."""
        return to_xyz(self._obj, comment=comment, units=units)

    @needs(dims={'atom', 'direction'}, groupable={'time'}, coords_or_vars={'atNames'})
    def traj_to_xyz(self, units='angstrom') -> str:
        """Wrapper for :py:func:`shnitsel.bridges.traj_to_xyz`."""
        return traj_to_xyz(self._obj, units=units)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def to_mol(self, charge: int | None=None, covFactor: float=1.2, to2D: bool=True, molAtomMapNumber: Union=None, atomNote: Union=None, atomLabel: Union=None) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.to_mol`."""
        return to_mol(self._obj, charge=charge, covFactor=covFactor, to2D=to2D, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def smiles_map(self, charge=0, covFactor=1.5) -> str:
        """Wrapper for :py:func:`shnitsel.bridges.smiles_map`."""
        return smiles_map(self._obj, charge=charge, covFactor=covFactor)

    def default_mol(self, molAtomMapNumber: Union=None, atomNote: Union=None, atomLabel: Union=None) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.default_mol`."""
        return default_mol(self._obj, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel)

    @needs(dims={'atom'})
    def pairwise_dists_pca(self, mean: bool=False, return_pca_object=False, **kwargs) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.pca.pairwise_dists_pca`."""
        return pairwise_dists_pca(self._obj, mean=mean, return_pca_object=return_pca_object, **kwargs)

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
        """Wrapper for :py:func:`shnitsel.data.multi_indices.mdiff`."""
        return mdiff(self._obj)

    def flatten_levels(self, idx_name: str, levels: Sequence[str], new_name: str | None=None, position: int=0, renamer: Callable | None=None) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.flatten_levels`."""
        return flatten_levels(self._obj, idx_name, levels, new_name=new_name, position=position, renamer=renamer)

    def expand_midx(self, midx_name: str, level_name: str, value) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.expand_midx`."""
        return expand_midx(self._obj, midx_name, level_name, value)

    def assign_levels(self, levels: dict[str, npt.ArrayLike] | None=None, **levels_kwargs: npt.ArrayLike) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.assign_levels`."""
        return assign_levels(self._obj, levels=levels, **levels_kwargs)

    def mgroupby(self, levels: Sequence[str]) -> DataArrayGroupBy | DatasetGroupBy:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.mgroupby`."""
        return mgroupby(self._obj, levels)

    def msel(self, **kwargs) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.msel`."""
        return msel(self._obj, **kwargs)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert=False) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajids(self, trajids: npt.ArrayLike, invert=False) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.sel_trajids`."""
        return sel_trajids(self._obj, trajids, invert=invert)

    def true_upto(self, dim: Hashable) -> DataArray:
        """Wrapper for :py:func:`shnitsel.clean.common.true_upto`."""
        return true_upto(self._obj, dim)

    @needs(dims={'atom'})
    def dihedral(self, i: int, j: int, k: int, l: int, deg: bool=False, full: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.dihedral`."""
        return dihedral(self._obj, i, j, k, l, deg=deg, full=full)

    @needs(dims={'atom'})
    def angle(self, i: int, j: int, k: int, deg: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.angle`."""
        return angle(self._obj, i, j, k, deg=deg)

    @needs(dims={'atom'})
    def distance(self, i: int, j: int) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.distance`."""
        return distance(self._obj, i, j)

    @needs(dims={'atom', 'direction'})
    def get_bond_lengths(self, matches_or_mol: dict | rdkit.Chem.rdchem.Mol | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_bond_lengths`."""
        return get_bond_lengths(self._obj, matches_or_mol=matches_or_mol)

    @needs(dims={'atom', 'direction'})
    def get_bond_angles(self, matches_or_mol: dict | rdkit.Chem.rdchem.Mol | None=None, mol: rdkit.Chem.rdchem.Mol | None=None, ang: Literal=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_bond_angles`."""
        return get_bond_angles(self._obj, matches_or_mol=matches_or_mol, mol=mol, ang=ang)

    @needs(dims={'atom', 'direction'})
    def get_bond_torsions(self, matches_or_mol: dict | None=None, signed: bool | None=None, ang: Literal=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_bond_torsions`."""
        return get_bond_torsions(self._obj, matches_or_mol=matches_or_mol, signed=signed, ang=ang)

    def get_pyramids(self, pyramid_idxs: dict[int, list[int]] | None=None, mol: rdkit.Chem.rdchem.Mol | None=None, deg: bool=False, signed=True) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_pyramids`."""
        return get_pyramids(self._obj, pyramid_idxs=pyramid_idxs, mol=mol, deg=deg, signed=signed)

    @needs(dims={'atom', 'direction'})
    def get_bla_chromophor(self, matches_or_mol: dict | rdkit.Chem.rdchem.Mol | None=None, mol: rdkit.Chem.rdchem.Mol | None=None, ang: Literal=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_bla_chromophor`."""
        return get_bla_chromophor(self._obj, matches_or_mol=matches_or_mol, mol=mol, ang=ang)

    @needs(dims={'atom', 'direction'})
    def get_bats(self, matches_or_mol: dict | rdkit.Chem.rdchem.Mol | None=None, signed: bool | None=None, ang: Literal=False, pyr=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_bats`."""
        return get_bats(self._obj, matches_or_mol=matches_or_mol, signed=signed, ang=ang, pyr=pyr)

    @needs(dims={'atom', 'direction'})
    def kabsch(self, reference_or_indexers: xarray.core.dataarray.DataArray | dict | None=None, **indexers_kwargs) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc.kabsch`."""
        return kabsch(self._obj, reference_or_indexers=reference_or_indexers, **indexers_kwargs)

    @needs(dims={'atom', 'direction'})
    def get_bats_matching(self, l_smarts: list, signed: bool | None=None, ang: Literal=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geomatch_exact.get_bats_matching`."""
        return get_bats_matching(self._obj, l_smarts, signed=signed, ang=ang)

    def FrameSelector(self, data_var=None, dim=None, xname=None, yname=None, title='', allowed_ws_origin=None, webgl=True):
        """Wrapper for :py:func:`shnitsel.vis.plot.select.FrameSelector`."""
        return FrameSelector(self._obj, data_var=data_var, dim=dim, xname=xname, yname=yname, title=title, allowed_ws_origin=allowed_ws_origin, webgl=webgl)

    def TrajSelector(self, data_var=None, dim=None, xname=None, yname=None, title='', allowed_ws_origin=None, webgl=True):
        """Wrapper for :py:func:`shnitsel.vis.plot.select.TrajSelector`."""
        return TrajSelector(self._obj, data_var=data_var, dim=dim, xname=xname, yname=yname, title=title, allowed_ws_origin=allowed_ws_origin, webgl=webgl)

    @needs(dims={'atom', 'direction'}, coords_or_vars={'atNames'}, not_dims={'frame'})
    def frame3D(self):
        """Wrapper for :py:func:`shnitsel.vis.plot.p3mhelpers.frame3D`."""
        return frame3D(self._obj)

    @needs(dims={'atom', 'direction'}, groupable={'frame'}, coords_or_vars={'atNames'})
    def frames3Dgrid(self):
        """Wrapper for :py:func:`shnitsel.vis.plot.p3mhelpers.frames3Dgrid`."""
        return frames3Dgrid(self._obj)

    @needs(dims={'atom', 'direction'}, groupable={'time'}, coords_or_vars={'atNames'})
    def traj3D(self):
        """Wrapper for :py:func:`shnitsel.vis.plot.p3mhelpers.traj3D`."""
        return traj3D(self._obj)

    @needs(dims={'atom', 'direction'}, coords={'trajid'}, groupable={'time'}, coords_or_vars={'atNames'})
    def trajs3Dgrid(self, trajids: list[int | str] | None=None, loop: str='forward'):
        """Wrapper for :py:func:`shnitsel.vis.plot.p3mhelpers.trajs3Dgrid`."""
        return trajs3Dgrid(self._obj, trajids=trajids, loop=loop)

    def traj_vmd(self, groupby='trajid'):
        """Wrapper for :py:func:`shnitsel.vis.vmd.traj_vmd`."""
        return traj_vmd(self._obj, groupby=groupby)

    def pca(self, dim: str, n_components: int=2, return_pca_object: bool=False) -> tuple[xarray.core.dataarray.DataArray, sklearn.decomposition._pca.PCA] | xarray.core.dataarray.DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.pca.pca`."""
        return pca(self._obj, dim, n_components=n_components, return_pca_object=return_pca_object)

    def lda(self, dim: str, cats: str | xarray.core.dataarray.DataArray, n_components: int=2) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.lda.lda`."""
        return lda(self._obj, dim, cats, n_components=n_components)

    def pls(self, ydata_array: DataArray, n_components: int=2, common_dim: str | None=None) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.pls.pls`."""
        return pls(self._obj, ydata_array, n_components=n_components, common_dim=common_dim)


class DatasetAccessor(DSManualAccessor):
    _methods = [
        'pca_and_hops',
        'validate',
        'assign_fosc',
        'ds_broaden_gauss',
        'get_per_state',
        'get_inter_state',
        'calc_classical_populations',
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
        'hops',
        'focus_hops',
        'assign_hop_time',
        'FrameSelector',
        'TrajSelector',
    ]

    @needs(coords_or_vars={'astate', 'atXYZ'})
    def pca_and_hops(self, mean: bool) -> tuple:
        """Wrapper for :py:func:`shnitsel.analyze.pca.pca_and_hops`."""
        return pca_and_hops(self._obj, mean)

    def validate(self) -> ndarray:
        """Wrapper for :py:func:`shnitsel.data.helpers.validate`."""
        return validate(self._obj)

    @needs(coords={'statecomb'}, data_vars={'dip_trans', 'energy_interstate'})
    def assign_fosc(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.spectra.assign_fosc`."""
        return assign_fosc(self._obj)

    @needs(data_vars={'energy_interstate', 'fosc'})
    def ds_broaden_gauss(self, width_in_eV: float=0.5, nsamples: int=1000, xmax: float | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.spectra.ds_broaden_gauss`."""
        return ds_broaden_gauss(self._obj, width_in_eV=width_in_eV, nsamples=nsamples, xmax=xmax)

    @needs(dims={'state'})
    def get_per_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.stats.get_per_state`."""
        return get_per_state(self._obj)

    @needs(dims={'state'}, coords={'state'})
    def get_inter_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.stats.get_inter_state`."""
        return get_inter_state(self._obj)

    @needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
    def calc_classical_populations(self) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.populations.calc_classical_populations`."""
        return calc_classical_populations(self._obj)

    def default_mol(self, molAtomMapNumber: Union=None, atomNote: Union=None, atomLabel: Union=None) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.default_mol`."""
        return default_mol(self._obj, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel)

    def flatten_levels(self, idx_name: str, levels: Sequence[str], new_name: str | None=None, position: int=0, renamer: Callable | None=None) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.flatten_levels`."""
        return flatten_levels(self._obj, idx_name, levels, new_name=new_name, position=position, renamer=renamer)

    def expand_midx(self, midx_name: str, level_name: str, value) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.expand_midx`."""
        return expand_midx(self._obj, midx_name, level_name, value)

    def assign_levels(self, levels: dict[str, npt.ArrayLike] | None=None, **levels_kwargs: npt.ArrayLike) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.assign_levels`."""
        return assign_levels(self._obj, levels=levels, **levels_kwargs)

    def mgroupby(self, levels: Sequence[str]) -> DataArrayGroupBy | DatasetGroupBy:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.mgroupby`."""
        return mgroupby(self._obj, levels)

    def msel(self, **kwargs) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.msel`."""
        return msel(self._obj, **kwargs)

    @needs(dims={'frame'}, coords_or_vars={'trajid'})
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert=False) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    def unstack_trajs(self) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.unstack_trajs`."""
        return unstack_trajs(self._obj)

    def stack_trajs(self) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.stack_trajs`."""
        return stack_trajs(self._obj)

    def write_shnitsel_file(self, savepath: str | os.PathLike, complevel: int=9):
        """Wrapper for :py:func:`shnitsel.io.shnitsel.write.write_shnitsel_file`."""
        return write_shnitsel_file(self._obj, savepath, complevel=complevel)

    @needs(coords={'frame', 'trajid'}, data_vars={'energy_interstate', 'fosc'})
    def spectra_all_times(self) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.spectra.spectra_all_times`."""
        return spectra_all_times(self._obj)

    def energy_filtranda(self, etot_drift: float | None=None, etot_step: float | None=None, epot_step: float | None=None, ekin_step: float | None=None, hop_epot_step: float | None=None, units='eV'):
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.energy_filtranda`."""
        return energy_filtranda(self._obj, etot_drift=etot_drift, etot_step=etot_step, epot_step=epot_step, ekin_step=ekin_step, hop_epot_step=hop_epot_step, units=units)

    def sanity_check(self, cut: Union='truncate', units='eV', etot_drift: float=nan, etot_step: float=nan, epot_step: float=nan, ekin_step: float=nan, hop_epot_step: float=nan, plot_thresholds: Union=False, plot_populations: Union=False):
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.sanity_check`."""
        return sanity_check(self._obj, cut=cut, units=units, etot_drift=etot_drift, etot_step=etot_step, epot_step=epot_step, ekin_step=ekin_step, hop_epot_step=hop_epot_step, plot_thresholds=plot_thresholds, plot_populations=plot_populations)

    def bond_length_filtranda(self, search_dict: dict[str, numbers.Number] | None=None, units='angstrom', mol: rdkit.Chem.rdchem.Mol | None=None):
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.bond_length_filtranda`."""
        return bond_length_filtranda(self._obj, search_dict=search_dict, units=units, mol=mol)

    def filter_by_length(self, cut: Union='truncate', search_dict: dict[str, numbers.Number] | None=None, units: str='angstrom', plot_thresholds: Union=False, plot_populations: Union=False, mol: rdkit.Chem.rdchem.Mol | None=None):
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.filter_by_length`."""
        return filter_by_length(self._obj, cut=cut, search_dict=search_dict, units=units, plot_thresholds=plot_thresholds, plot_populations=plot_populations, mol=mol)

    @needs(data_vars={'filtranda'})
    def omit(self):
        """Wrapper for :py:func:`shnitsel.clean.common.omit`."""
        return omit(self._obj)

    @needs(data_vars={'filtranda'})
    def truncate(self):
        """Wrapper for :py:func:`shnitsel.clean.common.truncate`."""
        return truncate(self._obj)

    @needs(data_vars={'filtranda'})
    def transect(self, cutoff: float):
        """Wrapper for :py:func:`shnitsel.clean.common.transect`."""
        return transect(self._obj, cutoff)

    @needs(data_vars={'atNames', 'atNums', 'atXYZ', 'energy'})
    def write_ase_db(self, db_path: str, db_format: Optional, keys_to_write: Optional=None, preprocess: bool=True):
        """Wrapper for :py:func:`shnitsel.io.ase.write.write_ase_db`."""
        return write_ase_db(self._obj, db_path, db_format, keys_to_write=keys_to_write, preprocess=preprocess)

    def pls_ds(self, xname: str, yname: str, n_components: int=2, common_dim: str | None=None) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.pls.pls_ds`."""
        return pls_ds(self._obj, xname, yname, n_components=n_components, common_dim=common_dim)

    def hops(self, hop_types: list[tuple[int, int]] | None=None):
        """Wrapper for :py:func:`shnitsel.analyze.hops.hops`."""
        return hops(self._obj, hop_types=hop_types)

    def focus_hops(self, hop_types: list[tuple[int, int]] | None=None, window: slice | None=None):
        """Wrapper for :py:func:`shnitsel.analyze.hops.focus_hops`."""
        return focus_hops(self._obj, hop_types=hop_types, window=window)

    def assign_hop_time(self, hop_types: list[tuple[int, int]] | None=None, which: Literal='last'):
        """Wrapper for :py:func:`shnitsel.analyze.hops.assign_hop_time`."""
        return assign_hop_time(self._obj, hop_types=hop_types, which=which)

    def FrameSelector(self, data_var=None, dim=None, xname=None, yname=None, title='', allowed_ws_origin=None, webgl=True):
        """Wrapper for :py:func:`shnitsel.vis.plot.select.FrameSelector`."""
        return FrameSelector(self._obj, data_var=data_var, dim=dim, xname=xname, yname=yname, title=title, allowed_ws_origin=allowed_ws_origin, webgl=webgl)

    def TrajSelector(self, data_var=None, dim=None, xname=None, yname=None, title='', allowed_ws_origin=None, webgl=True):
        """Wrapper for :py:func:`shnitsel.vis.plot.select.TrajSelector`."""
        return TrajSelector(self._obj, data_var=data_var, dim=dim, xname=xname, yname=yname, title=title, allowed_ws_origin=allowed_ws_origin, webgl=webgl)

