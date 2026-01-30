import collections
import numbers
import numpy
import numpy.typing as npt
import os
import pathlib
import rdkit
import shnitsel
import sklearn
import typing
import xarray
import xarray as xr
from ._accessors import DAManualAccessor, DSManualAccessor
from ._contracts import needs
from numpy import nan, ndarray
from rdkit.Chem.rdchem import Mol
from shnitsel.analyze.generic import keep_norming, norm, pwdists, subtract_combinations
from shnitsel.analyze.hops import assign_hop_time, filter_data_at_hops, focus_hops, hops_mask_from_active_state
from shnitsel.analyze.lda import lda
from shnitsel.analyze.pca import PCAResult, pca, pca_and_hops
from shnitsel.analyze.pls import pls, pls_ds
from shnitsel.analyze.populations import PopulationStatistics, calc_classical_populations
from shnitsel.analyze.spectra import get_spectra
from shnitsel.analyze.stats import calc_confidence_interval, get_inter_state, get_per_state, time_grouped_confidence_interval
from shnitsel.bridges import construct_default_mol, default_mol, smiles_map, to_mol, to_xyz, traj_to_xyz
from shnitsel.clean import sanity_check
from shnitsel.clean.common import TrajectoryOrFrames, omit, transect, true_upto, truncate
from shnitsel.clean.filter_energy import calculate_energy_filtranda, filter_by_energy
from shnitsel.clean.filter_geo import calculate_bond_length_filtranda, filter_by_length
from shnitsel.core.typedefs import DataArrayOrVar, DatasetOrArray
from shnitsel.data.helpers import validate
from shnitsel.data.multi_indices import assign_levels, expand_midx, flatten_levels, mdiff, mgroupby, msel, sel_trajs, stack_trajs, unstack_trajs
from shnitsel.geo.alignment import kabsch
from shnitsel.geo.geocalc import get_bats
from shnitsel.geo.geocalc_.angles import angle, get_angles
from shnitsel.geo.geocalc_.bla_chromophor import get_max_chromophor_BLA
from shnitsel.geo.geocalc_.dihedrals import dihedral, get_dihedrals
from shnitsel.geo.geocalc_.distances import distance, get_distances
from shnitsel.geo.geocalc_.pyramids import get_pyramidalization, pyramidalization_angle
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

default_mol = construct_default_mol
default_mol.__name__ = 'default_mol'


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
        'convert_energy',
        'convert_force',
        'convert_dipole',
        'convert_length',
        'convert_time',
        'convert_nacs',
        'mdiff',
        'flatten_levels',
        'expand_midx',
        'assign_levels',
        'mgroupby',
        'msel',
        'sel_trajs',
        'stack_trajs',
        'unstack_trajs',
        'true_upto',
        'distance',
        'angle',
        'dihedral',
        'pyramidalization_angle',
        'get_bats',
        'get_distances',
        'get_angles',
        'get_dihedrals',
        'get_pyramidalization',
        'get_max_chromophor_BLA',
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
        'hops_mask_from_active_state',
        'filter_data_at_hops',
        'focus_hops',
        'assign_hop_time',
    ]

    def norm(self, dim: str='direction', keep_attrs: bool | str | None=None) -> DataArrayOrVar:
        """Wrapper for :py:func:`shnitsel.analyze.generic.norm`."""
        return norm(self._obj, dim=dim, keep_attrs=keep_attrs)

    def subtract_combinations(self, dim: str, add_labels: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.generic.subtract_combinations`."""
        return subtract_combinations(self._obj, dim, add_labels=add_labels)

    def keep_norming(self, exclude: Optional=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.generic.keep_norming`."""
        return keep_norming(self._obj, exclude=exclude)

    def pwdists(self, center_mean: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.generic.pwdists`."""
        return pwdists(self._obj, center_mean=center_mean)

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

    def default_mol(self, to2D: bool=True, charge: int | float | None=None, molAtomMapNumber: Union=None, atomNote: Union=None, atomLabel: Union=None, silent_mode: bool=False) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.default_mol`."""
        return default_mol(self._obj, to2D=to2D, charge=charge, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel, silent_mode=silent_mode)

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

    @needs(dims={'frame'})
    def mdiff(self, dim: str | None=None) -> xr.DataArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.mdiff`."""
        return mdiff(self._obj, dim=dim)

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
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert: bool=False) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    def stack_trajs(self) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.stack_trajs`."""
        return stack_trajs(self._obj)

    def unstack_trajs(self, fill_value=shnitsel.data.multi_indices.dtype_NA) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.unstack_trajs`."""
        return unstack_trajs(self._obj, fill_value=fill_value)

    def true_upto(self, dim: str) -> DataArray:
        """Wrapper for :py:func:`shnitsel.clean.common.true_upto`."""
        return true_upto(self._obj, dim)

    @needs(dims={'atom'})
    def distance(self, i: int, j: int) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.distances.distance`."""
        return distance(self._obj, i, j)

    @needs(dims={'atom'})
    def angle(self, a_index: int, b_index: int, c_index: int, deg: bool=False) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.angles.angle`."""
        return angle(self._obj, a_index, b_index, c_index, deg=deg)

    @needs(dims={'atom'})
    def dihedral(self, a_index: int, b_index: int, c_index: int, d_index: int, deg: Union=True, full: bool=False) -> "xarray.core.dataarray.DataArray | tuple[xarray.core.dataarray.DataArray, xarray.core.dataarray.DataArray]":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.dihedrals.dihedral`."""
        return dihedral(self._obj, a_index, b_index, c_index, d_index, deg=deg, full=full)

    @needs(dims={'atom', 'direction'})
    def pyramidalization_angle(self, x_index: int, a_index: int, b_index: int, c_index: int, deg: Union=True) -> "xarray.core.dataarray.DataArray | tuple[xarray.core.dataarray.DataArray, xarray.core.dataarray.DataArray]":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.pyramids.pyramidalization_angle`."""
        return pyramidalization_angle(self._obj, x_index, a_index, b_index, c_index, deg=deg)

    @needs(dims={'atom', 'direction'})
    def get_bats(self, structure_selection: Union=None, default_features: Sequence=['bonds', 'angles', 'dihedrals'], signed: bool=False, deg: Union=True) -> "xarray.core.dataarray.DataArray | shnitsel.data.tree.node.TreeNode[DataArray]":
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_bats`."""
        return get_bats(self._obj, structure_selection=structure_selection, default_features=default_features, signed=signed, deg=deg)

    @needs(dims={'atom', 'direction'})
    def get_distances(self, structure_selection: Union=None) -> "shnitsel.data.tree.node.TreeNode[DataArray] | xarray.core.dataarray.DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.distances.get_distances`."""
        return get_distances(self._obj, structure_selection=structure_selection)

    @needs(dims={'atom', 'direction'})
    def get_angles(self, structure_selection: Union=None, deg: Union=True, signed: bool=True) -> "shnitsel.data.tree.node.TreeNode[DataArray] | xarray.core.dataarray.DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.angles.get_angles`."""
        return get_angles(self._obj, structure_selection=structure_selection, deg=deg, signed=signed)

    @needs(dims={'atom', 'direction'})
    def get_dihedrals(self, structure_selection: Union=None, deg: Union=True, signed: bool=True) -> "shnitsel.data.tree.node.TreeNode[DataArray] | xarray.core.dataarray.DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.dihedrals.get_dihedrals`."""
        return get_dihedrals(self._obj, structure_selection=structure_selection, deg=deg, signed=signed)

    def get_pyramidalization(self, structure_selection: Union=None, deg: Union=True, signed: bool=True) -> "shnitsel.data.tree.node.TreeNode[DataArray] | xarray.core.dataarray.DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.pyramids.get_pyramidalization`."""
        return get_pyramidalization(self._obj, structure_selection=structure_selection, deg=deg, signed=signed)

    @needs(dims={'atom', 'direction'})
    def get_max_chromophor_BLA(self, structure_selection: Union=None, SMARTS: str | None=None, num_double_bonds: int | None=None, allowed_chain_elements: str='#6,#7,#8,#15,#16', max_considered_BLA_double_bonds: int=50) -> "shnitsel.data.tree.node.TreeNode[DataArray] | xarray.core.dataarray.DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.bla_chromophor.get_max_chromophor_BLA`."""
        return get_max_chromophor_BLA(self._obj, structure_selection=structure_selection, SMARTS=SMARTS, num_double_bonds=num_double_bonds, allowed_chain_elements=allowed_chain_elements, max_considered_BLA_double_bonds=max_considered_BLA_double_bonds)

    @needs(dims={'atom', 'direction'})
    def kabsch(self, reference_or_indexers: xarray.core.dataarray.DataArray | dict | None=None, **indexers_kwargs) -> DataArray:
        """Wrapper for :py:func:`shnitsel.geo.alignment.kabsch`."""
        return kabsch(self._obj, reference_or_indexers=reference_or_indexers, **indexers_kwargs)

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

    def pca(self, structure_selection: Union=None, dim: Optional=None, n_components: int=2, center_mean: bool=False) -> Union:
        """Wrapper for :py:func:`shnitsel.analyze.pca.pca`."""
        return pca(self._obj, structure_selection=structure_selection, dim=dim, n_components=n_components, center_mean=center_mean)

    def lda(self, dim: str, cats: str | xarray.core.dataarray.DataArray, n_components: int=2) -> DataArray:
        """Wrapper for :py:func:`shnitsel.analyze.lda.lda`."""
        return lda(self._obj, dim, cats, n_components=n_components)

    def pls(self, ydata_array: DataArray, n_components: int=2, common_dim: str | None=None) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.pls.pls`."""
        return pls(self._obj, ydata_array, n_components=n_components, common_dim=common_dim)

    def hops_mask_from_active_state(self, hop_type_selection: Union=None, dim: str | None=None) -> "xarray.core.dataarray.DataArray | shnitsel.data.tree.node.TreeNode[DataArray]":
        """Wrapper for :py:func:`shnitsel.analyze.hops.hops_mask_from_active_state`."""
        return hops_mask_from_active_state(self._obj, hop_type_selection=hop_type_selection, dim=dim)

    def filter_data_at_hops(self, hop_type_selection: Union=None) -> "shnitsel.data.dataset_containers.data_series.DataSeries | xarray.core.dataarray.DataArray | shnitsel.data.tree.node.TreeNode[DataSeries] | shnitsel.data.tree.node.TreeNode[DataArray]":
        """Wrapper for :py:func:`shnitsel.analyze.hops.filter_data_at_hops`."""
        return filter_data_at_hops(self._obj, hop_type_selection=hop_type_selection)

    def focus_hops(self, hop_types: list[tuple[int, int]] | None=None, window: slice | None=None):
        """Wrapper for :py:func:`shnitsel.analyze.hops.focus_hops`."""
        return focus_hops(self._obj, hop_types=hop_types, window=window)

    def assign_hop_time(self, hop_types: list[tuple[int, int]] | None=None, which: Literal='last'):
        """Wrapper for :py:func:`shnitsel.analyze.hops.assign_hop_time`."""
        return assign_hop_time(self._obj, hop_types=hop_types, which=which)


class DatasetAccessor(DSManualAccessor):
    _methods = [
        'pca_and_hops',
        'validate',
        'get_spectra',
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
        'calculate_energy_filtranda',
        'filter_by_energy',
        'sanity_check',
        'calculate_bond_length_filtranda',
        'filter_by_length',
        'omit',
        'truncate',
        'transect',
        'write_ase_db',
        'pls_ds',
        'hops_mask_from_active_state',
        'filter_data_at_hops',
        'focus_hops',
        'assign_hop_time',
        'FrameSelector',
        'TrajSelector',
    ]

    @needs(coords_or_vars={'astate', 'atXYZ'})
    def pca_and_hops(self, structure_selection: Union=None, center_mean: bool=False, n_components: int=2) -> Union:
        """Wrapper for :py:func:`shnitsel.analyze.pca.pca_and_hops`."""
        return pca_and_hops(self._obj, structure_selection=structure_selection, center_mean=center_mean, n_components=n_components)

    def validate(self) -> ndarray:
        """Wrapper for :py:func:`shnitsel.data.helpers.validate`."""
        return validate(self._obj)

    @needs(coords={'statecomb', 'time'}, data_vars={'energy', 'fosc'})
    def get_spectra(self, state_selection: shnitsel.filtering.state_selection.StateSelection | None=None, times: Union=None, rel_cutoff: float=0.01) -> Union:
        """Wrapper for :py:func:`shnitsel.analyze.spectra.get_spectra`."""
        return get_spectra(self._obj, state_selection=state_selection, times=times, rel_cutoff=rel_cutoff)

    @needs(dims={'state'})
    def get_per_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.stats.get_per_state`."""
        return get_per_state(self._obj)

    @needs(dims={'state'}, coords={'state'})
    def get_inter_state(self) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.stats.get_inter_state`."""
        return get_inter_state(self._obj)

    @needs(dims={'frame', 'state'}, coords={'time'}, data_vars={'astate'})
    def calc_classical_populations(self) -> "shnitsel.analyze.populations.PopulationStatistics | shnitsel.data.tree.node.TreeNode[PopulationStatistics]":
        """Wrapper for :py:func:`shnitsel.analyze.populations.calc_classical_populations`."""
        return calc_classical_populations(self._obj)

    def default_mol(self, to2D: bool=True, charge: int | float | None=None, molAtomMapNumber: Union=None, atomNote: Union=None, atomLabel: Union=None, silent_mode: bool=False) -> Mol:
        """Wrapper for :py:func:`shnitsel.bridges.default_mol`."""
        return default_mol(self._obj, to2D=to2D, charge=charge, molAtomMapNumber=molAtomMapNumber, atomNote=atomNote, atomLabel=atomLabel, silent_mode=silent_mode)

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
    def sel_trajs(self, trajids_or_mask: Sequence[int] | Sequence[bool], invert: bool=False) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.sel_trajs`."""
        return sel_trajs(self._obj, trajids_or_mask, invert=invert)

    def unstack_trajs(self, fill_value=shnitsel.data.multi_indices.dtype_NA) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.unstack_trajs`."""
        return unstack_trajs(self._obj, fill_value=fill_value)

    def stack_trajs(self) -> DatasetOrArray:
        """Wrapper for :py:func:`shnitsel.data.multi_indices.stack_trajs`."""
        return stack_trajs(self._obj)

    def write_shnitsel_file(self, savepath: str | os.PathLike, complevel: int=9):
        """Wrapper for :py:func:`shnitsel.io.shnitsel.write.write_shnitsel_file`."""
        return write_shnitsel_file(self._obj, savepath, complevel=complevel)

    def calculate_energy_filtranda(self, energy_thresholds: dict[str, float] | shnitsel.clean.filter_energy.EnergyFiltrationThresholds | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.calculate_energy_filtranda`."""
        return calculate_energy_filtranda(self._obj, energy_thresholds=energy_thresholds)

    def filter_by_energy(self, filter_method: Union='truncate', energy_thresholds: dict[str, float] | shnitsel.clean.filter_energy.EnergyFiltrationThresholds | None=None, plot_thresholds: Union=False, plot_populations: Literal=False) -> Optional:
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.filter_by_energy`."""
        return filter_by_energy(self._obj, filter_method=filter_method, energy_thresholds=energy_thresholds, plot_thresholds=plot_thresholds, plot_populations=plot_populations)

    def sanity_check(self, filter_method: Union='truncate', energy_thresholds: dict[str, float] | shnitsel.clean.filter_energy.EnergyFiltrationThresholds | None=None, geometry_thresholds: dict[str, float] | shnitsel.clean.filter_geo.GeometryFiltrationThresholds | None=None, plot_thresholds: Union=False, plot_populations: Literal=False, mol: rdkit.Chem.rdchem.Mol | None=None, drop_empty_trajectories: bool=False) -> Union:
        """Wrapper for :py:func:`shnitsel.clean.sanity_check`."""
        return sanity_check(self._obj, filter_method=filter_method, energy_thresholds=energy_thresholds, geometry_thresholds=geometry_thresholds, plot_thresholds=plot_thresholds, plot_populations=plot_populations, mol=mol, drop_empty_trajectories=drop_empty_trajectories)

    def calculate_bond_length_filtranda(self, geometry_thresholds: dict[str, float] | shnitsel.clean.filter_geo.GeometryFiltrationThresholds | None=None, mol: rdkit.Chem.rdchem.Mol | None=None) -> DataArray:
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.calculate_bond_length_filtranda`."""
        return calculate_bond_length_filtranda(self._obj, geometry_thresholds=geometry_thresholds, mol=mol)

    def filter_by_length(self, filter_method: Union='truncate', geometry_thresholds: dict[str, float] | shnitsel.clean.filter_geo.GeometryFiltrationThresholds | None=None, mol: rdkit.Chem.rdchem.Mol | None=None, plot_thresholds: Union=False, plot_populations: Literal=False) -> Optional:
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.filter_by_length`."""
        return filter_by_length(self._obj, filter_method=filter_method, geometry_thresholds=geometry_thresholds, mol=mol, plot_thresholds=plot_thresholds, plot_populations=plot_populations)

    def omit(self) -> Optional:
        """Wrapper for :py:func:`shnitsel.clean.common.omit`."""
        return omit(self._obj)

    def truncate(self) -> Union:
        """Wrapper for :py:func:`shnitsel.clean.common.truncate`."""
        return truncate(self._obj)

    def transect(self, cutoff_time: float) -> "shnitsel.data.dataset_containers.trajectory.Trajectory | None":
        """Wrapper for :py:func:`shnitsel.clean.common.transect`."""
        return transect(self._obj, cutoff_time)

    @needs(data_vars={'atNames', 'atNums', 'atXYZ', 'energy'})
    def write_ase_db(self, db_path: str, db_format: Optional=None, keys_to_write: Optional=None, preprocess: bool=True, force: bool=False):
        """Wrapper for :py:func:`shnitsel.io.ase.write.write_ase_db`."""
        return write_ase_db(self._obj, db_path, db_format=db_format, keys_to_write=keys_to_write, preprocess=preprocess, force=force)

    def pls_ds(self, xname: str, yname: str, n_components: int=2, common_dim: str | None=None) -> Dataset:
        """Wrapper for :py:func:`shnitsel.analyze.pls.pls_ds`."""
        return pls_ds(self._obj, xname, yname, n_components=n_components, common_dim=common_dim)

    def hops_mask_from_active_state(self, hop_type_selection: Union=None, dim: str | None=None) -> "xarray.core.dataarray.DataArray | shnitsel.data.tree.node.TreeNode[DataArray]":
        """Wrapper for :py:func:`shnitsel.analyze.hops.hops_mask_from_active_state`."""
        return hops_mask_from_active_state(self._obj, hop_type_selection=hop_type_selection, dim=dim)

    def filter_data_at_hops(self, hop_type_selection: Union=None) -> "shnitsel.data.dataset_containers.data_series.DataSeries | xarray.core.dataarray.DataArray | shnitsel.data.tree.node.TreeNode[DataSeries] | shnitsel.data.tree.node.TreeNode[DataArray]":
        """Wrapper for :py:func:`shnitsel.analyze.hops.filter_data_at_hops`."""
        return filter_data_at_hops(self._obj, hop_type_selection=hop_type_selection)

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

