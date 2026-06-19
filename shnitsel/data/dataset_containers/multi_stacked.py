from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Self, Sequence, Literal, Any

import shnitsel
from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
from shnitsel.data.dataset_containers.trajectory import Trajectory
from shnitsel.data.dataset_containers.frames import Frames
import xarray as xr


if TYPE_CHECKING:
    from rdkit.Chem import Mol
    import numpy.typing as npt

    from .multi_layered import MultiSeriesLayered
    from shnitsel.analyze.populations import PopulationStatistics
    from shnitsel.data.tree.node import TreeNode
    from shnitsel.core.typedefs import DataArrayOrVar, DatasetOrArray

    from collections.abc import Collection, Hashable, Iterable, Mapping, Sequence
    from matplotlib.axes._axes import Axes
    from numpy import nan, ndarray
    from os import PathLike

    from shnitsel.analyze.dimred.lda import LDAResult, lda
    from shnitsel.analyze.dimred.pca import PCAResult, pca, pca_and_hops, pca_direct
    from shnitsel.analyze.dimred.pls import pls, pls_ds

    from shnitsel.analyze.populations import PopulationStatistics

    from shnitsel.clean import sanity_check
    from shnitsel.clean.common import TrajectoryOrFrames
    from shnitsel.clean.filter_energy import EnergyFiltrationThresholds
    from shnitsel.clean.filter_geo import GeometryFiltrationThresholds
    from shnitsel.core.typedefs import DataArrayOrVar, DatasetOrArray
    from shnitsel.data.dataset_containers.data_series import DataSeries
    from shnitsel.data.dataset_containers.frames import Frames
    from shnitsel.data.dataset_containers.inter_state import InterState
    from shnitsel.data.dataset_containers.multi_series import MultiSeriesDataset
    from shnitsel.data.dataset_containers.shared import ShnitselDataset
    from shnitsel.data.dataset_containers.trajectory import Trajectory
    from shnitsel.data.multi_indices import dtype_NA
    from shnitsel.data.tree.node import TreeNode
    from shnitsel.data.xr_io_compatibility import SupportsToXrConversion
    from shnitsel.filtering.state_selection import StateSelection
    from shnitsel.filtering.structure_selection import StructureSelection

    from shnitsel.vis.plot.p3mhelpers import frame3D, frames3Dgrid, traj3D, trajs3Dgrid
    from typing import Any, Callable, Dict, List, Literal, Optional, Union
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.groupby import DataArrayGroupBy, DatasetGroupBy

    from typing import TypeAlias

    StructureSelectionEquiv: TypeAlias = (
        StructureSelection
        | str
        | Literal["atoms", "bonds", "angles", "dihedrals", "pyramids", "pwdist", "BLA"]
        | int
        | tuple
        | Collection
        | None
    )


@dataclass
class MultiSeriesStacked(Frames, MultiSeriesDataset):
    """A version of the multi-series dataset where the data is indexed along a sahred `frame` (Multi-index) dimension.
    There is no padding necessary to make the trajectories the same length.
    """

    _layered_repr_cached: "MultiSeriesLayered | None" = None

    def __init__(
        self, framesets: Sequence[Frames | Trajectory | xr.Dataset] | xr.Dataset
    ):
        if isinstance(framesets, xr.Dataset):
            assert 'frame' in framesets.dims, (
                "Stacked dataset must have `frame` dimension"
            )
            assert 'atrajectory' in framesets.coords, (
                "Stacked dataset must have `atrajectory` coordinate"
            )
            assert 'trajectory' in framesets.dims, (
                "Stacked dataset must have `trajectory` dimension"
            )

            MultiSeriesDataset.__init__(self, framesets)
            self._is_multi_trajectory = framesets.sizes['trajectory'] > 1
        else:
            from shnitsel.data.traj_combiner_methods import concat_trajs

            # TODO: FIXME: Stack frames into one single big frameset.
            is_multi_trajectory = False
            if len(framesets) > 1:
                is_multi_trajectory = True
            elif len(framesets) == 1 and framesets[0].is_multi_trajectory:
                is_multi_trajectory = True

            # TODO: FIXME: Make sure that concatenation would work. Convert variables to same unit, etc.

            # Build the concatenated trajectory. May trigger exceptions
            combined_dataset = concat_trajs(framesets)

            Frames.__init__(self, combined_dataset)
            MultiSeriesDataset.__init__(self, framesets, combined_dataset)
            self._is_multi_trajectory = is_multi_trajectory

    @property
    def grouping_dimension(self) -> str:
        return 'atrajectory'

    @cached_property
    def as_layered(self) -> "MultiSeriesLayered":
        """Get a layered representation of the stacked datasets in this object

        Returns
        -------
        MultiSeriesLayered
            The converted (or extracted from cache) layered version of this multi-data dataset.
        """
        from .multi_layered import MultiSeriesLayered
        from shnitsel.data.dataset_containers import wrap_dataset

        if self._layered_repr_cached is not None and isinstance(
            self._layered_repr_cached, MultiSeriesLayered
        ):
            if self._layered_repr_cached._stacked_repr_cached is not self:
                self._layered_repr_cached._stacked_repr_cached = self

            return self._layered_repr_cached


        # Use multi-index stack/unstack logic
        from shnitsel.data.multi_indices import unstack_trajs
        ds: xr.Dataset = self.dataset
        ds_unstacked = unstack_trajs(ds)
        tmp_res = MultiSeriesLayered(ds_unstacked)
        tmp_res._basis_data = self._basis_data

        # if self._basis_data is not None:
        #     tmp_res = MultiSeriesLayered(self._basis_data)
        # else:
        #     ds: xr.Dataset = self.dataset
        #     datasets: Sequence[Frames | Trajectory] = [
        #         wrap_dataset(
        #             ds.sel(trajectory=id, atrajectory=id)
        #             .drop_dims(['trajectory', 'atrajectory'], errors="ignore")
        #             .drop_vars('atrajectory'),
        #             expected_types=Trajectory | Frames,
        #         )
        #         for id in ds.coords['trajectory'].values
        #     ]

        #     tmp_res = MultiSeriesLayered(datasets)
        # Set self as cached result of the inverse conversion
        tmp_res._stacked_repr_cached = self
        self._layered_repr_cached = tmp_res
        return tmp_res

    @cached_property
    def as_stacked(self) -> Self:
        return self

    @classmethod
    def get_type_marker(cls) -> str:
        return "shnitsel::MultiSeriesStacked"

    # Wrappers for functions that take Dataset-like input

    def get_spectra(
        self,
        state_selection: "(StateSelection | None)" = None,
        times: '(Iterable | Literal["all"] | None)' = None,
        rel_cutoff: float = 0.01,
    ) -> "(DataArray | TreeNode[Any, DataArray])":
        """Wrapper for :py:func:`shnitsel.analyze.spectra.get_spectra`."""
        return shnitsel.analyze.spectra.get_spectra(
            self, state_selection=state_selection, times=times, rel_cutoff=rel_cutoff
        )

    def calc_classical_populations(
        self,
    ) -> "(PopulationStatistics | TreeNode[Any, PopulationStatistics])":
        """Wrapper for :py:func:`shnitsel.analyze.populations.calc_classical_populations`."""

        return shnitsel.analyze.populations.calc_classical_populations(self)

    def default_mol(
        self,
        to2D: bool = True,
        charge: (int | float | None) = None,
        molAtomMapNumber: (list | Literal[True] | None) = None,
        atomNote: (list | Literal[True] | None) = None,
        atomLabel: (list | Literal[True] | None) = None,
        silent_mode: bool = False,
    ) -> "Mol":
        """Wrapper for :py:func:`shnitsel.bridges.default_mol`."""
        return shnitsel.bridges.default_mol(
            self,
            to2D=to2D,
            charge=charge,
            molAtomMapNumber=molAtomMapNumber,
            atomNote=atomNote,
            atomLabel=atomLabel,
            silent_mode=silent_mode,
        )

    def assign_levels(
        self,
        levels: dict[str, "npt.ArrayLike"] | None = None,
        **levels_kwargs: "npt.ArrayLike",
    ) -> "DatasetOrArray":
        """Wrapper for :py:func:`shnitsel.data.multi_indices.assign_levels`."""
        return shnitsel.data.multi_indices.assign_levels(
            self, levels=levels, **levels_kwargs
        )

    def mgroupby(self, levels: Sequence[str]) -> "DataArrayGroupBy | DatasetGroupBy":
        """Wrapper for :py:func:`shnitsel.data.multi_indices.mgroupby`."""
        return mgroupby(self, levels)

    def msel(self, **kwargs) -> "DatasetOrArray":
        """Wrapper for :py:func:`shnitsel.data.multi_indices.msel`."""
        return shnitsel.data.multi_indices.mselmsel(self, **kwargs)

    def sel_trajs(
        self, trajids_or_mask: Sequence[int] | Sequence[bool], invert: bool = False
    ) -> "DatasetOrArray":
        """Wrapper for :py:func:`shnitsel.data.multi_indices.sel_trajs`."""
        return shnitsel.data.multi_indices.sel_trajs(
            self, trajids_or_mask, invert=invert
        )

    def unstack_trajs(
        self, fill_value="dtype_NA"
    ) -> "DatasetOrArray | ShnitselDataset":
        """Wrapper for :py:func:`shnitsel.data.multi_indices.unstack_trajs`."""
        return shnitsel.data.multi_indices.unstack_trajs(self, fill_value=fill_value)

    def stack_trajs(self) -> "DatasetOrArray":
        """Wrapper for :py:func:`shnitsel.data.multi_indices.stack_trajs`."""
        return shnitsel.data.multi_indices.stack_trajs(self)

    def write_shnitsel_file(
        self,
        savepath: "(str | PathLike)",
        complevel: int = 9,
        output_engine: (str | Literal["netcdf4", "h5netcdf"]) = 'h5netcdf',
        output_format: (Literal["NETCDF4", "NETCDF4_CLASSIC"] | None) = None,
    ):
        """Wrapper for :py:func:`shnitsel.io.shnitsel.write.write_shnitsel_file`."""
        return shnitsel.io.shnitsel.write.write_shnitsel_file(
            self,
            savepath,
            complevel=complevel,
            output_engine=output_engine,
            output_format=output_format,
        )

    def calculate_energy_filtranda(
        self, energy_thresholds: "(dict | EnergyFiltrationThresholds | None)" = None
    ) -> "DataArray":
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.calculate_energy_filtranda`."""
        return shnitsel.clean.filter_energy.calculate_energy_filtranda(
            self, energy_thresholds=energy_thresholds
        )

    def filter_by_energy(
        self,
        filter_method: (Literal["truncate", "omit", "annotate"] | float) = 'truncate',
        energy_thresholds: "(dict | EnergyFiltrationThresholds | None)" = None,
        plot_thresholds: (bool | Sequence) = False,
        plot_populations: Literal["independent", "intersections", False] = False,
    ) -> "TreeNode | TrajectoryOrFrames | None":
        """Wrapper for :py:func:`shnitsel.clean.filter_energy.filter_by_energy`."""
        return shnitsel.clean.filter_energy.filter_by_energy(
            self._obj,
            filter_method=filter_method,
            energy_thresholds=energy_thresholds,
            plot_thresholds=plot_thresholds,
            plot_populations=plot_populations,
        )

    def sanity_check(
        self,
        filter_method: (Literal["truncate", "omit", "annotate"] | float) = 'truncate',
        energy_thresholds: "(dict | EnergyFiltrationThresholds | None)" = None,
        geometry_thresholds: "(dict | GeometryFiltrationThresholds | None)" = None,
        plot_thresholds: (bool | Sequence) = False,
        plot_populations: Literal["independent", "intersections", False] = False,
        mol: "(Mol | None)" = None,
        drop_empty_trajectories: bool = False,
    ) -> "TreeNode | TrajectoryOrFrames | None":
        """Wrapper for :py:func:`shnitsel.clean.sanity_check`."""
        return shnitsel.clean.sanity_check(
            self,
            filter_method=filter_method,
            energy_thresholds=energy_thresholds,
            geometry_thresholds=geometry_thresholds,
            plot_thresholds=plot_thresholds,
            plot_populations=plot_populations,
            mol=mol,
            drop_empty_trajectories=drop_empty_trajectories,
        )

    def calculate_bond_length_filtranda(
        self,
        geometry_thresholds: "(dict | GeometryFiltrationThresholds | None)" = None,
        mol: "(Mol | None)" = None,
    ) -> "DataArray":
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.calculate_bond_length_filtranda`."""
        return shnitsel.clean.filter_geo.calculate_bond_length_filtranda(
            self, geometry_thresholds=geometry_thresholds, mol=mol
        )

    def filter_by_length(
        self,
        filter_method: (Literal["truncate", "omit", "annotate"] | float) = 'truncate',
        geometry_thresholds: "(dict | GeometryFiltrationThresholds | None)" = None,
        mol: "(Mol | None)" = None,
        plot_thresholds: (bool | Sequence) = False,
        plot_populations: Literal["independent", "intersections", False] = False,
    ) -> "TreeNode | TrajectoryOrFrames | None":
        """Wrapper for :py:func:`shnitsel.clean.filter_geo.filter_by_length`."""
        return shnitsel.clean.filter_geo.filter_by_length(
            self,
            filter_method=filter_method,
            geometry_thresholds=geometry_thresholds,
            mol=mol,
            plot_thresholds=plot_thresholds,
            plot_populations=plot_populations,
        )

    def omit(self) -> "TrajectoryOrFrames | None":
        """Wrapper for :py:func:`shnitsel.clean.common.omit`."""
        return shnitsel.clean.omit(self)

    def truncate(self) -> "TrajectoryOrFrames | Trajectory | Frames | None":
        """Wrapper for :py:func:`shnitsel.clean.common.truncate`."""
        return shnitsel.clean.common.truncate(self)

    def transect(self, cutoff_time: float) -> "Trajectory | None":
        """Wrapper for :py:func:`shnitsel.clean.common.transect`."""
        return shnitsel.clean.common.transect(self, cutoff_time)

    def write_ase_db(
        self,
        db_path: str,
        db_format: (Literal["spainn", "schnet"] | None) = None,
        keys_to_write: "(Collection | None)" = None,
        preprocess: bool = True,
        force: bool = False,
    ):
        """Wrapper for :py:func:`shnitsel.io.ase.write.write_ase_db`."""
        return shnitsel.io.ase.write.write_ase_db(
            self,
            db_path,
            db_format=db_format,
            keys_to_write=keys_to_write,
            preprocess=preprocess,
            force=force,
        )

    def pls_ds(
        self,
        xname: str,
        yname: str,
        n_components: int = 2,
        common_dim: (str | None) = None,
    ) -> "Dataset":
        """Wrapper for :py:func:`shnitsel.analyze.dimred.pls.pls_ds`."""
        return shnitsel.analyze.dimred.pls.pls_ds(
            self, xname, yname, n_components=n_components, common_dim=common_dim
        )

    def hops_mask_from_active_state(
        self,
        hop_type_selection: "(StateSelection | Sequence | Sequence | str | None)" = None,
        dim: "(str | Hashable | None)" = None,
    ) -> "DataArray | TreeNode[Any, DataArray]":
        """Wrapper for :py:func:`shnitsel.analyze.hops.hops_mask_from_active_state`."""
        return shnitsel.analyze.hops.hops_mask_from_active_state(
            self, hop_type_selection=hop_type_selection, dim=dim
        )

    def filter_data_at_hops(
        self,
        hop_type_selection: "(StateSelection | Sequence | Sequence | str | None)" = None,
    ) -> (
        "DataSeries | DataArray | TreeNode[Any, DataSeries] | TreeNode[Any, DataArray]"
    ):
        """Wrapper for :py:func:`shnitsel.analyze.hops.filter_data_at_hops`."""
        return shnitsel.analyze.hops.filter_data_at_hops(
            self, hop_type_selection=hop_type_selection
        )

    def focus_hops(
        self,
        hop_type_selection: "(StateSelection | Sequence | Sequence | str | None)" = None,
        window: (slice | None) = None,
    ) -> "Dataset | DataArray | DataSeries | TreeNode":
        """Wrapper for :py:func:`shnitsel.analyze.hops.focus_hops`."""
        return shnitsel.analyze.hops.focus_hops(
            self, hop_type_selection=hop_type_selection, window=window
        )

    def assign_hop_time(
        self,
        hop_type_selection: "(StateSelection | Sequence | Sequence | str | None)" = None,
        which: Literal["first", "last"] = 'last',
    ) -> "Dataset | DataArray | DataSeries | TreeNode":
        """Wrapper for :py:func:`shnitsel.analyze.hops.assign_hop_time`."""
        return shnitsel.analyze.hops.assign_hop_time(
            self, hop_type_selection=hop_type_selection, which=which
        )

    def FrameSelector(
        self,
        data_var=None,
        dim=None,
        xname=None,
        yname=None,
        title='',
        allowed_ws_origin=None,
        webgl=True,
    ):
        """Wrapper for :py:func:`shnitsel.vis.plot.select.FrameSelector`."""
        return shnitsel.vis.plot.select.FrameSelector(
            self,
            data_var=data_var,
            dim=dim,
            xname=xname,
            yname=yname,
            title=title,
            allowed_ws_origin=allowed_ws_origin,
            webgl=webgl,
        )

    def TrajSelector(
        self,
        data_var=None,
        dim=None,
        xname=None,
        yname=None,
        title='',
        allowed_ws_origin=None,
        webgl=True,
    ):
        """Wrapper for :py:func:`shnitsel.vis.plot.select.TrajSelector`."""
        return shnitsel.vis.plot.select.TrajSelector(
            self,
            data_var=data_var,
            dim=dim,
            xname=xname,
            yname=yname,
            title=title,
            allowed_ws_origin=allowed_ws_origin,
            webgl=webgl,
        )

    # Wrappers for functions that act on atXYZ

    def distance(self, i_index: (int | list), j_index: (int | list)) -> "DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.distances.distance`."""
        return shnitsel.geo.geocalc_.distances.distance(self['atXYZ'], i_index, j_index)

    def angle(
        self,
        a_index: (int | list),
        b_index: (int | list),
        c_index: (int | list),
        angles: Literal["deg", "rad", "trig"] = 'deg',
    ) -> "DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.angles.angle`."""
        return shnitsel.geo.geocalc_.angles.angle(
            self['atXYZ'], a_index, b_index, c_index, angles=angles
        )

    def dihedral(
        self,
        a_index: (int | list),
        b_index: (int | list),
        c_index: (int | list),
        d_index: (int | list),
        angles: Literal["deg", "rad", "trig"] = 'deg',
        full: bool = False,
    ) -> " DataArray | tuple":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.dihedrals.dihedral`."""
        return shnitsel.geo.geocalc_.dihedrals.dihedral(
            self['atXYZ'], a_index, b_index, c_index, d_index, angles=angles, full=full
        )

    def pyramidalization_angle(
        self,
        x_index: (int | list),
        a_index: (int | list),
        b_index: (int | list),
        c_index: (int | list),
        angles: Literal["deg", "rad", "trig"] = 'deg',
    ) -> "DataArray | tuple":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.pyramids.pyramidalization_angle`."""
        return shnitsel.geo.geocalc_.pyramids.pyramidalization_angle(
            self['atXYZ'], x_index, a_index, b_index, c_index, angles=angles
        )

    def get_bats(
        self,
        structure_selection: "StructureSelectionEquiv" = None,
        default_features: Sequence = ['bonds', 'angles', 'dihedrals'],
        signed: bool = False,
        angles: Literal["deg", "rad", "trig"] = 'deg',
    ) -> "DataArray | TreeNode[Any, DataArray]":
        """Wrapper for :py:func:`shnitsel.geo.geocalc.get_bats`."""
        return shnitsel.geo.geocalc.get_bats(
            self['atXYZ'],
            structure_selection=structure_selection,
            default_features=default_features,
            signed=signed,
            angles=angles,
        )

    def get_distances(
        self,
        structure_selection: "StructureSelectionEquiv" = None,
    ) -> "TreeNode[Any, DataArray] | DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.distances.get_distances`."""
        return shnitsel.geo.geocalc_.distances.get_distances(
            self['atXYZ'], structure_selection=structure_selection
        )

    def get_angles(
        self,
        structure_selection: "StructureSelectionEquiv" = None,
        angles: Literal["deg", "rad", "trig"] = 'deg',
        signed: bool = True,
    ) -> "TreeNode[Any, DataArray] | DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.angles.get_angles`."""
        return shnitsel.geo.geocalc_.angles.get_angles(
            self['atXYZ'],
            structure_selection=structure_selection,
            angles=angles,
            signed=signed,
        )

    def get_dihedrals(
        self,
        structure_selection: "StructureSelectionEquiv" = None,
        angles: Literal["deg", "rad", "trig"] = 'deg',
        signed: bool = True,
    ) -> "TreeNode[Any, DataArray] | DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.dihedrals.get_dihedrals`."""
        return shnitsel.geo.geocalc_.dihedrals.get_dihedrals(
            self['atXYZ'],
            structure_selection=structure_selection,
            angles=angles,
            signed=signed,
        )

    def get_pyramidalization(
        self,
        structure_selection: "StructureSelectionEquiv" = None,
        angles: Literal["deg", "rad", "trig"] = 'deg',
        signed: bool = True,
    ) -> "TreeNode[Any, DataArray] | DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.pyramids.get_pyramidalization`."""
        return shnitsel.geo.geocalc_.pyramids.get_pyramidalization(
            self['atXYZ'],
            structure_selection=structure_selection,
            angles=angles,
            signed=signed,
        )

    def get_max_chromophor_BLA(
        self,
        structure_selection: "StructureSelectionEquiv" = None,
        SMARTS: (str | None) = None,
        num_double_bonds: (int | None) = None,
        allowed_chain_elements: str = '#6,#7,#8,#15,#16',
        max_considered_BLA_double_bonds: int = 50,
    ) -> "TreeNode[Any, DataArray] | DataArray":
        """Wrapper for :py:func:`shnitsel.geo.geocalc_.bla_chromophor.get_max_chromophor_BLA`."""
        return shnitsel.geo.geocalc_.bla_chromophor.get_max_chromophor_BLA(
            self['atXYZ'],
            structure_selection=structure_selection,
            SMARTS=SMARTS,
            num_double_bonds=num_double_bonds,
            allowed_chain_elements=allowed_chain_elements,
            max_considered_BLA_double_bonds=max_considered_BLA_double_bonds,
        )

    def trajs3Dgrid(self, trajids: (list | None) = None, loop: str = 'forward'):
        """Wrapper for :py:func:`shnitsel.vis.plot.p3mhelpers.trajs3Dgrid`."""
        return shnitsel.vis.plot.p3mhelpers.trajs3Dgrid(
            self['atXYZ'], trajids=trajids, loop=loop
        )

    def traj_vmd(self, groupby: str = 'atrajectory'):
        """Wrapper for :py:func:`shnitsel.vis.vmd.traj_vmd`."""
        return shnitsel.vis.vmd.traj_vmd(self['atXYZ'], groupby=groupby)
