from dataclasses import dataclass
from functools import cached_property
import logging
from matplotlib.axes import Axes
from sklearn.decomposition import PCA as sk_PCA
from tqdm import tqdm
import xarray as xr
import numpy as np
import rdkit.Chem as rdchem
import matplotlib.pyplot as plt
from logging import info, warning
from timeit import default_timer as timer

from matplotlib.figure import Figure, SubFigure

from shnitsel.analyze.pca import PCAResult
from shnitsel.analyze.populations import (
    PopulationStatistics,
    calc_classical_populations,
)
from shnitsel.bridges import construct_default_mol
from shnitsel.analyze import stats
from shnitsel.core.typedefs import (
    AtXYZ,
    SpectraDictType,
    StateCombination,
)
from shnitsel.data.dataset_containers import (
    Trajectory,
    Frames,
    InterState,
    PerState,
    wrap_dataset,
)
from shnitsel.filtering.helpers import (
    _get_default_state_selection,
    _get_default_structure_selection,
)
from shnitsel.filtering.state_selection import StateSelection, StateSelectionDescriptor
from shnitsel.filtering.structure_selection import (
    StructureSelection,
    StructureSelectionDescriptor,
)
from shnitsel.vis.datasheet.figures.energy_bands import plot_energy_bands
from shnitsel.vis.datasheet.figures.soc_trans_hist import (
    plot_separated_spectra_and_soc_dip_hists,
    plot_separated_spectra_and_soc_dip_hists_groundstate,
    single_dip_trans_hist,
    single_soc_trans_hist,
)
from shnitsel.vis.colormaps import st_grey
from shnitsel.vis.plot.common import inlabel, outlabel

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from shnitsel.analyze.spectra import (
    get_fosc_gauss_broadened,
    get_spectra,
    get_spectra_groups,
)

from .figures.common import centertext, label_plot_grid
from .figures.per_state_hist import plot_per_state_histograms
from .figures.time import plot_timeplots

# from .figures.dip_trans_hist import (
#     plot_separated_spectra_and_hists,
#     plot_separated_spectra_and_hists_groundstate,
# )
from .figures.nacs_hist import plot_nacs_histograms
from ..plot.pca_biplot import plot_noodleplot
from .figures.structure import format_inchi, plot_pca_structure, plot_structure


@dataclass
class DatasheetPage:
    """Class encapsulating the data generation and plotting of a single page of a datasheet.

    Involves plotting of inter-state graphs like spectra, statistics of energy deltas and transitional dipoles,
    histograms of per-state and inter-state properties, etc.
    """

    state_selection: StateSelection
    state_selection_provided: bool = False
    structure_selection: StructureSelection | None = None
    structure_selection_provided: bool = False

    spectra_times: list[int | float] | np.ndarray | None = None
    charge: int = 0
    structure_skeletal: bool = False
    name: str = ''

    def __init__(
        self,
        data: xr.Dataset | Trajectory | Frames | Self,
        state_selection: StateSelection | StateSelectionDescriptor | None = None,
        structure_selection: StructureSelection
        | StructureSelectionDescriptor
        | None = None,
        *,
        spectra_times: list[int | float] | np.ndarray | None = None,
        col_state: list | None = None,
        col_inter: list | None = None,
    ):
        """Create a new DataSheet page, initializating the data of the class and preparing the plotting later on

        Parameters
        ----------
        data : ShnitselDataset | Self
            Either a Shnitsel dataset data object from which to generate the statistical information visualized in this DataSheet page
            or another DataSheetPage that should be copied.
        state_selection : StateSelection | StateSelectionDescriptor, optional
            Optional parameter to specify a subset of states and state combinations that may be considered for the dataset.
            Will be generated if not provided.
        structure_selection: StructureSelection | StructureSelectionDescriptor, optional
            Optional parameter to limit the PCA plot and analysis to a specific subset of the structure.
            Will be generated if not provided.
        spectra_times : list[int  |  float] | np.ndarray, optional
            Sequence of times to calculate spectra at. Defaults to None.
        col_state : list, optional
            A list of colors to use for the states. Defaults to default shnitsel colors.
        col_inter : list, optional
            A list of colors to use for state combinations. Defaults to default shnitsel colors.

        Raises
        ------
        TypeError
            If wrong type of `data` parameter is provided.
        ValueError
            If the wrong number of colors for the states is provided.
        ValueError
            If the wrong number of colors for state transitions is provided.

        Returns
        -------
        DatasheetPage
            The constructed (or copied) DataSheetPage
        """
        if isinstance(data, DatasheetPage):
            self._copy_data(old=data)
            return
        else:
            if isinstance(data, xr.Dataset):
                try:
                    data = Trajectory(data)
                except:
                    try:
                        data = Frames(data)
                    except:
                        raise TypeError(
                            "Data provided to datasheet page is neither of the right shape to be a trajectory nor of the type to be considered frames data."
                        )

            if isinstance(data, Frames):
                self.frames = data
            elif isinstance(data, Trajectory):
                self.frames = data.as_frames
            else:
                raise TypeError("Neither DatasheetPage nor frames/Trajectory given.")

        assert isinstance(self.frames, Frames)

        if spectra_times is not None:
            self.spectra_times = spectra_times
        elif 'time' not in self.frames.coords:
            logging.warning(
                "No 'time' variable found. Have ICONDs been passed as frames?"
            )
        elif self.frames is not None:
            max_time = self.frames.coords['time'].max().item()
            # self.spectra_times = [max_time * i / 40 for i in range(5)]
            # self.spectra_times += [max_time * i / 20 for i in range(5)]
            # self.spectra_times += [max_time * i / 3 for i in range(4)]
            self.spectra_times = [max_time * 1 / 20]
            self.spectra_times += [max_time * 1 / 10]
            self.spectra_times += [max_time * 1 / 5]
            self.spectra_times += [max_time * 1 / 3]

        # Initialize state selection or use provided selection
        if state_selection is not None:
            state_selection = _get_default_state_selection(
                state_selection, state_source=self.frames
            )
            self.state_selection_provided = True
            self.state_selection = state_selection
        else:
            self.state_selection_provided = False
            self.state_selection = StateSelection.init_from_dataset(self.frames)

        # Initialize feature selection or use provided selection
        if structure_selection is not None:
            position_data: xr.DataArray
            charge_info: int | None
            if isinstance(self.frames, xr.DataArray):
                position_data = self.frames
                charge_info = None
            else:
                wrapped_ds = wrap_dataset(self.frames, (Trajectory | Frames))
                position_data = wrapped_ds.atXYZ
                charge_info = int(wrapped_ds.charge)

            structure_selection = _get_default_structure_selection(
                structure_selection,
                atXYZ_source=position_data,
                default_levels=['atoms', 'bonds', 'angles', 'dihedrals', 'pyramids'],
                charge_info=charge_info,
            )
            self.structure_selection_provided = True
            self.structure_selection = structure_selection
        else:
            self.structure_selection = None

        # print(self.frames)

        # nstates = self.frames.sizes['state']
        # if col_state is not None:
        #     assert (ncols := len(col_state)) >= nstates, (
        #         f"`col_state` has {ncols} colors, "
        #         f"but should contain one color for each of the {nstates} states"
        #     )
        #     self.col_state = col_state[:nstates]
        # elif nstates <= 3:
        #     # SHNITSEL-colours
        #     self.col_state = ['#4DAD15', '#AD2915', '#7515AD'][:nstates]
        # elif nstates <= 10:
        #     cmap = plt.get_cmap('tab10')
        #     self.col_state = [mpl.colors.rgb2hex(c) for c in cmap.colors][:nstates]  # type: ignore
        # elif nstates <= 20:
        #     cmap = plt.get_cmap('tab20')
        #     self.col_state = [mpl.colors.rgb2hex(c) for c in cmap.colors][:nstates]  # type: ignore
        # else:
        #     raise ValueError(
        #         f"These data have {nstates} states. "
        #         "When passing data with more than 10 states, please "
        #         "also pass an appropriate colormap to `col_state`."
        #     )

        # ncombs = self.frames.sizes['statecomb']
        # if col_inter is not None:
        #     assert (ncols := len(col_inter)) == ncombs, (
        #         f"`col_inter` has {ncols} colors, "
        #         f"but should contain one color for each of the {ncombs} state combinations"
        #     )
        #     self.col_inter = col_inter
        # elif ncombs <= 3:
        #     self.col_inter = col_inter or ['#2c3e50', '#C4A000', '#7E5273'][:ncombs]
        # elif ncombs <= 10:
        #     # TODO: choose colours distinct from per_state colours
        #     cmap = plt.get_cmap('tab10')
        #     self.col_inter = [mpl.colors.rgb2hex(c) for c in cmap.colors][:ncombs]  # type: ignore
        # elif ncombs <= 20:
        #     cmap = plt.get_cmap('tab20')
        #     self.col_inter = [mpl.colors.rgb2hex(c) for c in cmap.colors][:ncombs]  # type: ignore
        # else:
        #     print(self.frames.statecomb.values)
        #     raise ValueError(
        #         f"These data have {ncombs} state combinations. "
        #         "When passing data with more than 10 state combinations, please "
        #         "also pass an appropriate colormap to `col_inter`."
        #     )

        # Automatically find colors for the states and state combinations
        self.state_selection.auto_assign_colors()

        self.can = {}

        # print(self.frames['state_charges'])
        self.charge = int(np.round(self.frames.charge))

        def check(*ks):
            return all(
                self.frames.has_variable(k) or self.frames.has_coordinate(k) for k in ks
            )

        self.can['per_state_histograms'] = check('energy', 'forces', 'dip_trans')
        self.can['separated_spectra_and_hists'] = self.frames.has_variable(
            'energy'
        )  # check('dip_trans', 'time')
        self.can['noodle'] = check('atXYZ', 'state', 'time')
        self.can['structure'] = ('smiles_map' in self.frames.attrs) or check('atXYZ')
        self.can['nacs_histograms'] = check('nacs') and (
            check('energy') or check('forces')
        )
        self.can['timeplots'] = check('time', 'astate')

        try:
            self.name = self.frames.attrs['long_name']
        except KeyError:
            pass

        return None

    def _copy_data(self, old: Self):
        """Copy data from the other DataSheetPage into this entity's fields.

        Parmaeters
        ----------
        old : Self
            The DataSheetPage to create a copy of.
        """
        self.spectra_times = old.spectra_times
        self.structure_selection = old.structure_selection
        self.structure_selection_provided = old.structure_selection_provided
        self.state_selection = old.state_selection
        self.state_selection_provided = old.state_selection_provided
        # self.col_state = old.col_state
        # self.col_inter = old.col_inter
        self.name = old.name
        self.charge = old.charge
        self.structure_skeletal = old.structure_skeletal
        self.per_state = old.per_state
        self.inter_state = old.inter_state
        self.pops = old.pops
        self.delta_E = old.delta_E
        self.fosc_time = old.fosc_time
        self.spectra = old.spectra
        self.spectra_groups = old.spectra_groups
        self.spectra_ground = old.spectra_ground
        self.spectra_excited = old.spectra_excited
        self.noodle = old.noodle
        self.hops = old.hops
        self.structure_atXYZ = old.structure_atXYZ
        self.mol = old.mol
        self.mol_skeletal = old.mol_skeletal
        self.smiles = old.smiles
        self.inchi = old.inchi

    @cached_property
    def per_state(self) -> PerState:
        """Get per-state data for the underlying dataset and cache it for repeated use.

        Returns
        -------
        PerState
            Per-state data of the self.frames object. (Energies, permanent dipoles)
        """
        start = timer()
        per_state = self.frames.per_state
        # per_state['_color'] = 'state', self.col_state
        end = timer()
        info(f"cached per_state in {end - start} s")
        return per_state

    @cached_property
    def inter_state(self) -> InterState:
        """Inter-state (state-transition) data of the underlying self.frames object

        Returns
        -------
        InterState
            Inter-state properties of the underlying data. (delta_energies, transition dipoles, SOCs, NACs, fosc)
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        inter_state = self.frames.inter_state
        # inter_state['_color'] = 'statecomb', self.col_inter

        # Calculate fosc if missing and conditions met
        # if (
        #     "fosc" not in inter_state
        #     and 'dip_trans' in inter_state
        #     and "energy_interstate" in inter_state
        # ):
        #     inter_state = assign_fosc(inter_state)

        for var, tex in [
            ('energy', r"$\Delta E$"),
            ('nacs', r"$\|\mathrm{NAC}_{i,j}\|_2$"),
            ('dip_trans', r"$\|\mathbf{\mu}_{i,j}\|_2$"),
            ('fosc', r"$f_\mathrm{osc}$"),
        ]:
            try:
                inter_state.dataset[var].attrs['tex'] = tex
            except KeyError:
                pass
        end = timer()
        info(f"cached inter_state in {end - start} s")
        return inter_state

    @cached_property
    def pops(self) -> PopulationStatistics:
        """Population data for the underlying self.frames

        Returns
        -------
        PopulationStatistics
            Population data encapsulated in an object holding absolute and relative population data
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        pops = calc_classical_populations(self.frames)
        # pops['_color'] = 'state', self.col_state
        end = timer()
        info(f"cached pops in {end - start} s")
        return pops

    @cached_property
    def delta_E(self) -> xr.Dataset:
        """Energy deltas between different states in a dataset.

        Returns
        -------
        xr.Dataset
            Dataset holding 'energy_interstate' variable.
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        res = stats.time_grouped_confidence_interval(self.inter_state.energy_interstate)
        # res['_color'] = 'statecomb', self.col_inter
        res.attrs['tex'] = r"$\Delta E$"
        end = timer()
        info(f"cached delta_E in {end - start} s")
        return res

    @cached_property
    def fosc_time(self) -> xr.Dataset | None:
        """Strength of oscillator/transition rate between states at different points in time.

        Returns
        -------
        xr.Dataset | None
            Either the f_osc data (with confidence intervals) or None if not sufficient data in self.frames to calculate it.
        """
        start = timer()
        if self.inter_state.has_variable('fosc'):
            res = stats.time_grouped_confidence_interval(self.inter_state.fosc)
            # res['_color'] = 'statecomb', self.col_inter
            res.attrs['tex'] = r"$f_\mathrm{osc}$"
        else:
            res = None
        end = timer()
        info(f"cached fosc_time in {end - start} s")
        return res

    @cached_property
    def spectra(self) -> SpectraDictType:
        """Spectral statistics of the self.frames object.

        Returns
        -------
        SpectraDictType
            The spectral information per state transition
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        if self.frames.has_variable('dip_trans'):
            res = get_spectra(self.inter_state, times=self.spectra_times)
        else:
            res = {}
        end = timer()
        info(f"cached spectra in {end - start} s")
        return res

    @cached_property
    def spectra_groups(
        self,
    ) -> tuple[
        SpectraDictType,
        SpectraDictType,
    ]:
        """Get different spectral groups for ground-state transitions and for excited-state transitions

        Returns
        -------
        tuple[SpectraDictType, SpectraDictType]
            One spectral dict per ground-state or excited-state transitions.
        """
        start = timer()
        # TODO: FIXME: Use state selection for split
        if self.frames.has_variable('dip_trans'):
            res = get_spectra_groups(self.spectra)
        else:
            res = ({}, {})

        end = timer()
        info(f"cached spectra_groups in {end - start} s")
        return res

    @cached_property
    def spectra_ground(self) -> SpectraDictType:
        """Extracted spectral information of only the ground-state transitions

        Returns
        -------
        SpectraDictType
            Extracted spectral information of only the ground-state transitions
        """
        return self.spectra_groups[0]

    @cached_property
    def spectra_excited(self) -> SpectraDictType:
        """Extracted spectral information of only the excited-state transitions

        Returns
        -------
        SpectraDictType
            Extracted spectral information of only the excited-state transitions
        """
        return self.spectra_groups[1]

    @cached_property
    def pca_full_data(self) -> PCAResult:
        """Get full PCA result with PCA detailed info.

        Returns
        -------
        xr.DataArray
            The pairwise distance PCA results
        """
        from shnitsel.analyze.pca import pca

        start = timer()
        res = pca(self.frames, structure_selection=self.structure_selection)
        end = timer()
        info(f"cached pca_data in {end - start} s")
        return res

    @cached_property
    def pca_data(self) -> xr.DataArray:
        """Noodle plot source data derived from principal component analysis (PCA) on the full data in self.frames using only pairwise distances.

        Returns
        -------
        xr.DataArray
            The pairwise distance PCA results
        """
        return self.pca_full_data.projected_inputs

    @cached_property
    def pca_info(self) -> sk_PCA:
        """Detailed PCA decomposition information

        Returns
        -------
        sk_PCA
            The sklearn.decomposition.PCA object containing all PCA information.
        """
        return self.pca_full_data.fitted_pca_object

    @cached_property
    def pca_explanation(self) -> dict[str, str]:
        """Provide an explanation for the components of the PCA components.

        Returns
        -------
        dict[str, str]
            Result of the explanation process.
        """
        res: dict[str, str] = {}

        pca_components = self.pca_full_data.principal_components

        logging.debug(f"{self.pca_full_data.inputs=}")

        for i, component in enumerate(pca_components):
            logging.debug(f"{i=}: {component=}")

        return res

    @cached_property
    def hops(self) -> xr.DataArray:
        """The PCA plots at the hopping points

        Returns:
            xr.DataArray: PCA data at the hopping points
        """
        from shnitsel.analyze.hops import hops_mask_from_active_state

        mask = hops_mask_from_active_state(self.frames)
        return mask  # self.pca_data[mask]

    @cached_property
    def structure_atXYZ(self) -> AtXYZ:
        """Structure/Position data in the first frame/timestep of the trajectory

        Returns
        -------
        AtXYZ
            Positional data.
        """
        leading_dim_name = self.frames.leading_dim
        return self.frames.atXYZ.isel({leading_dim_name: 0})

    @cached_property
    def mol(self) -> rdchem.Mol:
        """Property to get an rdkit Mol object from the structural data

        Returns
        -------
        rdkit.Chem.Mol
            Molecule object representing the structure in the first frame
        """
        # # TODO: FIXME: Shouldn't this be a private attribute prefixed with `__` ?
        # if 'smiles_map' in self.frames['atXYZ'].attrs:
        #     mol = numbered_smiles_to_mol(
        #         self.frames['atXYZ'].attrs['smiles_map']
        #     )
        #     for atom in mol.GetAtoms():
        #         atom.ClearProp("molAtomMapNumber")
        #         atom.SetProp("atomNote", str(atom.GetIdx()))
        #     return mol
        # else:
        return construct_default_mol(self.frames.dataset, charge=self.frames.charge)

    @cached_property
    def mol_skeletal(self) -> rdchem.Mol:
        """Skeletal representation of the the rdkit.Chem.Mol representation of the structure

        Returns
        -------
        rdkit.Chem.Mol
            Molecule object representing the skeletal structure (no H atoms) in the first frame
        """
        mol = rdchem.Mol(self.mol)
        return rdchem.RemoveHs(mol)

    @cached_property
    def smiles(self) -> str:
        """Smiles representation of the skeletal molecule structure.

        Returns
        -------
        str
            Smiles representation of the skeletal molecule structure
        """
        return rdchem.MolToSmiles(self.mol_skeletal)

    @cached_property
    def inchi(self) -> str:
        """InChI representation of the skeletal molecule structure.

        Returns
        -------
        str
            InChI representation of the skeletal molecule structure.
        """
        return rdchem.MolToInchi(self.mol_skeletal)

    # @cached_property
    # def axs(self):

    def calc_all(self):
        """Helper method to allow for precalculation of all cached properties"""
        self.per_state
        self.inter_state
        self.pops
        self.delta_E
        self.fosc_time
        self.spectra
        self.spectra_groups
        self.pca_data
        self.hops
        self.structure_atXYZ
        self.mol_skeletal
        self.smiles
        self.inchi

    def plot_per_state_histograms(
        self,
        state_selection: StateSelection,
        fig: Figure | SubFigure | None = None,
        shape: tuple[int, int] | None = None,
    ) -> dict[str, Axes]:
        """Plot histograms of forces, energies and permanent dipoles for each selected state.

        Parameters
        ----------
        fig : Figure | SubFigure | None, optional
            Figure to plot the graphs to. Defaults to None.
        shape : tuple[int,int], optional
            The shape (rows, cols) that the sub-plots should take. Defaults to one row and 3 columns.

        Returns
        -------
        dict[str, Axes]
            The axes to which the histograms have been plotted with their assigned names.
        """
        start = timer()
        res = plot_per_state_histograms(
            per_state=self.per_state,
            state_selection=state_selection,
            fig=fig,
            shape=shape,
        )
        end = timer()
        info(f"finished plot_per_state_histograms in {end - start} s")
        return res

    def plot_timeplots(
        self, state_selection: StateSelection, fig: Figure | SubFigure | None = None
    ) -> dict[str, Axes]:
        """Create the Time plots of populations and energy level errors of each state for this DataSheetPage.

        Args:
            fig (Figure | SubFigure | None, optional): The figure to plot to. Defaults to None.

        Returns:
            Axes: The axes that have been plotted to
        """
        start = timer()
        res = plot_timeplots(
            pops=self.pops,
            delta_E=self.delta_E,
            fosc_time=self.fosc_time,
            fig=fig,
            state_selection=state_selection,
        )
        end = timer()
        info(f"finished plot_timeplots in {end - start} s")
        return res

    def plot_separated_spectra_and_hists(
        self,
        state_selection: StateSelection,
        fig: Figure | SubFigure | None = None,
        current_multiplicity: int = 1,
    ) -> dict[str, Axes]:
        start = timer()
        res = plot_separated_spectra_and_soc_dip_hists(
            inter_state=self.inter_state,
            spectra_groups=self.spectra_groups,
            fig=fig,
            state_selection=state_selection,
            current_multiplicity=current_multiplicity,
        )
        end = timer()
        info(f"finished plot_separated_spectra_and_hists in {end - start} s")
        return res

    def plot_separated_spectra_and_hists_groundstate(
        self,
        state_selection: StateSelection,
        fig: Figure | SubFigure | None = None,
        scmap=plt.get_cmap('turbo'),
    ) -> dict[str, Axes]:
        start = timer()
        res = plot_separated_spectra_and_soc_dip_hists_groundstate(
            inter_state=self.inter_state,
            spectra_groups=self.spectra_groups,
            fig=fig,
            scmap=scmap,
            state_selection=state_selection,
        )
        end = timer()
        info(
            f"finished plot_separated_spectra_and_hists_groundstate in {end - start} s"
        )
        return res

    def plot_nacs_histograms(
        self,
        state_selection: StateSelection,
        fig: Figure | SubFigure | None = None,
    ) -> dict[str, Axes]:
        if not self.can['nacs_histograms']:
            return {}
        start = timer()
        res = plot_nacs_histograms(
            self.inter_state,
            self.hops,
            fig=fig,
            state_selection=state_selection,
        )
        end = timer()
        info(f"finished plot_nacs_histograms in {end - start} s")
        return res

    def plot_noodle(
        self,
        fig: Figure | SubFigure | None = None,
        state_selection: StateSelection | None = None,
    ) -> Axes:
        if not self.can['noodle']:
            raise RuntimeError(
                f"Cannot plot `noodle plot` on page with name {self.name}"
            )

        start = timer()
        res = plot_noodleplot(self.pca_data, self.hops, fig=fig)
        end = timer()
        info(f"finished plot_noodle in {end - start} s")
        return res

    def plot_energy_bands(
        self,
        state_selection: StateSelection,
        fig: Figure | SubFigure | None = None,
    ) -> dict[str, Axes]:
        start = timer()
        res = plot_energy_bands(
            self.per_state,
            self.pca_data,
            state_selection=state_selection,
            hops_mask=self.hops,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_energy_bands in {end - start} s")
        return res

    def plot_structure(
        self,
        fig: Figure | SubFigure | None = None,
        state_selection: StateSelection | None = None,
    ) -> Axes:
        start = timer()
        mol = self.mol_skeletal if self.structure_skeletal else self.mol
        res = plot_structure(
            mol,
            smiles=self.smiles,
            # inchi=self.inchi,
            ax=None,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_structure in {end - start} s")
        return res

    def plot_pca_structure(
        self,
        fig: Figure | SubFigure | None = None,
        state_selection: StateSelection | None = None,
    ) -> Axes:
        start = timer()
        res = plot_pca_structure(
            self.frames,
            self.pca_data,
            axs=None,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_pca_structure in {end - start} s")
        return res[(0, 'min')]

    def plot_coupling_page(
        self,
        figure: Figure | SubFigure,
        suplots: dict[StateCombination, Axes],
        state_selection: StateSelection,
        simple_mode: bool = False,
    ) -> dict[StateCombination, Axes]:
        """Plot coupling and state-info data on an axes grid.

        Either a plot of all state-coupling information or a matrix of plots with color-coded information of

        Parameters
        ----------
        figure : Figure | SubFigure
            Figure this is being plotted into. Used for some after-the fact manipulation like introducing
        suplots : dict[StateCombination, Axes]
            Axes to plot the data of individual state combinations into.
        state_selection : StateSelection
            The StateSelection object to limit the combinations to include
        simple_mode : bool, optional
            Flag to determine whether we want full coupling plots (NACs, SOCs, dip_trans) or just color plots. Defaults to False meaning permitted and unpermitted transitions will be color-coded.

        Returns
        -------
        dict[StateCombination, Axes]
            The dict of axes after plotting to them
        """
        start = timer()
        # mol = self.mol_skeletal if self.structure_skeletal else self.mol
        interstate = self.inter_state
        has_dip_trans = 'dip_trans_norm' in interstate.data_vars
        has_nacs = 'nacs_norm' in interstate.data_vars
        has_socs = 'socs_norm' in interstate.data_vars
        num_states = len(state_selection.states)
        logging.info(
            f"Rendering coupling page for {num_states} states with at least {len(state_selection.state_combinations)} relevant state transitions."
        )
        if not simple_mode:
            logging.info("This may take a while during saving or rendering.")
        elif isinstance(figure, Figure):
            # Make figure smaller for the simple mode figure.
            figure.set_size_inches(
                # 11.69 / 6 * ncols,
                # 8.27 / 4 * nrows,
                0.5 * num_states,
                0.5 * num_states,
            )  # landscape A4

        res: dict[StateCombination, Axes] = {}

        for sc in state_selection.state_combinations:
            ax1 = suplots[sc]
            res[sc] = ax1

            sccolor = state_selection.get_state_combination_color(sc)
            centertext(
                "",  # "NAC",
                ax1,
                clearticks="xy",
                background_color=sccolor,  # "green"
            )

        for sc in tqdm(state_selection.state_combinations):
            sc_label = state_selection.get_state_combination_tex_label(sc)
            s1, s2 = sc
            sc_r = (s2, s1)
            ax2 = suplots[sc_r]
            res[sc_r] = ax2

            s1label = state_selection.get_state_tex_label(s1)
            s2label = state_selection.get_state_tex_label(s2)
            sccolor = state_selection.get_state_combination_color(sc)
            if has_nacs or has_dip_trans:
                interstate_sc = InterState(
                    direct_interstate_data=interstate.dataset.sel(statecomb=sc)
                )
                if has_nacs:
                    if interstate_sc.nacs_norm.max() > 1e-9:
                        if simple_mode:
                            centertext(
                                r"$\checkmark$",  # "NAC",
                                ax2,
                                clearticks="xy",
                                background_color="green",
                            )
                        else:
                            single_dip_trans_hist(
                                interstate_sc,
                                sc_label,
                                (s2label, s1label),
                                sccolor,
                                ax=ax2,
                                plot_marginals=False,
                            )
                        continue
                if has_dip_trans:
                    if interstate_sc.dipole_transition_norm.max() > 1e-9:
                        centertext(
                            r"$\checkmark$",
                            # r"$\mathbf{\mu}_\mathrm{trans}$",
                            ax2,
                            clearticks="xy",
                            background_color="green",
                        )
                        continue

                if has_socs:
                    found_soc = False
                    interstate_sc = InterState(
                        direct_interstate_data=interstate.dataset.sel(
                            statecomb=sc, full_statecomb=sc_r
                        )
                    )
                    if interstate_sc.socs_norm.max() > 1e-9:
                        if simple_mode:
                            centertext(
                                "X",
                                # r"$SOC$",
                                ax2,
                                clearticks="xy",
                                background_color="brown",
                            )
                        else:
                            single_soc_trans_hist(
                                interstate_sc,
                                sc_label,
                                (s2label, s1label),
                                sccolor,
                                ax=ax2,
                                plot_marginals=False,
                            )
                        found_soc = True

                    if found_soc:
                        continue

            centertext(
                r"?",  # f"No Permitted Transition {sc_label}",
                ax2,
                clearticks="xy",
            )

        for state in state_selection.states:
            if simple_mode:
                centertext(
                    r"",
                    suplots[(state, state)],
                    clearticks="xy",
                    background_color=state_selection.get_state_color(state),  # "brown"
                )
            else:
                suplots[(state, state)].remove()

        if figure is not None:
            state_labels = [
                f"${state_selection.get_state_tex_label(s)}$"
                for s in state_selection.states
            ]
            label_plot_grid(figure, row_headers=state_labels, col_headers=state_labels)
        end = timer()
        info(f"finished plot_coupling_page in {end - start} s")
        return res

    def render_meta_page(
        self,
    ) -> tuple[Figure, dict[str, SubFigure]]:
        """Helper function to output the entire `meta` overview page to a Shnitsel Tools datasheet.

        Returns
        -------
        tuple[Figure, dict[str, SubFigure]]
            The pair of the entire figure and a dict of the subfigures involved in the meta-overview page.
        """
        letter_base = "abcdefghijk"
        letter_it = iter(letter_base)
        fig, subfigs = self.get_subfigures_meta_page()
        fig.suptitle(f'Datasheet:{self.name} [Page: Overview]', fontsize=16)

        ax = plot_structure(
            self.mol,
            ax=None,
            fig=subfigs['structure_plot'],
        )
        outlabel(ax, next(letter_it))

        metainfo = []
        if 'trajectory' in self.frames.sizes:
            num_trajs = self.frames.sizes['trajectory']
            trajectory_ids = list(set(self.frames.coords['trajectory'].values))
        elif 'atrajectory' in self.frames.coords:
            trajectory_ids = list(set(self.frames.coords['atrajectory'].values))

            num_trajs = len(trajectory_ids)
        else:
            num_trajs = 1
            trajectory_ids = ['1']

        time_unit = 'unknown'
        if self.frames.has_coordinate('time'):
            time_unit = self.frames.coords['time'].attrs['units']
        # metainfo.append(('$t$ unit', time_unit))

        def var_or_attr(ds, name, default=None):
            if name in self.frames.dataset:
                return self.frames.dataset[name].values
            elif name in self.frames.attrs:
                return self.frames.dataset.attrs[name]
            else:
                return default

        t_max = np.max(var_or_attr(self.frames, 't_max', -1))
        metainfo.append(('Source software', self.frames.attrs['input_format']))
        metainfo.append(
            ('Source software version', self.frames.attrs['input_format_version'])
        )
        metainfo.append(('Data type (dyn/stat)', self.frames.attrs['input_type']))
        metainfo.append(
            ('Basis set', self.frames.attrs.get('theory_basis_set', 'unknown'))
        )
        metainfo.append(('EST level', self.frames.attrs.get('est_level', 'unknown')))
        metainfo.append(('Compound [smile]', self.smiles))
        # metainfo.append(('Compound [InChI]', format_inchi(self.inchi)))

        metainfo.append(('Number of trajectories', str(num_trajs)))
        if 'frame' in self.frames.sizes:
            metainfo.append(('Number of frames', str(self.frames.sizes['frame'])))
        elif 'time' in self.frames.sizes:
            metainfo.append(('Number of frames', str(self.frames.sizes['time'])))
        metainfo.append(('Maximum $t$', np.max(var_or_attr(self.frames, 't_max', -1))))
        metainfo.append(
            (
                'Timestep $\\Delta t$',
                np.max(var_or_attr(self.frames, 'delta_t', -1)),
            )
        )
        metainfo.append(('$t$ unit', time_unit))
        metainfo.append(
            (
                'Forces in set',
                (
                    'all'
                    if self.frames.forces_format == True
                    else self.frames.forces_format
                ),
            )
        )
        metainfo.append(('Num Singlets', self.frames.num_singlets))

        if self.frames.num_doublets > 0:
            metainfo.append(('Num Doublets', self.frames.num_doublets))

        metainfo.append(('Num Triplets', self.frames.num_triplets))

        var_meta_info: list[tuple[str, tuple[str, str]]] = []
        for varname, val in self.frames.data_vars.items():
            unit = '-'
            orig_unit = '-'

            if 'units' in val.attrs:
                unit = str(val.attrs['units'])
            if 'original_units' in val.attrs:
                orig_unit = str(val.attrs['original_units'])

            var_meta_info.append((varname, (unit, orig_unit)))
        for varname, val in self.frames.coords.items():
            unit = '-'
            orig_unit = '-'

            if 'units' in val.attrs:
                unit = str(val.attrs['units'])
            else:
                continue
            if 'original_units' in val.attrs:
                orig_unit = str(val.attrs['original_units'])

            var_meta_info.append((str(varname), (unit, orig_unit)))

        ax = subfigs['length_plot'].add_subplot(1, 1, 1)
        if num_trajs > 1:
            ax.set_ylabel('Trajectory (sorted by length)')
            ax.set_xlabel(f'Traj. length (t / {time_unit})')
            ax.set_ylim((0, num_trajs - 1))
            ax.set_xlim((0, t_max))
            from matplotlib.ticker import MaxNLocator

            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

            cutoff_times = []
            if 'atrajectory' in self.frames.dataset and 'time' in self.frames.coords:
                for id, traj_data in self.frames.dataset.groupby('atrajectory'):
                    t_max_present = traj_data['time'].max()
                    cutoff_times.append(t_max_present)
            cutoff_times.sort()
            index = list(range(num_trajs))
            ax.plot(cutoff_times, index)

            ax.fill_between(
                [0] + cutoff_times + [t_max],
                [0] + index + [num_trajs - 1],
                [num_trajs - 1] + [num_trajs - 1] * len(cutoff_times) + [num_trajs - 1],
                color=st_grey,
            )

        else:
            centertext("Not enough data", ax, clearticks='xy')
        outlabel(ax, next(letter_it))
        ax = subfigs['metadata_a'].add_subplot(1, 1, 1)
        ax.axis('off')
        meta_table = ax.table(
            [[k, str(v)] for k, v in metainfo],
            colLabels=['Attribute', 'Value'],
            loc='center',
        )
        meta_table.auto_set_font_size(False)
        meta_table.set_fontsize(12)
        # centertext("Metadata", ax, clearticks='xy')
        outlabel(ax, next(letter_it))
        # Variables and units
        ax = subfigs['variables_a'].add_subplot(1, 1, 1)
        ax.axis('off')
        var_table = ax.table(
            [[k, *(v)] for k, v in var_meta_info],
            colLabels=['Variable', 'Unit', 'Orig unit'],
            loc='center',
        )
        var_table.auto_set_font_size(False)
        var_table.set_fontsize(12)
        # centertext("Metadata", ax, clearticks='xy')
        outlabel(ax, next(letter_it))
        hist_axs = self.plot_per_state_histograms(
            self.state_selection, subfigs['per_state_histograms'], shape=(3, 1)
        )
        # print(hist_axs.keys())
        outlabel(hist_axs['energy'], next(letter_it))

        return fig, subfigs

    @staticmethod
    def get_subfigures_main_page(
        include_per_state_hist: bool = False, borders: bool = False
    ) -> tuple[Figure, dict[str, SubFigure]]:
        """Helper function to prepare a figure to hold all subfigures in this DatasheetPage

        Parameters
        ----------
        include_per_state_hist : bool, optional
            Flag whether per state histograms will be included. Defaults to False.
        borders : bool, optional
            Flag whether figure borders should be drawn. Defaults to False.

        Returns
        -------
        tuple[Figure, dict[str, SubFigure]]
            The overall figure and a dict to access individual subfigures by their name.
        """
        nrows = 6 if include_per_state_hist else 5
        top_spacing = 1 if include_per_state_hist else 0

        fig, oaxs = plt.subplots(nrows, 3, layout='constrained')
        vscale = 1 if include_per_state_hist else 5 / 6
        fig.set_size_inches(8.27, 11.69 * vscale)  # portrait A4
        if borders:
            fig.set_facecolor('#0d0d0d')
        gs = oaxs[0, 0].get_subplotspec().get_gridspec()
        for ax in oaxs.ravel():
            ax.remove()
        gridspecs = dict(
            per_state_histograms=gs[0, :],
            timeplots=gs[top_spacing + 2 :, 2],
            noodle=gs[top_spacing + 0 : top_spacing + 2, 1:],
            separated_spectra_and_hists=gs[top_spacing + 0 :, 0],
            nacs_histograms=gs[top_spacing + 3 :, 1],
            structure=gs[top_spacing + 2, 1],
        )
        if not include_per_state_hist:
            del gridspecs['per_state_histograms']
        subfigures = {
            sub_name: fig.add_subfigure(sub_gridspec)
            for sub_name, sub_gridspec in gridspecs.items()
        }
        return fig, subfigures

    @staticmethod
    def get_subfigures_pca_page(
        borders: bool = False,
    ) -> tuple[Figure, dict[str, SubFigure]]:
        """Helper function to prepare a figure to hold all subfigures in this DatasheetPage covering all PCA information.

        Parameters
        ----------
        borders : bool, default=False
            Flag whether figure borders should be drawn. Defaults to False.

        Returns
        -------
        tuple[Figure, dict[str, SubFigure]]
            The overall figure and a dict to access individual subfigures by their name.
        """
        nrows = 8

        fig, oaxs = plt.subplots(nrows, 4, layout='constrained')
        fig.set_size_inches(8.27, 11.69)  # portrait A4

        if borders:
            fig.set_facecolor('#0d0d0d')

        gs = oaxs[0, 0].get_subplotspec().get_gridspec()

        for ax in oaxs.ravel():
            ax.remove()
        gridspecs = dict(
            pca_plot=gs[:2, :2],
            pca_extrema_plot=gs[:2, 2:],
            feature_selection=gs[2:4, :],
            feature_explanation=gs[4:, :],
        )
        subfigures = {
            sub_name: fig.add_subfigure(sub_gridspec)
            for sub_name, sub_gridspec in gridspecs.items()
        }
        return fig, subfigures

    @staticmethod
    def get_subfigures_meta_page(
        borders: bool = False,
    ) -> tuple[Figure, dict[str, SubFigure]]:
        """Helper function to prepare a figure to hold all subfigures in this DatasheetPage covering all Meta information.

        Parameters
        ----------
        borders : bool, default=False
            Flag whether figure borders should be drawn. Defaults to False.

        Returns
        -------
        tuple[Figure, dict[str, SubFigure]]
            The overall figure and a dict to access individual subfigures by their name.
        """
        nrows = 3

        fig, oaxs = plt.subplots(nrows, 6, layout='constrained')
        fig.set_size_inches(8.27, 11.69)  # portrait A4

        if borders:
            fig.set_facecolor('#0d0d0d')

        gs = oaxs[0, 0].get_subplotspec().get_gridspec()

        for ax in oaxs.ravel():
            ax.remove()

        gridspecs = dict(
            structure_plot=gs[0, :3],
            length_plot=gs[0, 3:],
            metadata_a=gs[1, :4],
            # metadata_b=gs[1, 1],
            variables_a=gs[2, :4],
            # variables_b=gs[2, 1],
            per_state_histograms=gs[1:, 4:],
        )
        subfigures = {
            sub_name: fig.add_subfigure(sub_gridspec)
            for sub_name, sub_gridspec in gridspecs.items()
        }
        return fig, subfigures

    @staticmethod
    def get_subfigures_coupling_page(
        state_selection: StateSelection, borders: bool = False
    ) -> tuple[Figure, dict[StateCombination, Axes]]:
        """Helper function to prepare a figure to hold all state interaction figures.

        Parameters
        ----------
        n_states : int, optional
            Number of states (will be square in the end.)
        borders : bool, optional
            Flag whether figure borders should be drawn. Defaults to False.

        Returns
        ----------
        tuple[Figure, dict[str, SubFigure]]
            The overall figure and a dict to access individual subfigures by their name.
        """
        n_states: int = len(state_selection.states)
        nrows = n_states
        ncols = n_states

        states = state_selection.states

        fig, state_grid = plt.subplots(nrows, ncols, layout='constrained')
        fig.set_size_inches(
            # 11.69 / 6 * ncols,
            # 8.27 / 4 * nrows,
            2 * ncols,
            2 * nrows,
        )  # landscape A4
        if borders:
            fig.set_facecolor('#0d0d0d')
        # gs = state_grid[0, 0].get_subplotspec().get_gridspec()
        # gridspecs = {
        #     (states[a], states[b]): gs[a, b]
        #     for a in range(n_states)
        #     for b in range(n_states)
        # }
        # subfigures = {
        #     comb_label: fig.add_subfigure(sub_gridspec)
        #     for comb_label, sub_gridspec in gridspecs.items()
        # }

        subplots = {
            (states[a], states[b]): state_grid[a, b]
            for a in range(n_states)
            for b in range(n_states)
        }
        return fig, subplots

    def plot(
        self,
        include_per_state_hist: bool = False,
        include_coupling_page: bool = True,
        include_pca_page: bool = False,
        include_meta_page: bool = False,
        borders: bool = False,
        consistent_lettering: bool = True,
    ) -> Figure | list[Figure]:
        """Function to plot this Datasheet.

        Will generate all subplots and calculate necessary data if it has not yet been generated.

        Parameters
        ----------
        include_per_state_hist : bool, optional
            Flag whether per-state histograms should be included. Defaults to False.
        include_coupling_page : bool, optional
            Flag to create a full page with state-coupling plots. Defaults to True.
        include_pca_page : bool, optional
            Flag to create a PCA analysis page with details on PCA results. Defaults to False.
        include_meta_page : bool, optional
            Flag to add a page with meta-information about the trajectory data. Defaults to False
        borders : bool, optional
            Flag whether the figure should have borders or not. Defaults to False.
        consistent_lettering : bool, optional
            Flag whether consistent lettering should be used, i.e. whether the same plot should always have the same label letter. Defaults to True.

        Returns
        -------
        Figure
            The figure holding the entirety of plots in this Datasheet page.
        list[Figure]
            The list of figures holding the plots in this Datasheet page if there are multiple (depending on the configuration flags)
        """
        letter_base = 'abcdefghijkl'

        figures = []

        pages = []
        if self.state_selection_provided:
            pages = [(1, self.state_selection, "Selection")]

        pages = pages + [
            (
                1,
                self.state_selection.select_states(
                    multiplicity=1, min_states_in_selection=2
                ),
                "Singlets",
            ),
            (
                2,
                self.state_selection.select_states(
                    multiplicity=2, min_states_in_selection=1
                ),
                "Doublets",
            ),
            (
                3,
                self.state_selection.select_states(
                    multiplicity=3, min_states_in_selection=1
                ),
                "Triplets",
            ),
        ]

        pages_handled = []

        for page_index, (page_mult, page_selection, page_title) in enumerate(pages):
            if not page_selection.states:
                continue

            # print(mult, mult_name, page_selection.states)
            # print(mult, mult_name, page_selection.state_combinations)
            # continue

            pages_handled.append(page_index)

            fig, sfs = self.get_subfigures_main_page(
                include_per_state_hist=include_per_state_hist, borders=borders
            )

            letter_it = iter(letter_base)
            fig.suptitle(f'Datasheet:{self.name} [Page:{page_title}]', fontsize=16)

            # print(self.frames)

            # separated_spectra_and_hists
            if self.can['separated_spectra_and_hists']:
                axs = self.plot_separated_spectra_and_hists(
                    state_selection=page_selection,
                    fig=sfs['separated_spectra_and_hists'],
                    current_multiplicity=page_mult,
                )
                ax = axs['sg']
                outlabel(ax, next(letter_it))
            else:
                ax = sfs['separated_spectra_and_hists'].subplots(1, 1)
                centertext(r"No $\mathbf{\mu}_{ij}$ data", ax=ax)
                ax.get_yaxis().set_visible(False)
                ax.get_xaxis().set_visible(False)
                inlabel(ax, next(letter_it))
            # noodle
            if self.can['noodle']:
                if len(pages_handled) == 1:
                    ax = self.plot_noodle(
                        state_selection=page_selection, fig=sfs['noodle']
                    )
                    inlabel(ax, next(letter_it))
                else:
                    axs = self.plot_energy_bands(
                        state_selection=page_selection, fig=sfs['noodle']
                    )
                    outlabel(axs['pc1'], next(letter_it))
            elif consistent_lettering:
                next(letter_it, next(letter_it))
            # structure
            if self.can['structure']:
                if len(pages_handled) == 1 or not self.can['noodle']:
                    ax = self.plot_structure(
                        state_selection=page_selection, fig=sfs['structure']
                    )
                    outlabel(ax, next(letter_it))
                else:
                    if consistent_lettering:
                        next(letter_it)
                    # ax = self.plot_pca_structure(
                    #     state_selection=page_selection, fig=sfs['structure']
                    # )

            elif consistent_lettering:
                next(letter_it)
            # nacs_histograms
            if self.can['nacs_histograms']:
                axs = self.plot_nacs_histograms(
                    state_selection=page_selection, fig=sfs['nacs_histograms']
                )
                ax = axs.get('ntd', axs['nde'])
                outlabel(ax, next(letter_it))
            elif consistent_lettering:
                next(letter_it)
            # time plots
            if self.can['timeplots']:
                axs = self.plot_timeplots(
                    state_selection=page_selection.select_states(
                        multiplicity=page_mult, min_states_in_selection=2
                    ),
                    fig=sfs['timeplots'],
                )
                ax = axs['pop']
                outlabel(ax, next(letter_it))
            elif consistent_lettering:
                next(letter_it)
            if include_per_state_hist:
                axs = self.plot_per_state_histograms(
                    state_selection=page_selection, fig=sfs['per_state_histograms']
                )
                ax = axs['energy']
                outlabel(ax, next(letter_it))
            elif consistent_lettering:
                next(letter_it)

            figures.append(fig)

        if include_coupling_page:
            fig, sfs = self.get_subfigures_coupling_page(
                state_selection=self.state_selection, borders=borders
            )
            fig.suptitle(
                f'Datasheet:{self.name} [Page: Inter-state Couplings]', fontsize=16
            )
            res_coupling_axes = self.plot_coupling_page(
                fig, sfs, self.state_selection, simple_mode=True
            )
            figures.append(fig)

        if include_pca_page:
            fig, subfigs = self.get_subfigures_pca_page()

            fig.suptitle(f'Datasheet:{self.name} [Page: PCA]', fontsize=16)
            ax = self.plot_noodle(
                state_selection=self.state_selection, fig=subfigs['pca_plot']
            )
            ax = self.plot_pca_structure(
                state_selection=self.state_selection, fig=subfigs['pca_extrema_plot']
            )
            ax = subfigs['feature_selection'].subplots(1, 1)
            centertext("Missing", ax=ax)
            ax = subfigs['feature_explanation'].subplots(1, 1)
            centertext("Missing", ax=ax)
            figures.append(fig)
            self.pca_explanation

        if include_meta_page:
            fig, subfigs = self.render_meta_page()
            figures.append(fig)

        # TODO: FIXME: Add Spearman's rank correlation coefficient analysis of energy.
        # TODO: FIXME: Full PCA with internal coordinates?
        # TODO: FIXME: Explain parts of PCA components?
        # TODO: FIXME: Electron bond matrix?
        # TODO: FIXME: State-shift changes in Spearman's analysis? Identify most important parts of molecule that change?

        return figures if len(figures) != 1 else figures[0]

    def _test_subfigures(
        self, include_per_state_hist: bool = False, borders: bool = False
    ):
        """Helper method to test whether subfigures are successfully plotted

        Parameters
        ----------
        include_per_state_hist : bool, default=False
            Flag to include per-state histograms. Defaults to False.
        borders : bool, default=False
            Whether the figures should have borders. Defaults to False.
        """
        fig, sfs = self.get_subfigures_main_page(
            include_per_state_hist=include_per_state_hist, borders=borders
        )
        for sf in sfs.values():
            sf.subplots(2, 2)
        if include_per_state_hist:
            sfs['per_state_histograms'].set_facecolor('blue')
        sfs['timeplots'].set_facecolor('green')
        sfs['separated_spectra_and_hists'].set_facecolor('orange')
        sfs['nacs_histograms'].set_facecolor('yellow')
        sfs['noodle'].set_facecolor('red')
        sfs['structure'].set_facecolor('purple')
