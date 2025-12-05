from dataclasses import dataclass
from functools import cached_property
from matplotlib.axes import Axes
import xarray as xr
import numpy as np
import rdkit.Chem as rdchem
import matplotlib as mpl
import matplotlib.pyplot as plt
from logging import info, warning
from timeit import default_timer as timer

from matplotlib.figure import Figure, SubFigure

import shnitsel
from shnitsel.analyze.populations import calc_classical_populations
import shnitsel.bridges
from shnitsel.analyze import stats
from shnitsel.core.typedefs import AtXYZ, InterState, PerState, SpectraDictType
from shnitsel.data.trajectory_format import Trajectory
from shnitsel.filtering.state_selection import StateSelection

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from shnitsel.analyze.spectra import assign_fosc, calc_spectra, get_spectra_groups

from .figures.common import centertext
from .figures.per_state_hist import plot_per_state_histograms
from .figures.time import plot_timeplots
from .figures.dip_trans_hist import (
    plot_separated_spectra_and_hists,
    plot_separated_spectra_and_hists_groundstate,
)
from .figures.nacs_hist import plot_nacs_histograms
from ..plot.pca_biplot import plot_noodleplot
from .figures.structure import plot_structure


@dataclass
class DatasheetPage:
    """Class encapsulating the data generation and plotting of a single page of a datasheet.

    Involves plotting of inter-state graphs like spectra, statistics of energy deltas and transitional dipoles,
    histograms of per-state and inter-state properties, etc.
    """

    state_selection: StateSelection
    feature_selection: None = None
    spectra_times: list[int | float] | np.ndarray | None = None
    charge: int = 0
    structure_skeletal: bool = False
    name: str = ''

    def __init__(
        self,
        data: Trajectory | Self,
        state_selection: StateSelection | None = None,
        feature_selection: None = None,
        *,
        spectra_times: list[int | float] | np.ndarray | None = None,
        col_state: list | None = None,
        col_inter: list | None = None,
    ):
        """Create a new DataSheet page, initializating the data of the class and preparing the plotting later on

        Args:
            data (Trajectory | Self): Either a Shnitsel Trajectory data object from which to generate the statistical information visualized in this DataSheet page
                or another DataSheetPage that should be copied.
            state_selection (StateSelection, optional): Optional parameter to specify a subset of states and state combinations that may be considered for the dataset.
                Will be generated if not provided.
            # TODO: FIXME: Add feature selection option
            feature_selection (optional): Optional parameter to limit the PCA plot and analysis to a specific subset of the structure. Will be generated if not provided.
            spectra_times (list[int  |  float] | np.ndarray | None, optional): Sequence of times to calculate spectra at. Defaults to None.
            col_state (list | None, optional): A list of colors to use for the states. Defaults to default shnitsel colors.
            col_inter (list | None, optional): A list of colors to use for state combinations. Defaults to default shnitsel colors.

        Raises:
            TypeError: If wrong type of `data` parameter is provided.
            ValueError: If the wrong number of colors for the states is provided.
            ValueError: If the wrong number of colors for state transitions is provided.

        Returns:
            DatasheetPage: The constructed (or copied) DataSheetPage
        """
        if isinstance(data, DatasheetPage):
            self._copy_data(old=data)
            return
        elif isinstance(data, Trajectory):
            self.frames = data
        else:
            raise TypeError("Neither DatasheetPage nor frames/Trajectory given.")

        assert isinstance(self.frames, xr.Dataset)

        if spectra_times is not None:
            self.spectra_times = spectra_times
        elif 'time' not in self.frames:
            warning("No 'time' variable found. Have ICONDs been passed as frames?")
        elif self.frames is not None:
            max_time = self.frames.coords['time'].max().item()
            self.spectra_times = [max_time * i / 40 for i in range(5)]
            self.spectra_times += [max_time * i / 20 for i in range(5)]
            self.spectra_times += [max_time * i / 3 for i in range(4)]

        # Initialize state selection or use provided selection
        if state_selection is not None:
            self.state_selection = state_selection
        else:
            self.state_selection = StateSelection.init_from_dataset(self.frames)

        # Initialize feature selection or use provided selection
        if feature_selection is not None:
            self.feature_selection = feature_selection
        else:
            self.feature_selection = None

        # print(self.frames)

        nstates = self.frames.sizes['state']
        if col_state is not None:
            assert (ncols := len(col_state)) == nstates, (
                f"`col_state` has {ncols} colors, "
                f"but should contain one color for each of the {nstates} states"
            )
            self.col_state = col_state
        elif nstates <= 3:
            # SHNITSEL-colours
            self.col_state = ['#4DAD15', '#AD2915', '#7515AD'][:nstates]
        elif nstates <= 10:
            cmap = plt.get_cmap('tab10')
            self.col_state = [mpl.colors.rgb2hex(c) for c in cmap.colors][:nstates]  # type: ignore
        elif nstates <= 20:
            cmap = plt.get_cmap('tab20')
            self.col_state = [mpl.colors.rgb2hex(c) for c in cmap.colors][:nstates]  # type: ignore
        else:
            raise ValueError(
                f"These data have {nstates} states. "
                "When passing data with more than 10 states, please "
                "also pass an appropriate colormap to `col_state`."
            )

        ncombs = self.frames.sizes['statecomb']
        if col_inter is not None:
            assert (ncols := len(col_inter)) == ncombs, (
                f"`col_inter` has {ncols} colors, "
                f"but should contain one color for each of the {ncombs} state combinations"
            )
            self.col_inter = col_inter
        elif ncombs <= 3:
            self.col_inter = col_inter or ['#2c3e50', '#C4A000', '#7E5273'][:ncombs]
        elif ncombs <= 10:
            # TODO: choose colours distinct from per_state colours
            cmap = plt.get_cmap('tab10')
            self.col_inter = [mpl.colors.rgb2hex(c) for c in cmap.colors][:ncombs]  # type: ignore
        elif ncombs <= 20:
            cmap = plt.get_cmap('tab20')
            self.col_inter = [mpl.colors.rgb2hex(c) for c in cmap.colors][:ncombs]  # type: ignore
        else:
            raise ValueError(
                f"These data have {ncombs} state combinations. "
                "When passing data with more than 10 state combinations, please "
                "also pass an appropriate colormap to `col_inter`."
            )

        self.can = {}

        # print(self.frames['state_charges'])
        if 'state_charges' in self.frames:
            self.charge = int(self.frames.state_charges.isel(state=0))

        def check(*ks):
            return all(k in self.frames for k in ks)

        self.can['per_state_histograms'] = check('energy', 'forces', 'dip_trans')
        self.can['separated_spectra_and_hists'] = check('dip_trans', 'time')
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

        Args:
            old (Self): The DataSheet to create a copy of.
        """
        self.spectra_times = old.spectra_times
        self.state_selection = old.state_selection
        self.feature_selection = old.feature_selection
        self.col_state = old.col_state
        self.col_inter = old.col_inter
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

        Returns:
            PerState: Per-state data of the self.frames object. (Energies, permanent dipoles)
        """
        start = timer()
        per_state = stats.get_per_state(self.frames)
        per_state['_color'] = 'state', self.col_state
        end = timer()
        info(f"cached per_state in {end - start} s")
        return per_state

    @cached_property
    def inter_state(self) -> InterState:
        """Inter-state (state-transition) data of the underlying self.frames object

        Returns:
            InterState: Inter-state properties of the underlying data. (delta_energies, transition dipoles, SOCs, NACs, fosc)
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        inter_state = stats.get_inter_state(self.frames)
        inter_state['_color'] = 'statecomb', self.col_inter

        # Calculate fosc if missing and conditions met
        if (
            "fosc" not in inter_state
            and 'dip_trans' in inter_state
            and "energy_interstate" in inter_state
        ):
            inter_state = assign_fosc(inter_state)

        for var, tex in [
            ('energy', r"$\Delta E$"),
            ('nacs', r"$\|\mathrm{NAC}_{i,j}\|_2$"),
            ('dip_trans', r"$\|\mathbf{\mu}_{i,j}\|_2$"),
            ('fosc', r"$f_\mathrm{osc}$"),
        ]:
            try:
                inter_state[var].attrs['tex'] = tex
            except KeyError:
                pass
        end = timer()
        info(f"cached inter_state in {end - start} s")
        return inter_state

    @cached_property
    def pops(self) -> xr.DataArray:
        """Population data for the underlying self.frames

        Returns:
            xr.DataArray: Population data encapsulated in a DataArray.
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        pops = calc_classical_populations(self.frames)
        pops['_color'] = 'state', self.col_state
        end = timer()
        info(f"cached pops in {end - start} s")
        return pops

    @cached_property
    def delta_E(self) -> xr.Dataset:
        """Energy deltas between different states in a dataset.

        Returns:
            xr.Dataset: Dataset holding 'energy_interstate' variable.
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        res = stats.time_grouped_confidence_interval(
            self.inter_state['energy_interstate']
        )
        res['_color'] = 'statecomb', self.col_inter
        res.attrs['tex'] = r"$\Delta E$"
        end = timer()
        info(f"cached delta_E in {end - start} s")
        return res

    @cached_property
    def fosc_time(self) -> xr.Dataset | None:
        """Strength of oscillator/transition rate between states at different points in time.

        Returns:
            xr.Dataset | None: Either the f_osc data (with confidence intervals) or None if not sufficient data in self.frames to calculate it.
        """
        start = timer()
        if 'fosc' in self.inter_state:
            res = stats.time_grouped_confidence_interval(self.inter_state['fosc'])
            res['_color'] = 'statecomb', self.col_inter
            res.attrs['tex'] = r"$f_\mathrm{osc}$"
        else:
            res = None
        end = timer()
        info(f"cached fosc_time in {end - start} s")
        return res

    @cached_property
    def spectra(self) -> SpectraDictType:
        """Spectral statistics of the self.frames object.

        Returns:
            SpectraDictType: The spectral information per state transition
        """
        start = timer()
        # TODO: FIXME: Use state selection for limit on which to calculate
        res = calc_spectra(self.inter_state, times=self.spectra_times)
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

        Returns:
            tuple[ SpectraDictType, SpectraDictType, ]: One spectral dict per ground-state or excited-state transitions.
        """
        start = timer()
        # TODO: FIXME: Use state selection for split
        res = get_spectra_groups(self.spectra)
        end = timer()
        info(f"cached spectra_groups in {end - start} s")
        return res

    @cached_property
    def spectra_ground(self) -> SpectraDictType:
        """Extracted spectral information of only the ground-state transitions

        Returns:
            SpectraDictType: Extracted spectral information of only the ground-state transitions
        """
        return self.spectra_groups[0]

    @cached_property
    def spectra_excited(self) -> SpectraDictType:
        """Extracted spectral information of only the excited-state transitions

        Returns:
            SpectraDictType: Extracted spectral information of only the excited-state transitions
        """
        return self.spectra_groups[1]

    @cached_property
    def noodle(self) -> xr.DataArray:
        """Noodle plot source data derived from principal component analysis (PCA) on the full data in self.frames using only pairwise distances.

        Returns:
            xr.DataArray: The pairwise distance PCA results
        """
        from shnitsel.analyze.pca import pairwise_dists_pca

        start = timer()
        res = pairwise_dists_pca(self.frames.atXYZ)
        end = timer()
        info(f"cached noodle in {end - start} s")
        return res

    @cached_property
    def hops(self) -> xr.DataArray:
        """The PCA plots at the hopping points

        Returns:
            xr.DataArray: PCA data at the hopping points
        """
        from shnitsel.data.multi_indices import mdiff

        mask = mdiff(self.frames.astate) != 0
        return self.noodle[mask]

    @cached_property
    def structure_atXYZ(self) -> AtXYZ:
        """Structure/Position data in the first frame/timestep of the trajectory

        Returns:
            AtXYZ: Positional data.
        """
        if "frame" in self.frames.sizes:
            return self.frames.atXYZ.isel(frame=0)
        else:
            return self.frames.atXYZ.isel(time=0)

    @cached_property
    def mol(self) -> rdchem.Mol:
        """Property to get an rdkit Mol object from the structural data

        Returns:
            rdkit.Chem.Mol: Molecule object representing the structure in the first frame
        """
        # TODO: FIXME: Shouldn't this be a private attribute prefixed with `__` ?
        if 'smiles_map' in self.frames['atXYZ'].attrs:
            mol = shnitsel.bridges.numbered_smiles_to_mol(
                self.frames['atXYZ'].attrs['smiles_map']
            )
            for atom in mol.GetAtoms():
                atom.ClearProp("molAtomMapNumber")
                atom.SetProp("atomNote", str(atom.GetIdx()))
            return mol
        else:
            return shnitsel.bridges.to_mol(self.structure_atXYZ, charge=self.charge)

    @cached_property
    def mol_skeletal(self) -> rdchem.Mol:
        """Skeletal representation of the the rdkit.Chem.Mol representation of the structure

        Returns:
            rdkit.Chem.Mol: Molecule object representing the skeletal structure (no H atoms) in the first frame
        """
        mol = rdchem.Mol(self.mol)
        return rdchem.RemoveHs(mol)

    @cached_property
    def smiles(self) -> str:
        """Smiles representation of the skeletal molecule structure.

        Returns:
            str: Smiles representation of the skeletal molecule structure
        """
        return rdchem.MolToSmiles(self.mol_skeletal)

    @cached_property
    def inchi(self) -> str:
        """InChI representation of the skeletal molecule structure.

        Returns:
            str: InChI representation of the skeletal molecule structure.
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
        self.noodle
        self.hops
        self.structure_atXYZ
        self.mol_skeletal
        self.smiles
        self.inchi

    def plot_per_state_histograms(
        self, fig: Figure | SubFigure | None = None
    ) -> dict[str, Axes]:
        """Plot histograms of forces, energies and permanent dipoles for each selected state.

        Args:
            fig (Figure | SubFigure | None, optional): Figure to plot the graphs to. Defaults to None.

        Returns:
            Axes: _description_
        """
        start = timer()
        res = plot_per_state_histograms(
            per_state=self.per_state,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_per_state_histograms in {end - start} s")
        return res

    def plot_timeplots(self, fig: Figure | SubFigure | None = None) -> dict[str, Axes]:
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
            state_selection=self.state_selection,
        )
        end = timer()
        info(f"finished plot_timeplots in {end - start} s")
        return res

    def plot_separated_spectra_and_hists(
        self, fig: Figure | SubFigure | None = None
    ) -> dict[str, Axes]:
        start = timer()
        res = plot_separated_spectra_and_hists(
            inter_state=self.inter_state,
            sgroups=self.spectra_groups,
            fig=fig,
            state_selection=self.state_selection,
        )
        end = timer()
        info(f"finished plot_separated_spectra_and_hists in {end - start} s")
        return res

    def plot_separated_spectra_and_hists_groundstate(
        self, fig: Figure | SubFigure | None = None, scmap=plt.get_cmap('turbo')
    ) -> dict[str, Axes]:
        start = timer()
        res = plot_separated_spectra_and_hists_groundstate(
            inter_state=self.inter_state,
            spectra_groups=self.spectra_groups,
            fig=fig,
            scmap=scmap,
        )
        end = timer()
        info(
            f"finished plot_separated_spectra_and_hists_groundstate in {end - start} s"
        )
        return res

    def plot_nacs_histograms(
        self, fig: Figure | SubFigure | None = None
    ) -> dict[str, Axes]:
        start = timer()
        res = plot_nacs_histograms(self.inter_state, self.hops.frame, fig=fig)
        end = timer()
        info(f"finished plot_nacs_histograms in {end - start} s")
        return res

    def plot_noodle(self, fig: Figure | SubFigure | None = None) -> Axes:
        start = timer()
        res = plot_noodleplot(self.noodle, self.hops, fig=fig)
        end = timer()
        info(f"finished plot_noodle in {end - start} s")
        return res

    def plot_structure(self, fig: Figure | SubFigure | None = None) -> Axes:
        start = timer()
        mol = self.mol_skeletal if self.structure_skeletal else self.mol
        res = plot_structure(
            mol,
            name=self.name,
            smiles=self.smiles,
            inchi=self.inchi,
            ax=None,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_structure in {end - start} s")
        return res

    @staticmethod
    def get_subfigures(
        include_per_state_hist: bool = False, borders: bool = False
    ) -> tuple[Figure, dict[str, SubFigure]]:
        """Helper function to prepare a figure to hold all subfigures in this DatasheetPage

        Args:
            include_per_state_hist (bool, optional): Flag whether per state histograms will be included. Defaults to False.
            borders (bool, optional): Flag whether figure borders should be drawn. Defaults to False.

        Returns:
            tuple[Figure, dict[str, SubFigure]]: The overall figure and a dict to access individual subfigures by their name.
        """
        nrows = 6 if include_per_state_hist else 5
        s = 1 if include_per_state_hist else 0

        fig, oaxs = plt.subplots(nrows, 3, layout='constrained')
        vscale = 1 if include_per_state_hist else 5 / 6
        fig.set_size_inches(8.27, 11.69 * vscale)  # portrait A4
        if borders:
            fig.set_facecolor('#ddd')
        gs = oaxs[0, 0].get_subplotspec().get_gridspec()
        for ax in oaxs.ravel():
            ax.remove()
        gridspecs = dict(
            per_state_histograms=gs[0, :],
            timeplots=gs[s + 2 :, 2],
            noodle=gs[s + 0 : s + 2, 1:],
            separated_spectra_and_hists=gs[s + 0 :, 0],
            nacs_histograms=gs[s + 3 :, 1],
            structure=gs[s + 2, 1],
        )
        if not include_per_state_hist:
            del gridspecs['per_state_histograms']
        sfs = {name: fig.add_subfigure(sgs) for name, sgs in gridspecs.items()}
        return fig, sfs

    def plot(
        self,
        include_per_state_hist: bool = False,
        borders: bool = False,
        consistent_lettering: bool = True,
    ) -> Figure:
        """Function to plot this Datasheet.

        Will generate all subplots and calculate necessary data if it has not yet been generated.

        Args:
            include_per_state_hist (bool, optional): Flag whether per-state histograms should be included. Defaults to False.
            borders (bool, optional): Flag whether the figure should have borders or not. Defaults to False.
            consistent_lettering (bool, optional): Flag whether consistent lettering should be used, i.e. whether the same plot should always have the same label letter. Defaults to True.

        Returns:
            Figure: The figure holding the entirety of plots in this Datasheet page.
        """
        letters = iter('abcdef')

        def outlabel(ax):
            nonlocal letters
            fixedtrans = mpl.transforms.ScaledTranslation(
                -20 / 72, +7 / 72, ax.figure.dpi_scale_trans
            )
            transform = ax.transAxes + fixedtrans
            return ax.text(
                0.0,
                1.0,
                next(letters) + ")",
                transform=transform,
                va='bottom',
                fontweight='bold',
                bbox=dict(facecolor='0.9', edgecolor='none', pad=3.0),
            )

        def inlabel(ax):
            nonlocal letters
            return ax.annotate(
                next(letters) + ")",
                xy=(0, 1),
                xycoords='axes fraction',
                xytext=(+0.5, -0.5),
                textcoords='offset fontsize',
                va='top',
                fontweight='bold',
                bbox=dict(facecolor='0.9', edgecolor='none', pad=3.0),
            )

        fig, sfs = self.get_subfigures(
            include_per_state_hist=include_per_state_hist, borders=borders
        )

        # print(self.frames)

        # separated_spectra_and_hists
        if self.can['separated_spectra_and_hists']:
            axs = self.plot_separated_spectra_and_hists(
                fig=sfs['separated_spectra_and_hists']
            )
            ax = axs['sg']
            outlabel(ax)
        else:
            ax = sfs['separated_spectra_and_hists'].subplots(1, 1)
            centertext(r"No $\mathbf{\mu}_{ij}$ data", ax=ax)
            ax.get_yaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)
            inlabel(ax)
        # noodle
        if self.can['noodle']:
            ax = self.plot_noodle(fig=sfs['noodle'])
            inlabel(ax)
        elif consistent_lettering:
            next(letters)
        # structure
        if self.can['structure']:
            ax = self.plot_structure(fig=sfs['structure'])
            inlabel(ax)
        elif consistent_lettering:
            next(letters)
        # nacs_histograms
        if self.can['nacs_histograms']:
            axs = self.plot_nacs_histograms(fig=sfs['nacs_histograms'])
            ax = axs.get('ntd', axs['nde'])
            outlabel(ax)
        elif consistent_lettering:
            next(letters)
        # time plots
        if self.can['timeplots']:
            axs = self.plot_timeplots(fig=sfs['timeplots'])
            ax = axs['pop']
            outlabel(ax)
        elif consistent_lettering:
            next(letters)
        if include_per_state_hist:
            axs = self.plot_per_state_histograms(fig=sfs['per_state_histograms'])
            ax = axs['energy']
            outlabel(ax)
        elif consistent_lettering:
            next(letters)
        return fig

    def _test_subfigures(
        self, include_per_state_hist: bool = False, borders: bool = False
    ):
        """Helper method to test whether subfigures are successfully plotted

        Args:
            include_per_state_hist (bool, optional): Flag to include per-state histograms. Defaults to False.
            borders (bool, optional): Whether the figures should have borders. Defaults to False.
        """
        fig, sfs = self.get_subfigures(
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
