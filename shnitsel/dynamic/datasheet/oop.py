from functools import cached_property
import xarray as xr
import numpy as np
import rdkit.Chem as rc
import matplotlib as mpl
import matplotlib.pyplot as plt
from logging import info, warning
from timeit import default_timer as timer

from matplotlib.figure import Figure
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from .. import postprocess as P
from .. import xrhelpers as xh

from .calc import calc_spectra, get_sgroups

# from .plot import create_subfigures
from .plot.per_state_hist import plot_per_state_histograms
from .plot.time import plot_timeplots
# import
from .plot.dip_trans_hist import plot_separated_spectra_and_hists
from .plot.nacs_hist import plot_nacs_histograms
from ..pca_biplot import plot_noodleplot  # , xyz_to_mol
from .plot.structure import plot_structure, xyz_to_mol


class Datasheet:
    def __init__(
        self,
        *,
        path: str | None = None,
        frames: xr.Dataset | None = None,
        copy_data: Self | None = None,
        spectra_times: list[int | float] | None = None,
        col_state: list | None = None,
        col_inter: list | None = None,
    ):
        if copy_data is not None:
            if path is not None or frames is not None:
                raise ValueError(
                    "if `copy` is set, neither `path` nor `frames` should be set"
                )
            self._copy_data(old=copy_data)
            return

        if path is not None and frames is not None:
            raise ValueError("`path` and `frames` should not both be set")
        elif frames is not None:
            self.frames = frames
        elif path is not None:
            self.frames = xh.open_frames(path)
        else:
            warning("Neither path nor frames given, please set frames manually")

        if spectra_times is not None:
            self.spectra_times = spectra_times
        elif self.frames is not None:
            max_time = self.frames.coords['time'].max().item()
            self.spectra_times = [max_time * i / 3 for i in range(3)]

        if col_state is not None:
            assert (ncols := len(col_state)) == (nstates := frames.sizes['state']), (
                f"`col_state` has {ncols} colors,"
                f"but should contain one color for each of the {nstates} states"
            )
            self.col_state = col_state
        elif (s := self.frames.sizes['state']) <= 3:
            self.col_state = ['#4DAD15', '#AD2915', '#7515AD'][:s]  # SHNITSEL-colours
        elif s <= 10:
            cmap = plt.get_cmap('tab10')
            self.col_state = [mpl.colors.rgb2hex(c) for c in cmap.colors][:s]
        else:
            raise ValueError(
                f"These data have {nstates} states. "
                "When passing data with more than 10 states, please"
                "also pass an appropriate colormap to `col_state`."
            )

        if col_inter is not None:
            assert (ncols := len(col_inter)) == (ncombs := frames.sizes['statecomb']), (
                f"`col_inter` has {ncols} colors, "
                f"but should contain one color for each of the {ncombs} state combinations"
            )
        elif (s := self.frames.sizes['statecomb']) <= 3:
            self.col_inter = col_inter or ['#2c3e50', '#C4A000', '#7E5273'][:s]
        elif s <= 10:
            # TODO: choose colours distinct from per_state colours
            cmap = plt.get_cmap('tab10')
            self.col_inter = [mpl.colors.rgb2hex(c) for c in cmap.colors][:s]
        else:
            raise ValueError(
                f"These data have {ncombs} state combinations. "
                "When passing data with more than 10 state combinations, please "
                "also pass an appropriate colormap to `col_inter`."
            )

        try:
            self.name = self.frames.attrs['long_name']
        except KeyError:
            pass

        return None

    spectra_times: np.ndarray
    charge: int = 0
    structure_skeletal: bool = False
    name: str = ''

    def _copy_data(self, old: Self):
        self.spectra_times = old.spectra_times
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
    def per_state(self):
        start = timer()
        per_state = P.get_per_state(self.frames)
        per_state['_color'] = 'state', self.col_state
        end = timer()
        info(f"cached per_state in {end - start} s")
        return per_state

    @cached_property
    def inter_state(self):
        start = timer()
        inter_state = P.get_inter_state(self.frames)
        inter_state['_color'] = 'statecomb', self.col_inter
        if 'dip_trans' in inter_state:
            inter_state = P.assign_fosc(inter_state)
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
    def pops(self):
        start = timer()
        pops = P.calc_pops(self.frames)
        pops['_color'] = 'state', self.col_state
        end = timer()
        info(f"cached pops in {end-start} s")
        return pops

    @cached_property
    def delta_E(self):
        start = timer()
        res = P.time_grouped_ci(self.inter_state['energy'])
        res['_color'] = 'statecomb', self.col_inter
        res.attrs['tex'] = r"$\Delta E$"
        end = timer()
        info(f"cached delta_E in {end-start} s")
        return res

    @cached_property
    def fosc_time(self):
        start = timer()
        res = P.time_grouped_ci(self.inter_state['fosc'])
        res['_color'] = 'statecomb', self.col_inter
        res.attrs['tex'] = r"$f_\mathrm{osc}$"
        end = timer()
        info(f"cached fosc_time in {end-start} s")
        return res

    @cached_property
    def spectra(self):
        start = timer()
        res = calc_spectra(self.inter_state, times=self.spectra_times)
        end = timer()
        info(f"cached spectra in {end-start} s")
        return res

    @cached_property
    def spectra_groups(self):
        start = timer()
        res = get_sgroups(self.spectra)
        end = timer()
        info(f"cached spectra_groups in {end-start} s")
        return res

    @cached_property
    def spectra_ground(self):
        return self.spectra_groups[0]

    @cached_property
    def spectra_excited(self):
        return self.spectra_groups[1]

    @cached_property
    def noodle(self):
        start = timer()
        res = P.pairwise_dists_pca(self.frames.atXYZ)
        end = timer()
        info(f"cached noodle in {end-start} s")
        return res

    @cached_property
    def hops(self):
        mask = P.hop_indices(self.frames.astate)
        return self.noodle[mask]

    @cached_property
    def structure_atXYZ(self):
        return self.frames.atXYZ.isel(frame=0)

    @cached_property
    def mol(self):
        if 'smiles_map' in self.frames['atXYZ'].attrs:
            mol = P.numbered_smiles_to_mol(self.frames['atXYZ'].attrs['smiles_map'])
            for atom in mol.GetAtoms():
                atom.ClearProp("molAtomMapNumber")
                atom.SetProp("atomNote", str(atom.GetIdx()))
            return mol
        else:
            return xyz_to_mol(self.structure_atXYZ, charge=self.charge)

    @cached_property
    def mol_skeletal(self):
        mol = rc.Mol(self.mol)
        return rc.RemoveHs(mol)

    @cached_property
    def smiles(self):
        return rc.MolToSmiles(self.mol_skeletal)

    @cached_property
    def inchi(self):
        return rc.MolToInchi(self.mol_skeletal)

    # @cached_property
    # def axs(self):

    def calc_all(self):
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

    def plot_per_state_histograms(self, fig: Figure | None = None):
        start = timer()
        res = plot_per_state_histograms(
            per_state=self.per_state,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_per_state_histograms in {end - start} s")
        return res

    def plot_timeplots(self, fig: Figure | None = None):
        start = timer()
        res = plot_timeplots(
            pops=self.pops,
            delta_E=self.delta_E,
            fosc_time=self.fosc_time,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_timeplots in {end - start} s")
        return res

    def plot_separated_spectra_and_hists(self, fig: Figure | None = None):
        start = timer()
        res = plot_separated_spectra_and_hists(
            inter_state=self.inter_state,
            sgroups=self.spectra_groups,
            fig=fig,
        )
        end = timer()
        info(f"finished plot_separated_spectra_and_hists in {end - start} s")
        return res

    def plot_nacs_histograms(self, fig: Figure | None = None):
        start = timer()
        res = plot_nacs_histograms(self.inter_state, self.hops.frame, fig=fig)
        end = timer()
        info(f"finished plot_nacs_histograms in {end - start} s")
        return res

    def plot_noodle(self, fig: Figure | None = None):
        start = timer()
        res = plot_noodleplot(self.noodle, self.hops, fig=fig)
        end = timer()
        info(f"finished plot_noodle in {end - start} s")
        return res

    def plot_structure(self, fig: Figure | None = None):
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
    def get_subfigures(include_per_state_hist: bool = False, borders: bool = False):
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

    def plot(self, include_per_state_hist: bool = False, borders: bool = False):
        def label(ax, x, y, ha, va='baseline'):
            nonlocal letters
            return ax.text(
                x,
                y,
                f"{next(letters)})",
                fontweight='bold',
                transform=ax.transAxes,
                ha=ha,
                va=va,
            )

        def outlabel(ax):
            return label(ax, -0.05, 1.05, ha='right')

        def inlabel(ax):
            return label(ax, 0.05, 0.95, ha='left', va='top')

        fig, sfs = self.get_subfigures(
            include_per_state_hist=include_per_state_hist, borders=borders
        )
        letters = iter('abcdef')
        if include_per_state_hist:
            self.plot_per_state_histograms(fig=sfs['per_state_histograms'])
        ## separated_spectra_and_hists
        if 'dip_trans' in self.frames.data_vars:
            axs = self.plot_separated_spectra_and_hists(
                fig=sfs['separated_spectra_and_hists']
            )
            ax = axs['sg']
            outlabel(ax)
        ## noodle
        ax = self.plot_noodle(fig=sfs['noodle'])
        inlabel(ax)
        ## structure
        ax = self.plot_structure(fig=sfs['structure'])
        inlabel(ax)
        ## nacs_histograms
        axs = self.plot_nacs_histograms(fig=sfs['nacs_histograms'])
        ax = axs.get('ntd', axs['nde'])
        outlabel(ax)
        ## time plots
        axs = self.plot_timeplots(fig=sfs['timeplots'])
        ax = axs['pop']
        outlabel(ax)
        return fig

    def _test_subfigures(
        self, include_per_state_hist: bool = False, borders: bool = False
    ):
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
