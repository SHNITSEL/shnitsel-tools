from functools import cached_property
import xarray as xr
import numpy as np
import rdkit.Chem as rc
import matplotlib.pyplot as plt

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
from ..pca_biplot import plot_noodleplot, xyz_to_mol
from .plot.structure import plot_structure


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
            self.frames = xh.get_frames(path)
        else:
            print("DEBUG -- new!")
            print("Neither path nor frames given, please set frames manually")

        if spectra_times is not None:
            self.spectra_times = spectra_times
        elif self.frames is not None:
            max_time = self.frames.coords['time'].max().item()
            self.spectra_times = [max_time * i / 3 for i in range(3)]

        self.col_state = col_state or ['#4DAD15', '#AD2915', '#7515AD']
        self.col_inter = col_inter or ['#2c3e50', '#C4A000', '#7E5273']
        try:
            self.name = self.frames.attrs['longname']
        except KeyError:
            pass

        return None

    spectra_times: np.ndarray
    charge: int = 0
    structure_skeletal: bool = False
    name: str = ''

    def _copy_data(self, old: Self):
        print("10:48")
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
        per_state = P.get_per_state(self.frames)
        per_state['_color'] = 'state', self.col_state
        return per_state

    @cached_property
    def inter_state(self):
        inter_state = P.get_inter_state(self.frames)
        inter_state['_color'] = 'statecomb', self.col_inter
        if 'dip_trans' in inter_state:
            inter_state = P.assign_fosc(inter_state)
        for var, tex in [
            ('energies', r"$\Delta E$"),
            ('nacs', r"$\|\mathrm{NAC}_{i,j}\|_2$"),
            ('dip_trans', r"$\|\mathbf{\mu}_{i,j}\|_2$"),
            ('fosc', r"$f_\mathrm{osc}$"),
        ]:
            try:
                inter_state[var].attrs['tex'] = tex
            except KeyError:
                pass
        return inter_state

    @cached_property
    def pops(self):
        pops = P.calc_pops(self.frames)
        pops['_color'] = 'state', self.col_state
        return pops

    @cached_property
    def delta_E(self):
        res = P.time_grouped_ci(self.inter_state['energies'])
        res['_color'] = 'statecomb', self.col_inter
        res.attrs['tex'] = r"$\Delta E$"
        return res

    @cached_property
    def fosc_time(self):
        res = P.time_grouped_ci(self.inter_state['fosc'])
        res['_color'] = 'statecomb', self.col_inter
        res.attrs['tex'] = r"$f_\mathrm{osc}$"
        return res

    @cached_property
    def spectra(self):
        return calc_spectra(self.inter_state, times=self.spectra_times)

    @cached_property
    def spectra_groups(self):
        return get_sgroups(self.spectra)

    @cached_property
    def spectra_ground(self):
        return self.spectra_groups[0]

    @cached_property
    def spectra_excited(self):
        return self.spectra_groups[1]

    @cached_property
    def noodle(self):
        return P.pairwise_dists_pca(self.frames.atXYZ)

    @cached_property
    def hops(self):
        mask = P.hop_indices(self.frames.astate)
        return self.noodle[mask]

    @cached_property
    def structure_atXYZ(self):
        return self.frames.atXYZ.isel(frame=0)

    @cached_property
    def mol(self):
        return xyz_to_mol(self.structure_atXYZ, charge=self.charge)

    @cached_property
    def mol_skeletal(self):
        return rc.RemoveHs(self.mol)

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
        self.pops
        self.delta_E
        self.fosc_time
        self.spectra  # -> inter_state
        self.spectra_groups
        self.noodle
        self.hops
        self.structure_atXYZ
        self.mol_skeletal
        self.smiles
        self.inchi

    def plot_per_state_histograms(self, fig: Figure | None = None):
        return plot_per_state_histograms(
            per_state=self.per_state,
            fig=fig,
        )

    def plot_timeplots(self, fig: Figure | None = None):
        return plot_timeplots(
            pops=self.pops,
            delta_E=self.delta_E,
            fosc_time=self.fosc_time,
        )

    def plot_separated_spectra_and_hists(self, fig: Figure | None = None):
        return plot_separated_spectra_and_hists(
            inter_state=self.inter_state,
            sgroups=self.spectra_groups,
        )

    def plot_nacs_histograms(self, fig: Figure | None = None):
        return plot_nacs_histograms(self.inter_state, self.hops.frame)

    def plot_noodle(self, fig: Figure | None = None):
        return plot_noodleplot(self.noodle, self.hops)

    def plot_structure(self, fig: Figure | None = None):
        mol = self.mol_skeletal if self.structure_skeletal else self.mol
        return plot_structure(
            mol,
            name=self.name,
            smiles=self.smiles,
            inchi=self.inchi,
            ax=None,
        )

    # TODO: Make whole datasheet and modularize creation process
    # 1. Get it working with existing `create_subfigures` function
    # 2. Adjust all the other functions to take a `figure` argument
    # 3. Modify `create_subfigures` so that it doesn't make all the axes at once

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
        fig, sfs = self.get_subfigures(
            include_per_state_hist=include_per_state_hist, borders=borders
        )
        if include_per_state_hist:
            self.plot_per_state_histograms(fig=sfs['per_state_histograms'])
        self.plot_timeplots(fig=sfs['timeplots'])
        self.plot_separated_spectra_and_hists(fig=sfs['separated_spectra_and_hists'])
        self.plot_nacs_histograms(fig=sfs['nacs_histograms'])
        self.plot_noodle(fig=sfs['noodle'])
        self.plot_structure(fig=sfs['structure'])
        return fig
