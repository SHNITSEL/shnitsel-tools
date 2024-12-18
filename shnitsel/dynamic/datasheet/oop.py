from functools import cached_property
import xarray as xr
import numpy as np
import rdkit.Chem as rc

from .. import postprocess as P
from .. import xrhelpers as xh

from .calc import calc_spectra

# from .plot import create_subfigures
from .plot.per_state_hist import plot_per_state_histograms
from .plot.time import plot_timeplots
# import
from .plot.nacs_hist import plot_nacs_histograms
from ..pca_biplot import plot_noodleplot, xyz_to_mol
from .plot.structure import plot_structure


class Datasheet:
    def __init__(
        self,
        *,
        path: str | None = None,
        frames: xr.Dataset | None = None,
        col_state: list | None = None,
        col_inter: list | None = None,
    ):
        if path is not None and frames is not None:
            raise ValueError("`path` and `frames` should not both be set")
        elif frames is not None:
            self.frames = frames
        elif path is not None:
            self.frames = xh.get_frames(path)
        else:
            print("Neither path nor frames given, please set frames manually")

        self.col_state = col_state or ['#4DAD15', '#AD2915', '#7515AD']
        self.col_inter = col_inter or ['#2c3e50', '#C4A000', '#7E5273']
        try:
            self.name = self.frames.attrs['longname']
        except KeyError:
            pass

        return None

    times: np.ndarray
    charge: int = 0
    structure_skeletal: bool = False
    name: str = ''

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
        return calc_spectra(self.inter_state, times=...)

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

    def plot_per_state_histograms(self):
        return plot_per_state_histograms(
            per_state=self.per_state,
        )

    def plot_timeplots(self):
        return plot_timeplots(
            pops=self.pops,
            delta_E=self.delta_E,
            fosc_time=self.fosc_time,
        )

    def plot_separated_spectra_and_hists(self):
        return NotImplemented

    def plot_nacs_histograms(self):
        return plot_nacs_histograms(self.inter_state, self.hops.frame)

    def plot_noodle(self):
        return plot_noodleplot(self.noodle, self.hops)

    def plot_structure(self):
        mol = self.mol_skeletal if self.structure_skeletal else self.mol
        return plot_structure(
            mol,
            name=self.name,
            smiles=self.smiles,
            inchi=self.inchi,
            ax=None,
        )