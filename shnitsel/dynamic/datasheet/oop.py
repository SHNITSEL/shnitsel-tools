from functools import cached_property
import xarray as xr
import numpy as np

from .. import postprocess as P
from .. import xrhelpers as xh

from .calc import calc_spectra

# from .plot import create_subfigures
from .plot.time import plot_timeplots


class Datasheet:
    def __init__(self, *, path: str | None = None, frames: xr.Dataset | None = None):
        if path is not None and frames is not None:
            raise ValueError("`path` and `frames` should not both be set")
        elif frames is not None:
            self.frames = frames
        elif path is not None:
            self.frames = xh.get_frames(path)
        else:
            print("Neither path nor frames given, please set frames manually")

        self.col_state = ['#4DAD15', '#AD2915', '#7515AD']
        self.col_inter = ['#2c3e50', '#C4A000', '#7E5273']

    times: np.ndarray

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
        return res

    @cached_property
    def fosc_time(self):
        res = P.time_grouped_ci(self.inter_state['fosc'])
        res['_color'] = 'statecomb', self.col_inter
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

    # @cached_property
    # def axs(self):

    def plot_timeplots(self):
        return plot_timeplots(
            pops=self.pops,
            delta_E=self.delta_E,
            fosc_time=self.fosc_time,
        )
