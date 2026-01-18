from functools import cached_property
from typing import Sequence

import numpy as np
import xarray as xr

from shnitsel.data.dataset_containers import wrap_dataset, ShnitselDataset
from shnitsel.data.multi_indices import ensure_unstacked, sel_trajs

class Filtration:
    def __init__(self, subject, filtranda=None):
        # Setting subjects to filter
        self.subject_original = subject
        if isinstance(subject, xr.Dataset):
            self.subject_wrapped = wrap_dataset(subject, expected_types=None)
            self.subject_dataset = subject
        elif isinstance(subject, ShnitselDataset):
            self.subject_wrapped = subject
            self.subject_dataset = subject.dataset
        else:
            raise ValueError(
                "Please pass an xarray.Dataset or a subclass of ShnitselDataset"
            )

        self.filtranda = subject['filtranda'] if filtranda is None else filtranda
        if 'thresholds' not in self.filtranda.coords:
            raise ValueError("The filtranda object should contain a 'thresholds' coord")
        self.thresholds = filtranda['thresholds']

        self.leading_dim = self.subject_wrapped.leading_dim
        self.trajectory_groupable = (
            'atrajectory'
            if 'atrajectory' in self.subject_dataset.coords
            # If unstacked:
            else 'trajectory'
            if {'time', 'trajectory'}.issubset(self.subject_dataset.coords)
            else 'trajid'
            if 'trajid' in self.subject_dataset.coords
            else ''
        )

        # FIXME: This may not generalize:
        self.unstacked_trajectory_dim = self.trajectory_groupable

    @cached_property
    def noncumulative_mask(self):
        return self.filtranda < self.thresholds

    @cached_property
    def cumulative_mask(self):
        if self.trajectory_groupable:
            res = self.noncumulative_mask.groupby(self.trajectory_groupable).cumprod()
        else:
            res = self.noncumulative_mask.cumprod(self.leading_dim)
        return res.astype(bool)

    @cached_property
    def good_throughout(self):
        if self.trajectory_groupable:
            res = self.noncumulative_mask.groupby(self.trajectory_groupable).all(
                self.leading_dim
            )
        else:
            res = self.noncumulative_mask.all(self.leading_dim)
        return res.astype(bool)

    # NOTE (thevro): I'm hesitant to make this a property, because unstacking can take time
    def get_unstacked_dataset(self):
        unstacked, _ = ensure_unstacked(self.subject_original)
        return unstacked

    def get_unstacked_filtranda(self):
        unstacked_filtranda, _ = ensure_unstacked(self.filtranda)
        return unstacked_filtranda

    def truncate(self):
        """Perform a truncation, i.e. cut off the trajectory at its last continuously valid frame from the begining."""
        filter_mask_all_criteria = self.cumulative_mask.all('criterion')
        return type(self.subject_original)(
            self.subject_dataset.isel(
                {self.subject_wrapped.leading_dim: filter_mask_all_criteria}
            )
        )

    def omit(self):
        all_critera_fulfilled = self.good_throughout.all('criterion')
        if self.trajectory_groupable:
            return type(self.subject_original)(
                sel_trajs(self.subject_dataset, all_critera_fulfilled)
            )

        # Otherwise we have a single trajectory.
        return self.subject_original if all_critera_fulfilled.item() else None

    def transect(self, cutoff_time):
        indexer = {"time": slice(float(cutoff_time))}
        sliced = self.get_unstacked_dataset().loc[indexer]
        sliced_filtranda = self.get_unstacked_filtranda().loc[indexer]
        return Filtration(sliced, sliced_filtranda).omit()

    def get_unstacked_cumulative_maxima(self):
        return xr.apply_ufunc(
            np.maximum.accumulate,
            self.get_unstacked_filtranda(),
            input_core_dims=[['time']],
            output_core_dims=[['time']],
            kwargs={'axis': -1},
        )

    def cum_max_quantiles(
        self, quantiles: Sequence[float] | None = None
    ) -> xr.DataArray:
        if quantiles is None:
            quantiles = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]
        return self.get_unstacked_cumulative_maxima().quantile(
            quantiles, self.unstacked_trajectory_dim
        )

    def get_unstacked_cumulative_mask(self):
        unstacked_mask, _ = ensure_unstacked(self.cumulative_mask)
        return unstacked_mask

    def get_independent_populations(self):
        mask = (
            self.get_unstacked_cumulative_mask()
            .to_dataset('criterion')
            .assign({'total_population': mask.coords['is_frame']})
            .to_dataarray('criterion')
        )
        return mask  # FIXME

    def plot_thresholds(self): ...

    def plot_populations(self): ...

    # NOTE (thevro): The following is not a property because it can fail if there is no time dimension
    # def good_upto(self):
    #     from shnitsel.clean.common import true_upto

    #     return true_upto(self.cumulative_mask)