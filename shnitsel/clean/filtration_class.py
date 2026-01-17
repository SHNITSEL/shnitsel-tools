from functools import cached_property
import xarray as xr

from shnitsel.data.dataset_containers import wrap_dataset, ShnitselDataset


class Filtration:
    def __init__(self, subject, filtranda=None):
        # Setting subjects to filter
        self.subject_original = subject
        if isinstance(subject, xr.Dataset):
            self.subject_wrapped = wrap_dataset(subject)
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

        self.trajectory_dim = (
            'atrajectory'
            if 'atrajectory' in self.subject_dataset.dims
            # If unstacked:
            else 'trajectory'
            if {'time', 'trajectory'}.issubset(self.subject_dataset.dims)
            else 'trajid'
            if 'trajid' in self.subject_dataset.dims
            else ''
        )

    @cached_property
    def noncumulative_mask(self):
        return self.filtranda < self.thresholds

    @cached_property
    def cumulative_mask(self):
        if self.trajectory_dim:
            res = self.noncumulative_mask.groupby(self.trajectory_dim).cumprod()
        else:
            res = self.noncumulative_mask.cumprod(self.subject_wrapped.leading_dim)
        return res.astype(bool)
    
    # @cached_property
    # def good_throughout(self): ...

    def truncate(self):
        """Perform a truncation, i.e. cut off the trajectory at its last continuously valid frame from the begining."""
        filter_mask_all_criteria = self.cumulative_mask.all('criterion')
        return type(self.subject_original)(
            self.subject_dataset.isel(
                {self.subject_wrapped.leading_dim: filter_mask_all_criteria}
            )
        )

    # def omit(self):


    def transect(self, cutoff_time): ...

    # NOTE (thevro): The following is not a property because it can fail if there is no time dimension
    # def good_upto(self):
    #     from shnitsel.clean.common import true_upto

    #     return true_upto(self.cumulative_mask)


    def plot_thresholds(self): ...

    def plot_populations(self): ...
