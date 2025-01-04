import xarray as xr

from typing import Hashable, Tuple, List


class TrajIndex(xr.Index):
    def __init__(self, traj_sizes: List[Tuple[Hashable, int]]):
        # self.size_dict = {trajid: size for trajid, size in traj_sizes}
        # self.offsets = np.cumsum([size for trajid, size in traj_sizes])
        ...

    @classmethod
    def from_variables(cls, variables):
        assert len(variables) == 1

    def sel(self, labels): ...