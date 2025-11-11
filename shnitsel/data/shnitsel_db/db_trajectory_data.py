from collections.abc import Callable
from typing import List, Self, TypeVar
import xarray as xr

from shnitsel.data.trajectory_format import Trajectory
from .datatree_level import _datatree_level_attribute_key

T = TypeVar("T")


class TrajectoryData(xr.DataTree):
    """DataTree node to keep track of a single trajectory entry"""

    def __init__(
        self,
        dataset: xr.Dataset | Trajectory | None = None,
        name: str | None = None,
    ):
        super().__init__(dataset=dataset, children=None, name=name)

        self.attrs[_datatree_level_attribute_key] = "TrajectoryData"

    def collect_trajectories(self) -> List[Self]:
        """Function to retrieve all trajectories in this subtree

        Returns:
            List[TrajectoryData]: List of all nodes with TrajectoryData type
        """
        return [self.copy()]

    def map_over_trajectories(
        self, map_func: Callable[[Trajectory], T], result_as_dict=False, result_var_name:str='result'
    ) -> Self | dict:
        """Method to apply a function to all trajectories in this subtree.

        Args:
            map_func (Callable[[Trajectory], T]): Function to be applied to each individual trajectory in this database structure.
            result_as_dict (bool, optional): Whether to return the result as a dict or as a TrajectoryData structure. Defaults to False which yields a CompoundGroup.
            result_var_name (str,optional): The name of the result variable to be assigned in either the result dataset or in the result dict.

        Returns:
            TrajectoryData|dict: The result, either again as a TrajectoryData structure or as a layered dict structure.
        """
        res = map_func(self.dataset)
        if result_as_dict:
            # res_dict["_group_info"] = self.get_group_info()
            return {result_var_name: res}
        else:
            if not isinstance(res, xr.Dataset):
                if isinstance(res, xr.DataArray):
                    res = res.to_dataset(name=result_var_name)
                else:
                    da = xr.DataArray(res)
                    res = da.to_dataset(name=result_var_name)

            return type(self)(
                res,
                None,
            ) # type: ignore
