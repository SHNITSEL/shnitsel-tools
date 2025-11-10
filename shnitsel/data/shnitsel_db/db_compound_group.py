from dataclasses import dataclass
import logging
from typing import Callable, List, Mapping, Self

from .db_trajectory_group import TrajectoryGroup

from .db_trajectory_data import TrajectoryData
import xarray as xr

from .datatree_level import _datatree_level_attribute_key


@dataclass
class CompoundInfo:
    """Class to hold identifying and auxiliary info of a compound type in ShnitselDB"""

    compound_name: str = "unknown"
    smile_repr: str | None = None


class CompoundGroup(xr.DataTree):
    """DataTree node to keep track of all data associated with a common compound within the datatree"""

    def __init__(
        self,
        compound_info: CompoundInfo | None = None,
        children: Mapping[
            str,
            TrajectoryGroup | TrajectoryData,
        ]
        | None = None,
    ):
        super().__init__(
            None,
            children,
            compound_info.compound_name if compound_info is not None else None,
        )
        self.attrs[_datatree_level_attribute_key] = "CompoundGroup"
        if compound_info is not None:
            self.attrs["compound_info"] = compound_info.__dict__

    def get_compound_info(self) -> CompoundInfo:
        """Get the store compound info of this Compound group.

        Returns:
            CompoundInfo: _description_
        """
        if "compound_info" in self.attrs:
            compound_data = self.attrs["compound_info"]
            return CompoundInfo(**compound_data)
        else:
            return CompoundInfo()

    def collect_trajectories(self):
        """Function to retrieve all trajectories in this subtree

        Returns:
            List[TrajectoryData]: List of all nodes with TrajectoryData type
        """
        res = []

        for x in self.children.values():
            res += x.collect_trajectories()

        return res

    def merge_with(self, other: Self) -> Self:
        """Function to merge to compound groups into one.

        Called when merging two database states. Will fail if compound_info differs between compounds to avoid loss of information.

        Args:
            other (CompoundGroup): The other CompoundGroup to be merged

        Raises:
            ValueError: Raised if the compound_info differs.

        Returns:
            CompoundGroup: A CompoundGroup object holding the entire merged subtree
        """
        own_info = self.attrs["compound_info"]
        other_info = other.attrs["compound_info"]
        if other_info != own_info:
            message = f"Cannot merge compounds with conflicting compound information: {other_info} vs. {own_info}"
            logging.error(message)
            raise ValueError(message)

        return CompoundGroup(
            own_info, {f"{i}": v for i, v in enumerate(self.collect_trajectories())}
        )  # type: ignore

    def filter_trajectories(
        self,
        filter_func: Callable[[xr.Dataset], bool] | None = None,
        est_level: str | List[str] | None = None,
        basis_set: str | List[str] | None = None,
        **kwargs,
    ) -> Self | None:
        """Function to filter trajectories based on their attributes.

        Args:
            filter_func (Callable[[xr.Dataset], bool] | None, optional): A function to evaluate whether a trajectory should be retained. Should return True if the trajectory should stay in the filtered set. Defaults to None.
            est_level (str | List[str] | None, optional): Option to filter for a certain level of electronic structure theory/calculation method. Can be a single key value or a set of values to retain. Defaults to None.
            basis_set (str | List[str] | None, optional): Option to filter for a certain basis set. Can be a single key value or a set of values to retain. Defaults to None.
            **kwargs: Key-value pairs, where the key denotes an attribute

        Returns:
            Compoundgroup|None: Either returns the CompoundGroup with the remaining set of trajectories or None if the group would be empty.
        """

        if filter_func is None:
            filter_func = lambda x: True

        if isinstance(est_level, str):
            est_level = list(est_level)

        if isinstance(basis_set, str):
            basis_set = list(basis_set)

        filter_vals = {
            k: list(v) if isinstance(v, str) else v for k, v in kwargs.items()
        }
        filter_vals["est_level"] = est_level
        filter_vals["basis_set"] = basis_set

        def composed_filter(data: xr.Dataset) -> bool:
            if not filter_func(data):
                return False

            for k, v in filter_vals.items():
                if k not in data.attrs or data.attrs[k] not in v:
                    return False

            return True

        filtered_traj = [
            t
            for t in self.collect_trajectories()
            if t.dataset is not None and composed_filter(t.dataset)
        ]

        if len(filtered_traj) > 0:
            return type(self)(
                self.attrs["compound_info"],
                {f"{i}": v for i, v in enumerate(filtered_traj)},
            )
        else:
            return None
