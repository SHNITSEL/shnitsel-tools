from dataclasses import dataclass
import logging
from typing import (
    Callable,
    List,
    Literal,
    Mapping,
    Self,
    TypeAlias,
    TypeVar,
    overload,
)
import numpy as np
import xarray as xr

from shnitsel.data.shnitsel_db.db_trajectory_group import GroupInfo

from .shnitsel_db import (
    CompoundGroup,
    TrajectoryData,
    TrajectoryGroup,
    CompoundInfo,
    _datatree_level_attribute_key,
)
from shnitsel.data.trajectory_format import Trajectory

T = TypeVar("T")


@dataclass
class MetaInformation:
    input_format: Literal["sharc", "newtonx", "ase", "pyrai2md"] | None = None
    input_type: Literal['static', 'dynamic'] | None = None
    input_format_version: str | None = None
    theory_basis_set: str | None = None
    est_level: str | None = None


class ShnitselDBRoot(xr.DataTree):
    """DataTree root node to keep track of Shnitsel data in a predefined structural hierarchy"""

    def __init__(self, compounds: Mapping[str, CompoundGroup] | None = None):
        super().__init__(dataset=None, children=compounds, name="ROOT")
        self.attrs[_datatree_level_attribute_key] = "ShnitselDBRoot"

    def is_level(self, target_level: str) -> bool:
        """Check whether we are at a certain level

        Args:
            target_level (str): Desired level to check for

        Returns:
            bool: True if this level satisfies the requirements
        """
        return target_level == "ShnitselDBRoot"

    def add_trajectory_group(
        self,
        group_name: str,
        filter_func_compound: Callable[[CompoundInfo], bool] | None = None,
        filter_func_trajectories: Callable[[Trajectory | TrajectoryGroup], bool]
        | None = None,
        flatten_compound_trajectories=False,
        **kwargs,
    ) -> Self:
        """Function to group ungrouped trajectories into a new group based on filter conditions.

        `filter_func_compound` can be used to only generate the group for certain compounds.
        This parameter should be a function that only returns True if the group should be created underneath this comound.
        `filter_func_trajectories` can be used to select only specific trajectories and groups underneath a compound to be part of this group.
        `flatten_compound_trajectories` can be set to `True` if existing groups within a compound are supposed to be dissolved (i.e. all trajectories gathered and put directly as children of the Compound)

        Args:
            group_name (str): The name of the group to be created.
            filter_func_compound (Callable[[CompoundInfo], bool] | None, optional): Filter function that should return True if the group should be created for this compound. Defaults to None.
            filter_func_trajectories (Callable[[Trajectory|TrajectoryGroup], bool] | None, optional): Filter function to determine whether a group or trajectory should be included in the new group. Defaults to None.
            flatten_compound_trajectories (bool, optional): Flag to determine whether all trajectories under selected compounds should be ungrouped before selecting for the new group. Defaults to False.
            **kwargs: Key-value pairs for the group that should be set in the group's key-value dict.

        Returns:
            ShnitselDB: A resulting ShnitselDB structure with the grouping applied.
        """
        new_children = {}

        for child_key, child in self.children.items():
            if filter_func_compound is None or filter_func_compound(
                child.get_compound_info()
            ):
                new_children[child_key] = child.add_trajectory_group(
                    group_name,
                    filter_func_trajectories,
                    flatten_compound_trajectories,
                    **kwargs,
                )
        return type(self)(new_children)

    def set_compound_info(
        self, compound_info: CompoundInfo, apply_to_all=False
    ) -> Self:
        """Function to set the compound information on either all unknown compounds or on all trajectories (if `apply_to_all=True`).

                By default, the compound info will only be applied to trajectories with unknown compounds.

                Args:
                    compound_info (CompoundInfo): The compound information.
                    apply_to_all (bool, optional): A flag to
        ShnitselDB: TypeAlias = ShnitselDBRoot specify that the compound information should be set to all trajectories in this tree. Use with caution! Defaults to False.

                Raises:
                    ValueError: _description_

                Returns:
                    ShnitselDBRoot: The updated database
        """
        # TODO: FIXME: Deal with clashes when flattening
        if apply_to_all:
            if len(self.children.keys()) > 1:
                logging.warning(
                    "Forcing the compound to be applied to all compounds and not just the `unknown` group is considered unsafe."
                )

            new_trajectories = {}
            traj_counter = 0

            for k, c in self.children.items():
                # iterate over compounds
                for _, t in c.children.items():
                    # Iterate over trajectories
                    key = f"{traj_counter}"
                    traj_counter += 1
                    new_trajectories[key] = t.copy()

            return type(self)(
                {
                    compound_info.compound_name: CompoundGroup(
                        compound_info, new_trajectories
                    )
                }
            )
        else:
            if compound_info.compound_name not in self.children:
                res_group = CompoundGroup(compound_info, None)
            else:
                res_group = self.children[compound_info.compound_name].copy()

            if "unknown" in self.children:
                unknown_trajectories = {
                    k: v.copy() for k, v in self.children["unknown"].children.items()
                }
                # traj_counter = (
                #    max(0, np.max([int(x) for x in res_group.children.keys()])) + 1
                # )
                res_group.update(unknown_trajectories)

            new_children: dict[str, CompoundGroup] = {
                k: v.copy() for k, v in self.children.items() if str(k) != "unknown"
            }  # type: ignore
            new_children[compound_info.compound_name] = res_group  # type: ignore

            return type(self)(new_children)

    def apply_trajectory_setup_properties(self, properties: MetaInformation) -> None:
        """Method to apply attributes on all trajectories in this tree

        Args:
            properties (MetaInformation): The attributes to set with their respective values.
        """
        if len(self.children.keys()) > 1:
            logging.error(
                "Cannot set trajectory setup properties if more than one compound is registered in database. Please set properties on each compound individually"
            )
            return

        props = {}

        for k, v in properties.__dict__.items():
            if v is not None:
                props[k] = v

        for leaf in self.leaves:
            leaf.dataset.attrs.update(props)

    def merge_with(self, other: Self) -> Self:
        """Function to merge two ShnitselDBRoots into one

        Called when merging two database states. Will fail if compounds in the two states have different compound_info.
        Attention: Trajectories that occur in both states will result in duplicate trajectories in the final result.

        Args:
            other (CompoundGroup): The other CompoundGroup to be merged

        Raises:
            ValueError: Raised if the compound_info differs on some CompoundGroup.

        Returns:
            ShnitselDBRoot: A ShnitselDBRoot object representing the merged database state"""

        total_compound_keys = set(self.children.keys())
        total_compound_keys.update(other.children.keys())
        shared_keys = set(self.children.keys())
        shared_keys = shared_keys.intersection(other.children.keys())

        total_compounds = {}

        for k in total_compound_keys:
            if k in shared_keys:
                total_compounds[k] = self.children[k].merge_with(other.children[k])
            elif k in self.children:
                total_compounds[k] = self.children[k].copy()
            else:
                total_compounds[k] = other.children[k].copy()

        return type(self)(total_compounds)

    def filter_compounds(
        self, compounds: str | List[str] | Callable[[CompoundInfo], bool]
    ) -> Self:
        """Function to filter out only a subset of compounds.

        Either retains exactly the requested compound if provided a string, only the set of provided compounds if provided a list of keys and all compounds for which the selector function returns `True` otherwise.
        The return value will never be None, but it may be a ShnitselDB tree without any trajectories in it.

        Args:
            compounds (str | Iterable[str] | Callable[[str],bool], optional): The name of a single compound, a list of names of compounds to retain or a function that says for each compound whether it will be retained.

        Returns:
            ShnitselDBRoot: A new database tree where only the filtered compounds are retained
        """

        filter_func = lambda x: True
        if callable(compounds):
            filter_func = compounds
        else:
            if isinstance(compounds, str):
                compounds = list(compounds)

            if isinstance(compounds, list):
                filter_func = lambda x: x.compound_name in compounds

        new_compounds = {}

        for k, v in self.children.items():
            if filter_func(v.attrs["compound_info"]):
                new_compounds[k] = v.copy()

        return type(self)(new_compounds)

    def filter_trajectories(
        self,
        filter_func: Callable[[xr.Dataset], bool] | None = None,
        est_level: str | List[str] | None = None,
        basis_set: str | List[str] | None = None,
        **kwargs,
    ) -> Self:
        """Function to filter trajectories based on their attributes.

        Args:
            filter_func (Callable[[xr.Dataset], bool] | None, optional): A function to evaluate whether a trajectory should be retained. Should return True if the trajectory should stay in the filtered set. Defaults to None.
            est_level (str | List[str] | None, optional): Option to filter for a certain level of electronic structure theory/calculation method. Can be a single key value or a set of values to retain. Defaults to None.
            basis_set (str | List[str] | None, optional): Option to filter for a certain basis set. Can be a single key value or a set of values to retain. Defaults to None.
            **kwargs: Key-value pairs, where the key denotes an attribute

        Returns:
            Self: _description_
        """
        new_compounds = {}

        for k, v in self.children.items():
            tmp_res = v.filter_trajectories(filter_func, est_level, basis_set, **kwargs)
            if tmp_res is not None:
                new_compounds[k] = tmp_res

        return type(self)(new_compounds)

    @overload
    def map_over_trajectories(
        self,
        map_func: Callable[[Trajectory], T],
        result_as_dict: Literal[False] = False,
        result_var_name: str = 'result',
        parallel: bool = False,
    ) -> Self:
        ...
        
    @overload
    def map_over_trajectories(
        self,
        map_func: Callable[[Trajectory], T],
        result_as_dict: Literal[True],
        result_var_name: str = 'result',
        parallel: bool = False,
    ) -> dict:
        ...
        
    def map_over_trajectories(
        self,
        map_func: Callable[[Trajectory], T],
        result_as_dict: bool = False,
        result_var_name: str = 'result',
        parallel: bool = False,
    ) -> Self | dict:
        """Method to apply a function to all trajectories

        Args:
            map_func (Callable[[Trajectory], T]): Function to be applied to each individual trajectory in this database structure.
            result_var_name (str,optional): The name of the result variable to be assigned in either the result dataset or in the result dict.
            result_as_dict (bool, optional): Whether to return the result as a dict or as a ShnitselDB structure. Defaults to False which yields a ShnitselDB.
            parallel (bool, optional): Whether the application of the function map_func can be evaluated in parallel.

        Returns:
            ShnitselDB|dict: The result, either again as a ShnitselDB structure or as a layered dict structure.
        """
        # TODO: FIXME: Implement parallel processing of map
        if result_as_dict:
            return {
                k: v.map_over_trajectories(map_func, result_as_dict, result_var_name)
                for k, v in self.children.items()
            }
        else:
            return type(self)(
                {
                    k: v.map_over_trajectories(
                        map_func, result_as_dict, result_var_name
                    )
                    for k, v in self.children.items()
                }
            )  # type: ignore


def collect_trajectories(
    db: xr.DataTree, only_direct_children: bool = False
) -> List[Trajectory]:
    """Function to retrieve all Trajectories in a subtree of a ShnitselDB

    Args:
        db (xr.DataTree): The subtree of the database in question
        only_direct (bool, optional): Whether to only gather trajectories from direct children of this subtree.

    Returns:
        List[Trajectory]: List of Trajectory entries in the Database
    """
    res: List[Trajectory] = []
    if only_direct_children:
        for key, child in db.children.items():
            if child.has_data:
                res.append(child.dataset)
    else:
        for leaf in db.leaves:
            if leaf.has_data:
                res.append(leaf.dataset)
    return res


def build_shnitsel_db(
    data: Trajectory | xr.DataTree | List,
) -> ShnitselDBRoot:
    """Function to generate a full -- i.e. up to ShnitselDBRoot -- Shnitsel DB structure.

    Wraps trajectories in xr.DataTree structures, converts xr.DataTree structures into Shnitsel-subclasses and converts lists of elements to either
    TrajectoryDate, TrajectoryGroup or Compoundgroup levels in a full ShnitselDB Tree

    Args:
        data (Trajectory | xr.DataTree | List[Trajectory  |  xr.DataTree]): Input data to be wrapped in a ShnitselDB format

    Raises:
        ValueError: If a list of ShnitselDBRoot objects was provided
        ValueError: If a list of xr.DataTree objects on different Levels of the ShnitselDB hierarchy was provided, e.g. a mix of TrajectoryData and CompoundGroup nodes.
        ValueError: If conversion of an xr.DataTree object yields an unknown ShnitselDB type.
        ValueError: If the provided data is of no ShnitselDB format type.

    Returns:
        ShnitselDB: The resulting ShnitselDB dataset structure.
    """
    if isinstance(data, ShnitselDBRoot):
        # If we already are a root node, return the structure.
        return data
    elif isinstance(data, Trajectory):
        # If we are provided a single trajectory wrap in a TrajectoryData node and build rest of tree
        return build_shnitsel_db(TrajectoryData(data, str(data.attrs.get("trajid", 0))))
    elif isinstance(data, list):
        # Convert individual nodes
        res_set = {}

        logging.debug(f"Building Shnitsel DB from {len(data)} data entries")
        shared_level_id = None
        for i, d in enumerate(data):
            if isinstance(d, xr.DataTree):
                d_conv = convert_shnitsel_tree(d)
                key = None

                if isinstance(d_conv, ShnitselDBRoot):
                    raise ValueError(
                        "Cannot build a Shnitsel DB structure from multipe DB roots."
                    )
                elif isinstance(d_conv, CompoundGroup):
                    if shared_level_id is None:
                        shared_level_id = "CompoundGroup"
                    elif shared_level_id != "CompoundGroup":
                        raise ValueError(
                            f"Cannot build Shnitsel DB structure from types {shared_level_id} and `CompoundGroup` on the same level of hierarchy."
                        )

                    key = d_conv.name if d_conv.name is not None else "unknown"
                elif isinstance(d_conv, TrajectoryData) or isinstance(
                    d_conv, TrajectoryGroup
                ):
                    if shared_level_id is None:
                        shared_level_id = "TrajectoryData|TrajectoryGroup"
                    elif shared_level_id != "TrajectoryData|TrajectoryGroup":
                        raise ValueError(
                            f"Cannot build Shnitsel DB structure from types {shared_level_id} and `TrajectoryData|TrajectoryGroup` on the same level of hierarchy."
                        )
                    key = (
                        d_conv.name
                        if d_conv.name is not None
                        else (
                            d_conv.dataset.attrs["trajid"]
                            if d_conv.has_data
                            and "trajid" in d_conv.dataset.attrs
                            else f"{i}"
                        )
                    )
                else:
                    raise ValueError(f"Found unsupported type: {type(d_conv)}")
            elif isinstance(d, Trajectory):
                # Map trajectory instance:
                if shared_level_id is None:
                    shared_level_id = "TrajectoryData|TrajectoryGroup"
                elif shared_level_id != "TrajectoryData|TrajectoryGroup":
                    raise ValueError(
                        f"Cannot build Shnitsel DB structure from types {shared_level_id} and `TrajectoryData|TrajectoryGroup` on the same level of hierarchy."
                    )
                key = f"{i}"
                d_conv = TrajectoryData(d, key)
            else:
                raise ValueError(f"Found unsupported type: {type(d)}")

            res_set[key] = d_conv

        logging.debug(f"Finalize Shnitsel DB from {len(res_set)} data entries")

        if shared_level_id == "TrajectoryData|TrajectoryGroup":
            tmp_compound = CompoundGroup(CompoundInfo(), children=res_set)
            return build_shnitsel_db(tmp_compound)
        elif shared_level_id == CompoundGroup:
            return ShnitselDBRoot(res_set)
        else:
            raise ValueError(
                "Could not find an appropriate level of the ShnitselDB hierarchy to build a full database."
            )

    elif isinstance(data, xr.DataTree):
        # We have a datatree instance, check for the DataTree_Level attribute to convert to right type.

        tmp_res = convert_shnitsel_tree(data)
        if isinstance(tmp_res, ShnitselDBRoot):
            return tmp_res
        elif isinstance(tmp_res, CompoundGroup):
            return ShnitselDBRoot(
                {tmp_res.name if tmp_res.name is not None else "Unknown": tmp_res}
            )
        elif isinstance(tmp_res, TrajectoryData) or isinstance(
            tmp_res, TrajectoryGroup
        ):
            tmp_compound = CompoundGroup(
                CompoundInfo(),
                {tmp_res.name if tmp_res.name is not None else "0": tmp_res},
            )
            return build_shnitsel_db(tmp_compound)

    raise ValueError(
        "Could not convert data to Shnitsel DataTree because data does not constitute a valid Shnitsel-style format."
    )


def convert_shnitsel_tree(
    data: xr.DataTree, restrict_levels=None
) -> ShnitselDBRoot | CompoundGroup | TrajectoryData | TrajectoryGroup:
    """Function to convert a generic xarray.DataTree into a Shnitsel-specific structure.

    Args:
        data (xr.DataTree): The generic DataTree to convert to Shnitsel-levels

    Raises:
        ValueError: If the "DataTree_Level" attribute is not set for a generic xr.DataTree node, we cannot convert.
        ValueError: If the "DataTree_Level" attribute has an invalid value which is not recognized by Shnitsel.

    Returns:
        ShnitselDBRoot|CompoundGroup|TrajectoryData|TrajectoryGroup: _description_
    """
    if (
        isinstance(data, TrajectoryData)
        or isinstance(data, TrajectoryGroup)
        or isinstance(data, CompoundGroup)
        or isinstance(data, ShnitselDBRoot)
    ):
        return data

    if _datatree_level_attribute_key not in data.attrs:
        raise ValueError("DataTree is not of valid ShnitselDB format.")

    if (
        restrict_levels is not None
        and data.attrs[_datatree_level_attribute_key] not in restrict_levels
    ):
        raise ValueError(
            f"Invalid level in hierarchy of ShnitselDB: found {data.attrs[_datatree_level_attribute_key]} but expected {restrict_levels}."
        )

    db_level = data.attrs[_datatree_level_attribute_key]
    if db_level == "ShnitselDBRoot":
        children_c: Mapping[str, CompoundGroup] = {
            k: convert_shnitsel_tree(v, restrict_levels=["CompoundGroup"])
            for k, v in data.children.items()
        }  # type: ignore
        res = ShnitselDBRoot(children_c)
    elif db_level == "CompoundGroup":
        if "compound_info" not in data.attrs:
            raise ValueError("Compound group level node has no compound_info set.")

        children_t: Mapping[str, TrajectoryGroup | TrajectoryData] = {
            k: convert_shnitsel_tree(
                v, restrict_levels=["TrajectoryGroup", "TrajectoryData"]
            )
            for k, v in data.children.items()
        }  # type: ignore
        res = CompoundGroup(CompoundInfo(**data.attrs["compound_info"]), children_t)
    elif db_level == "TrajectoryGroup":
        if "group_info" not in data.attrs:
            raise ValueError("TrajectoryGroup level node has no group_info set.")

        children_tg: Mapping[str, TrajectoryGroup | TrajectoryData] = {
            k: convert_shnitsel_tree(
                v, restrict_levels=["TrajectoryGroup", "TrajectoryData"]
            )
            for k, v in data.children.items()
        }  # type: ignore
        res = TrajectoryGroup(
            GroupInfo(
                data.name,
                data.attrs["group_info"] if "group_info" in data.attrs else None,
            ),
            children_tg,
        )
    elif db_level == "TrajectoryData":
        dataset = data.dataset
        res = TrajectoryData(dataset, name=data.name)
    else:
        logging.error(
            f"Unknown Level type in construction of Shnitsel DataTree: {data.attrs[_datatree_level_attribute_key]}."
        )
        raise ValueError(
            f"Unknown Shnitsel DataTree level: {data.attrs[_datatree_level_attribute_key]}"
        )

    return res


ShnitselDB: TypeAlias = ShnitselDBRoot
