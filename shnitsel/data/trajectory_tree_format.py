from dataclasses import dataclass
import logging
from typing import Any, Dict, List, Mapping, Self, TypeAlias
import xarray as xr

from shnitsel.data.trajectory_format import Trajectory


@dataclass
class GroupInfo:
    """Class to hold auxiliaryt info of a group of Trajectories in ShnitselDB"""
    grouped_settings: Dict[str, Any]


@dataclass
class CompoundInfo:
    """Class to hold identifying and auxiliary info of a compound type in ShnitselDB"""

    compound_name: str = "unknown"
    smile_repr: str | None = None


_datatree_level_attribute_key = "DataTree_Level"


class TrajectoryData(xr.DataTree):
    """DataTree node to keep track of a single trajectory entry"""

    def __init__(
        self,
        dataset: xr.Dataset | Trajectory | None = None,
        name: str | None = None,
    ):
        super().__init__(dataset=dataset, children=None, name=name)

        self.attrs[_datatree_level_attribute_key] = "TrajectoryData"


class TrajectoryGroup(xr.DataTree):
    """DataTree node to keep track of a group of trajectories where properties defining the group can be set"""

    def __init__(
        self,
        group_info: GroupInfo,
        children: Mapping[str, TrajectoryData | Self] | None = None,
        name: str | None = None,
    ):
        super().__init__(None, children, name)

        self.attrs[_datatree_level_attribute_key] = "TrajectoryGroup"
        self.attrs["group_info"] = group_info


class CompoundGroup(xr.DataTree):
    """DataTree node to keep track of all data associated with a common compound within the datatree"""

    def __init__(
        self,
        coumpound_info: CompoundInfo,
        children: Mapping[
            str,
            TrajectoryGroup | TrajectoryData,
        ]
        | None = None,
    ):
        super().__init__(None, children, coumpound_info.compound_name)
        self.attrs[_datatree_level_attribute_key] = "CompoundGroup"
        self.attrs["compound_info"] = coumpound_info


class ShnitselDBRoot(xr.DataTree):
    """DataTree root node to keep track of Shnitsel data in a predefined structural hierarchy"""

    def __init__(self, compounds: Mapping[str, CompoundGroup] | None = None):
        super().__init__(dataset=None, children=compounds, name="ROOT")
        self.attrs[_datatree_level_attribute_key] = "ShnitselDBRoot"


ShnitselDB: TypeAlias = ShnitselDBRoot


def build_shnitsel_db(
    data: Trajectory | xr.DataTree | List[Trajectory | xr.DataTree],
) -> ShnitselDB:
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
        ShnitselDBRoot: _description_
    """
    if isinstance(data, ShnitselDBRoot):
        return data
    elif isinstance(data, list):
        res_set = {}

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
                            if d_conv.dataset is not None
                            and "trajid" in d_conv.dataset.attrs
                            else f"{i}"
                        )
                    )
                else:
                    raise ValueError(f"Found unsupported type: {type(d_conv)}")

                res_set[key] = d_conv
            if shared_level_id == "TrajectoryData|TrajectoryGroup":
                tmp_compound = CompoundGroup(CompoundInfo(), children=res_set)
                return build_shnitsel_db(tmp_compound)
            else:
                return ShnitselDBRoot(res_set)
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
            tmp_compound_name = "unknown"
            tmp_compound = CompoundGroup(
                CompoundInfo(tmp_compound_name, None),
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
        raise ValueError(f"DataTree is not of valid ShnitselDB format.")

    if (
        restrict_levels is not None
        and data.attrs[_datatree_level_attribute_key] not in restrict_levels
    ):
        raise ValueError(
            f"Invalid level in hierarchy of ShnitselDB: found {data.attrs[_datatree_level_attribute_key]} but expected {restrict_levels}."
        )

    match data.attrs[_datatree_level_attribute_key]:
        case "ShnitselDBRoot":
            children_c: Mapping[str, CompoundGroup] = {
                k: convert_shnitsel_tree(v, restrict_levels=["CompoundGroup"])
                for k, v in data.children.items()
            }
            res = ShnitselDBRoot(children_c)
        case "CompoundGroup":
            if not "compound_info" in data.attrs:
                raise ValueError(f"Compound group level node has no compound_info set.")

            children_t: Mapping[str, TrajectoryGroup | TrajectoryData] = {
                k: convert_shnitsel_tree(
                    v, restrict_levels=["TrajectoryGroup", "TrajectoryData"]
                )
                for k, v in data.children.items()
            }
            res = CompoundGroup(data.attrs["compound_info"], children_t)
        case "TrajectoryGroup":
            if not "group_info" in data.attrs:
                raise ValueError(f"TrajectoryGroup level node has no group_info set.")

            children_tg: Mapping[str, TrajectoryGroup | TrajectoryData] = {
                k: convert_shnitsel_tree(
                    v, restrict_levels=["TrajectoryGroup", "TrajectoryData"]
                )
                for k, v in data.children.items()
            }
            res = TrajectoryGroup(data.attrs["group_info"], children_tg, name=data.name)
        case "TrajectoryData":
            dataset = data.dataset
            res = TrajectoryData(dataset, name=data.name)
        case _:
            logging.error(
                f"Unknown Level type in construction of Shnitsel DataTree: {data.attrs[_datatree_level_attribute_key]}."
            )
            raise ValueError(
                f"Unknown Shnitsel DataTree level: {data.attrs[_datatree_level_attribute_key]}"
            )

    return res
