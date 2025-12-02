from typing import NamedTuple
from typing_extensions import Literal

import xarray as xr

from .db_trajectory_group import GroupInfo, TrajectoryGroup


class ShMetadataKeys(NamedTuple):
    """Helper class to retain certain key metadata to differentiate different trajectories"""

    delta_t: float
    software: str
    version: str
    est_level: str | None
    theory_basis: str | None
    has_forces: bool | Literal['all', 'active_only'] | None
    input_type: Literal['static', 'dynamic', 'unknown']

    def to_groupname(self):
        return f"s:{self.software}|v{self.version}({self.input_type},{self.est_level}/{self.theory_basis})->dt={self.delta_t}[forces:{self.has_forces}]"


def group_subtree_by_metadata(
    subtree: xr.DataTree,
) -> xr.DataTree | tuple[ShMetadataKeys, xr.DataTree]:
    """If there are Trajectories below this level, introduces new groups, grouping them by specific metadata except when all of them already have the same metadata.

    Args:
        subtree (xr.DataTree): The subtree for which to recursively perform the grouping.

    Returns:
        xr.DataTree | tuple[ShMetadataKeys, xr.DataTree]: Either the grouped subtree or a pair of the metadata of the subtree and the subtree to be grouped in the Trajectory group of that metadata combination.
    """
    if len(subtree.children) > 0:
        res_children: dict[str, xr.DataTree] = {}

        meta_groups: dict[ShMetadataKeys, dict[str, xr.DataTree]] = {}
        grouped_keys = []

        for key, child in subtree.children.items():
            child_res = group_subtree_by_metadata(child)

            if isinstance(child_res, xr.DataTree):
                res_children[key] = child_res
            else:
                meta, new_child = child_res

                if meta not in meta_groups:
                    meta_groups[meta] = {}

                meta_groups[meta][key] = child
                grouped_keys.append(key)

        res_tree = subtree.copy()

        if len(res_children) == 0 and len(meta_groups) == 1:
            # Unwrap all grouped children again if we don't have other groups and children have same metadata
            for meta, children in meta_groups.items():
                res_children = children
                # Store new group info in the parent
                if "group_info" in res_tree.attrs:
                    if "group_attributes" in res_tree.attrs["group_info"]:
                        res_tree.attrs["group_info"]["group_attributes"].update(
                            meta._asdict()
                        )
                    else:
                        res_tree.attrs["group_info"]["group_attributes"] = (
                            meta._asdict()
                        )
                else:
                    groupname = "grouped:" + meta.to_groupname()
                    res_groupinfo = GroupInfo(
                        group_name=groupname, group_attributes=meta._asdict()
                    )
                    res_tree.attrs["group_info"] = res_groupinfo.__dict__
        else:
            for meta, children in meta_groups.items():
                groupname = "grouped:" + meta.to_groupname()
                res_groupinfo = GroupInfo(
                    group_name=groupname, group_attributes=meta._asdict()
                )
                new_grouping = TrajectoryGroup(res_groupinfo, children=children)
                res_children[groupname] = new_grouping

        res_tree = res_tree.drop_nodes(grouped_keys)
        res_tree = res_tree.assign(res_children)

        return res_tree
    else:
        if subtree.has_data:
            metadata = ShMetadataKeys(
                delta_t=subtree.dataset.attrs.get("delta_t", -1),
                software=subtree.dataset.attrs.get("input_format", "unknown"),
                input_type=subtree.dataset.attrs.get("input_type", "unknkown"),
                version=subtree.dataset.attrs.get("input_format_version", "unknown"),
                est_level=subtree.dataset.attrs.get("est_level", "unknown"),
                theory_basis=subtree.dataset.attrs.get("theory_basis", "unknown"),
                has_forces=subtree.dataset.attrs.get("has_forces", None),
            )

            return metadata, subtree.copy()
        else:
            return subtree.copy()
