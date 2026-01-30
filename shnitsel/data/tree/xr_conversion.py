import logging
from typing import Any, TypeVar
from typing_extensions import TypeForm
import xarray as xr
from shnitsel.data.helpers import dataclass_from_dict
from shnitsel.data.tree.node import TreeNode
from ..xr_io_compatibility import SupportsToXrConversion

from ..dataset_containers.xr_conversion import (
    data_to_xarray_dataset,
    xr_dataset_to_shnitsel_format,
)
from .tree import ShnitselDBRoot
from .compound import CompoundGroup, CompoundInfo
from .data_group import DataGroup, GroupInfo
from .data_leaf import DataLeaf
from dataclasses import asdict

from .datatree_level import _datatree_level_attribute_key, DataTreeLevelMap

DataType = TypeVar("DataType", bound=xr.Dataset | xr.DataArray | SupportsToXrConversion)


def tree_to_xarray_datatree(
    node: (TreeNode[Any, DataType]),
) -> xr.DataTree | None:
    """Helper function to convert a ShnitselDB tree format to xarray.DataTree format
    so that we can use the xarray functions to write a netcdf file.

    Will recursively convert the tree from the current `node` starting from the leaves upwards.
    If the type of the node is not supported or the datatype in leaves is not supported for being stored via
    the xarray functions, the conversion will fail.

    Parameters
    ----------
    node : TreeNode[Any, DataType]
        The root node of a subtree to be converted to a `xr.DataTree` structure.

    Returns
    -------
    xr.DataTree | None
        Either the converted tree or None if this subtree is not supported.
    """
    node_attrs = dict(node.attrs)
    node_attrs[_datatree_level_attribute_key] = node._level_name
    node_attrs["_shnitsel_tree_indicator"] = "TREE"
    node_name = str(node.name) if node.name is not None else None

    if isinstance(node, DataLeaf):
        tree_data = None
        if node.has_data:
            raw_data = node.data
            metadata = {}
            tree_data, metadata = data_to_xarray_dataset(
                raw_data=raw_data, metadata=metadata
            )
            node_attrs["_shnitsel_io_meta"] = metadata

        res_tree = xr.DataTree(dataset=tree_data, name=node_name)
        res_tree.attrs.update(node_attrs)
        return res_tree
    else:
        children_as_trees = {
            str(k): res
            for k, c in node.children.items()
            if (res := tree_to_xarray_datatree(c)) is not None
        }
        if isinstance(node, ShnitselDBRoot):
            # No further updates of properties for root node
            pass
        elif isinstance(node, CompoundGroup):
            compound_info = node.compound_info
            node_attrs["_shnitsel_compound_info"] = asdict(compound_info)
            if node._group_info:
                node_attrs["_shnitsel_group_info"] = asdict(node._group_info)
        elif isinstance(node, DataGroup):
            if node._group_info:
                node_attrs["_shnitsel_group_info"] = asdict(node._group_info)
        else:
            logging.error(
                "Currently unsupported node type %s found in tree to be converted to xarray datatree. Quietly skipping.",
                type(node),
            )
            return None

        res_tree = xr.DataTree(
            dataset=None,
            children=children_as_trees,
            name=node_name,
        )
        res_tree.attrs.update(node_attrs)
        return res_tree


def xarray_datatree_to_shnitsel_tree(
    node: xr.DataTree, dtype: type[DataType] | TypeForm[DataType] | None = None
) -> ShnitselDBRoot | CompoundGroup | DataGroup | DataLeaf | None:
    """
    Helper function to invert the operation of `tree_to_xarray_datastree` and deserialize the
    shnitsel tree/ShnitselDB from a stored xarray DataTree:

    Parameters
    ----------
    node : xr.DataTree
        The root node of a xarray subtree.
        will convert this subtree recursively.
    dtype: type[DataType] | TypeForm[DataType], optional
        Optional argument to specify the desired target type of data in the shnitsel tree structure.

    Returns
    -------
    ShnitselDBRoot | CompoundGroup | DataGroup | DataLeaf | None
        The converted type or `None` if the tree could not be converted.
    """
    # print(f"tree conversion {node=} {type(node)=}")
    # print(f"attrs: {node.attrs}")
    if (
        "_shnitsel_tree_indicator" not in node.attrs
        and _datatree_level_attribute_key not in node.attrs
    ):
        # Conversion of arbitrary tree without type hints.
        # We do not need to support this, but we will do our best.
        if node.has_data:
            assert len(node.children) == 0, (
                "no children must be provided at `data` level."
            )

            metadata = node.attrs.get("_shnitsel_io_meta", {})
            remaining_node_meta = dict(node.attrs)
            if "_shnitsel_io_meta" in remaining_node_meta:
                del remaining_node_meta["_shnitsel_io_meta"]
            return DataLeaf(
                name=node.name,
                data=xr_dataset_to_shnitsel_format(node.dataset, metadata),
                attrs=remaining_node_meta,
                dtype=dtype,
            )
        else:
            child_types = set()
            converted_children = {
                str(k): child_res
                for k, v in node.children.items()
                if (child_res := xarray_datatree_to_shnitsel_tree(v)) is not None
            }

            child_types = set(type(c) for c in converted_children.values())

            is_likely_root = node.parent is None and all(
                issubclass(ct, CompoundGroup) for ct in child_types
            )
            is_likely_compound = any(
                str(attr_name).find("compound") for attr_name in node.attrs
            ) and all(
                issubclass(ct, DataGroup) or issubclass(ct, DataLeaf)
                for ct in child_types
            )
            is_likely_group = any(
                str(attr_name).find("compound") for attr_name in node.attrs
            ) and all(
                issubclass(ct, DataGroup) or issubclass(ct, DataLeaf)
                for ct in child_types
            )
            is_likely_data = len(node.children) == 0

            if is_likely_root:
                return ShnitselDBRoot(
                    compounds=converted_children,  # type: ignore # The above assertion checks the correct type of the children
                    attrs=node.attrs,
                    dtype=dtype,
                )
            if is_likely_compound:
                compound_info = CompoundInfo()
                if "_shnitsel_compound_info" in node.attrs:
                    compound_info = dataclass_from_dict(
                        CompoundInfo, node.attrs["_shnitsel_compound_info"]
                    )
                elif "_compound_info" in node.attrs:
                    compound_info = dataclass_from_dict(
                        CompoundInfo, node.attrs["_compound_info"]
                    )
                group_info = None
                if "_shnitsel_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_shnitsel_group_info"]
                    )
                elif "_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_group_info"]
                    )

                return CompoundGroup(
                    name=node.name,
                    children=converted_children,  # type: ignore # The above assertion checks the correct type of the children
                    compound_info=compound_info,
                    group_info=group_info,
                    attrs={str(k): v for k, v in node.attrs.items()},
                    dtype=dtype,
                )

            if is_likely_group:
                group_info = None
                if "_shnitsel_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_shnitsel_group_info"]
                    )
                elif "_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_group_info"]
                    )

                return DataGroup(
                    name=node.name,
                    children=converted_children,  # type: ignore # The likelihood flag checks for correct child types
                    group_info=group_info,
                    attrs={str(k): v for k, v in node.attrs.items()},
                    dtype=dtype,
                )

            if is_likely_data:
                return DataLeaf(
                    name=node.name,
                    data=None,
                    attrs={str(k): v for k, v in node.attrs.items()},
                    dtype=dtype,
                )
        raise ValueError(
            "Provided tree did not have type hints and structure could not be mapped to shnitsel data structure."
        )
    else:
        # print("annotated tree")
        datatree_level = node.attrs.get(_datatree_level_attribute_key)
        if datatree_level not in DataTreeLevelMap:
            raise ValueError(
                "Unsupported DataTree level indicator %s. Tree was potentially generated with a later version of shnitsel tools."
                % datatree_level
            )
        mapped_level = DataTreeLevelMap[datatree_level]
        if mapped_level == DataTreeLevelMap["data"]:
            assert len(node.children) == 0, (
                "no children must be provided at `data` level."
            )

            metadata = node.attrs.get("_shnitsel_io_meta", {})
            remaining_node_meta = dict(node.attrs)
            if "_shnitsel_io_meta" in remaining_node_meta:
                del remaining_node_meta["_shnitsel_io_meta"]

            return DataLeaf(
                name=node.name,
                data=xr_dataset_to_shnitsel_format(node.dataset, metadata),
                attrs=remaining_node_meta,
                dtype=dtype,
            )
        else:
            converted_children = {
                str(k): child_res
                for k, v in node.children.items()
                if (child_res := xarray_datatree_to_shnitsel_tree(v)) is not None
            }

            if mapped_level == DataTreeLevelMap["root"]:
                assert all(
                    isinstance(child, CompoundGroup)
                    for child in converted_children.values()
                ), "Malformed tree provided as input for tree root."

                root_res = ShnitselDBRoot(
                    name=node.name,
                    compounds=converted_children,  # type: ignore # The above assertion checks the correct type of the children
                    attrs=node.attrs,
                    dtype=dtype,
                )
                return root_res

            if mapped_level == DataTreeLevelMap["compound"]:
                assert all(
                    isinstance(child, (DataGroup, DataLeaf))
                    for child in converted_children.values()
                ), "Malformed tree provided as input for compound."

                compound_info = CompoundInfo()
                if "_shnitsel_compound_info" in node.attrs:
                    compound_info = dataclass_from_dict(
                        CompoundInfo, node.attrs["_shnitsel_compound_info"]
                    )
                elif "_compound_info" in node.attrs:
                    compound_info = dataclass_from_dict(
                        CompoundInfo, node.attrs["_compound_info"]
                    )
                group_info = None
                if "_shnitsel_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_shnitsel_group_info"]
                    )
                elif "_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_group_info"]
                    )

                return CompoundGroup(
                    name=node.name,
                    children=converted_children,  # type: ignore # The above assertion checks the correct type of the children
                    compound_info=compound_info,
                    group_info=group_info,
                    attrs={str(k): v for k, v in node.attrs.items()},
                    dtype=dtype,
                )

            if mapped_level == DataTreeLevelMap["group"]:
                assert all(
                    isinstance(child, (DataGroup, DataLeaf))
                    for child in converted_children.values()
                ), "Malformed tree provided as input for group."

                group_info = None
                if "_shnitsel_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_shnitsel_group_info"]
                    )
                elif "_group_info" in node.attrs:
                    group_info = dataclass_from_dict(
                        GroupInfo, node.attrs["_group_info"]
                    )

                return DataGroup(
                    name=node.name,
                    children=converted_children,  # type: ignore # The above assertion checks the correct type of the children
                    group_info=group_info,
                    attrs={str(k): v for k, v in node.attrs.items()},
                    dtype=dtype,
                )
    logging.warning("Could not infer shnitsel node type for node: %s", node)
    return None
