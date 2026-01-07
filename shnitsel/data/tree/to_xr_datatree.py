import logging
from typing import Any, TypeVar
import xarray as xr
from .tree import ShnitselDBRoot
from .compound import CompoundGroup
from .data_group import DataGroup
from .data_leaf import DataLeaf
from ..dataset_containers.trajectory import Trajectory
from ..dataset_containers.frames import Frames
from dataclasses import asdict

from .datatree_level import _datatree_level_attribute_key

DataType = TypeVar("DataType")


def data_to_xarray_dataset(
    raw_data: xr.DataArray | xr.Dataset | Trajectory | Frames | None | Any,
    metadata: dict[str, Any],
) -> tuple[xr.Dataset | None, dict[str, Any]]:
    """Support function to convert data in trees (or that could be in trees) to xarray
    datasets.

    Parameters
    ----------
    raw_data : xr.DataArray | xr.Dataset | Trajectory | Frames | None
        Any type of raw data that should be converted. See type hints for supported types.
        Everything caught by `Any` and nothing else will be ignored and results in a `ValueError` being raised.
    metadata : dict[str, Any]
        Metadata attributes to update according to the performed conversion.

    Returns
    -------
    xr.Dataset | None
        The converted dataset result or none if no conversion was possible
    dict[str, Any]
        The updated metadata dictionary

    Raises
    ------
    ValueError
        If an unsupported type is provided for conversion.
    """
    if raw_data is None:
        tree_data = None
    elif isinstance(raw_data, xr.DataArray):
        metadata["_shnitsel_data_type"] = 'xarray::DataAray'
        metadata["_shnitsel_data_var"] = 'data'
        tree_data = raw_data.to_dataset(name='data')
    elif isinstance(raw_data, xr.Dataset):
        metadata["_shnitsel_data_type"] = 'xarray::Dataset'
        tree_data = raw_data
    elif isinstance(raw_data, Trajectory):
        metadata["_shnitsel_data_type"] = 'shnitsel::Trajectory'
        tree_data = raw_data.dataset
    elif isinstance(raw_data, Frames):
        metadata["_shnitsel_data_type"] = 'shnitsel::Frames'
        tree_data = raw_data.dataset
    else:
        logging.error(
            "Currently unsupported type %s found in tree to be converted to xarray datatree. Quietly skipping.",
            type(raw_data),
        )
        raise ValueError(
            "Currently unsupported type %s found in tree to be converted to xarray datatree. Quietly skipping."
            % type(raw_data)
        )
    return (tree_data, metadata)


def tree_to_xarray_datatree(
    node: ShnitselDBRoot[DataType]
    | CompoundGroup[DataType]
    | DataGroup[DataType]
    | DataLeaf[DataType],
) -> xr.DataTree | None:
    """Helper function to convert a ShnitselDB tree format to xarray.DataTree format
    so that we can use the xarray functions to write a netcdf file.

    Will recursively convert the tree from the current `node` starting from the leaves upwards.
    If the type of the node is not supported or the datatype in leaves is not supported for being stored via
    the xarray functions, the conversion will fail.

    Parameters
    ----------
    node : ShnitselDBRoot[DataType] | CompoundGroup[DataType] | DataGroup[DataType] | DataLeaf[DataType]
        The root node of a subtree to be converted to a `xr.DataTree` structure.

    Returns
    -------
    xr.DataTree | None
        Either the converted tree or None if this subtree is not supported.
    """
    node_attrs = dict(node.attrs)
    node_attrs[_datatree_level_attribute_key] = node._level_name
    node_attrs['_shnitsel_tree_indicator'] = "TREE"

    if isinstance(node, DataLeaf):
        tree_data = None
        if node.has_data:
            raw_data = node.data
            tree_data, node_attrs = data_to_xarray_dataset(
                raw_data=raw_data, metadata=node_attrs
            )

        res_tree = xr.DataTree(dataset=tree_data, name=node.name)
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
            node_attrs['_shnitsel_compound_info'] = asdict(compound_info)
            if node._group_info:
                node_attrs['_shnitsel_group_info'] = asdict(node._group_info)
        elif isinstance(node, DataGroup):
            if node._group_info:
                node_attrs['_shnitsel_group_info'] = asdict(node._group_info)
        else:
            logging.error(
                "Currently unsupported node type %s found in tree to be converted to xarray datatree. Quietly skipping.",
                type(node),
            )
            return None

        res_tree = xr.DataTree(dataset=None, children=children_as_trees, name=node.name)
        res_tree.attrs.update(node_attrs)
        return res_tree
