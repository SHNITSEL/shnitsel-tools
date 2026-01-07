import logging
from typing import Any
from . import Trajectory, Frames, InterState, PerState
import xarray as xr


def data_to_xarray_dataset(
    raw_data: xr.DataArray
    | xr.Dataset
    | Trajectory
    | Frames
    | InterState
    | PerState
    | None
    | Any,
    metadata: dict[str, Any],
) -> tuple[xr.Dataset | None, dict[str, Any]]:
    """Support function to convert simple, non-hierarchical data structures (that could be in trees)
    to xarray datasets.

    Parameters
    ----------
    raw_data : xr.DataArray | xr.Dataset | Trajectory | Frames | InterState | PerState | None
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
    elif isinstance(raw_data, Frames):
        metadata["_shnitsel_data_type"] = 'shnitsel::Frames'
        tree_data = raw_data.dataset
    elif isinstance(raw_data, Trajectory):
        metadata["_shnitsel_data_type"] = 'shnitsel::Trajectory'
        tree_data = raw_data.dataset
    elif isinstance(raw_data, InterState):
        metadata["_shnitsel_data_type"] = 'shnitsel::InterState'
        tree_data = raw_data.dataset
    elif isinstance(raw_data, PerState):
        metadata["_shnitsel_data_type"] = 'shnitsel::PerState'
        tree_data = raw_data.dataset
    else:
        logging.error(
            "Currently unsupported type %s found in data to be converted to xarray dataset.",
            type(raw_data),
        )
        raise ValueError(
            "Currently unsupported type %s found in data to be converted to xarray dataset."
            % type(raw_data)
        )
    return (tree_data, metadata)


def xr_dataset_to_shnitsel_format(
    raw_data: xr.Dataset,
) -> Trajectory | Frames | InterState | PerState | xr.Dataset | xr.DataArray:
    if '_shnitsel_data_type' in raw_data.attrs:
        shnitsel_type_hint = raw_data.attrs["_shnitsel_data_type"]
        if shnitsel_type_hint == 'xarray::Dataset':
            return raw_data
        elif shnitsel_type_hint == 'xarray::DataAray':
            array_var_name = raw_data.attrs.get("_shnitsel_data_var", 'data')
            return raw_data[array_var_name]
        elif shnitsel_type_hint == 'shnitsel::Frames':
            return Frames(raw_data)
        elif shnitsel_type_hint == 'shnitsel::Trajectory':
            return Trajectory(raw_data)
        elif shnitsel_type_hint == 'shnitsel::Interstate':
            return InterState(direct_interstate_data=raw_data)
        elif shnitsel_type_hint == 'shnitsel::PerState':
            return PerState(direct_perstate_data=raw_data)
        else:
            raise ValueError(
                "Unknown shnitsel deserialization type: %s" % shnitsel_type_hint
            )
    else:
        # Make best effort guess for the type:
        if (
            'energy_interstate' in raw_data
            or 'dip_trans_norm' in raw_data
            or 'fosc' in raw_data
            or 'nacs_norm' in raw_data
            or 'socs_norm' in raw_data
        ):
            try:
                return InterState(direct_interstate_data=raw_data)
            except AssertionError:
                pass

        if 'forces_norm' in raw_data or 'dip_perm_norm' in raw_data:
            try:
                return PerState(direct_perstate_data=raw_data)
            except AssertionError:
                pass

        if 'frame' in raw_data.dims:
            try:
                return Frames(raw_data)
            except AssertionError:
                pass

        if 'time' in raw_data.dims:
            try:
                return Trajectory(raw_data)
            except AssertionError:
                pass

        # Did not find matching shnitsel type to wrap
        return raw_data
