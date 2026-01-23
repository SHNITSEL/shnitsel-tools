from html import escape
import xarray as xr
from .shared import ShnitselDataset
from xarray.core.formatting_html import (
    dim_section,
    coord_section,
    datavar_section,
    index_section,
    attr_section,
    _obj_repr,
)

from xarray.core.formatting import (
    filter_nondefault_indexes,
)
from xarray.core.options import OPTIONS, _get_boolean_with_default


def _get_indexes_dict(indexes):
    return {
        tuple(index_vars.keys()): idx for idx, index_vars in indexes.group_by_index()
    }


def shnitsel_dataset_repr(ds: ShnitselDataset) -> str:
    """Helper function to represent Shnitsel Datasets
    similarly to xarray representation but indicating the contents of the Dataset
    through inclusion of the specific wrapper type.

    Parameters
    ----------
    ds : ShnitselDataset
        The Shnitsel-tools wrapped dataset.

    Returns
    -------
    str
        The html representation of the shnitsel dataset wrapper.
    """
    obj_type = f"{type(ds).__name__}[xarray.{type(ds._raw_dataset).__name__}]"

    header_components = [f"<div class='xr-obj-type'>{escape(obj_type)}</div>"]

    sections = []

    sections.append(dim_section(ds))

    if ds.coords:
        sections.append(coord_section(ds.coords))

    sections.append(datavar_section(ds.data_vars))

    display_default_indexes = _get_boolean_with_default(
        "display_default_indexes", False
    )
    xindexes = filter_nondefault_indexes(
        _get_indexes_dict(ds.xindexes), not display_default_indexes
    )
    if xindexes:
        sections.append(index_section(xindexes))

    if ds.attrs:
        sections.append(attr_section(ds.attrs))

    return _obj_repr(ds, header_components, sections)
