from typing import NamedTuple
import xarray as xr


class UnitInfo(NamedTuple):
    """Class to assign standard keys to the attributes of variables associated with dimensions"""

    unitdim: str
    units: str


class VariableInfoAttributes(NamedTuple):
    """Class to keep track of all variable attributes in xarrays acknowledged by the shnitsel standard"""

    long_name: str | None
    unit_info: UnitInfo | None


def set_variable_attributes(var: xr.DataArray, attrs: VariableInfoAttributes) -> None:
    """Function to set standardized variable attributes.

    This is meant to unify the keys under which information is stored in the shnitsel format.

    Parameters
    ----------
    var : xr.DataArray
        The variable/DataArray to set attributes on
    attrs : VariableInfoAttributes
        The standardized set of attributes to assign to the variable
    """

    if attrs.long_name is not None:
        var.attrs["long_name"] = attrs.long_name

    # TODO: We should add a check whether the units/unitdim combination is valid.
    if attrs.unit_info is not None:
        var.attrs["unitdim"] = attrs.unit_info.unitdim
        var.attrs["units"] = attrs.unit_info.units
