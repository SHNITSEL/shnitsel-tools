from typing import Callable, Dict
import xarray as xr
import shnitsel.units.definitions as definitions


class Converter:
    """Implements a generic Callable object to convert DataArrays between different units.

    See documentation of the ``__call__`` method for details on the implementation of the conversion.
    """

    def __init__(self, quantity_name: str, conversions: Dict[str, float]):
        self.__name__ = f"convert_{quantity_name}"
        self.quantity_name = quantity_name
        # Convert all keys to lower case to avoid some capitalization issues
        self.conversions = {k.lower(): v for k, v in conversions.items()}
        # For debugging and setting the output unit, keep the original capitalization:
        self.case_sensitive_mapping = {
            k.lower(): k for k in conversions.keys()}

        self.targets = list(conversions.keys())

    def __call__(self, da: xr.DataArray, to: str, convert_from: str | None = None) -> xr.DataArray:
        """Function to convert an xr.DataArray between two different units. 

        The function needs to be provided with a target unit `to` as well as with either an attribute `units` of the 
        input data array or by setting the input unit `convert_from`.
        A new datarray with converted units will be returned, if successful.
        The resulting array will have the new unit set as its `units` attribute and if no prior conversion has happened, 
        the previous unit will be set to the `original_units` attribute.

        Args:
            da (xr.DataArray): The input DataArray whose data should be converted
            to (str): The unit to which the data should be converted
            convert_from (str | None, optional): The Unit from which conversion should be started if `da` does not have a `units` attribute. Defaults to None.

        Raises:
            KeyError: Raised if the convert_from parameter is not set and the `da` input has no `units` attribute set.
            ValueError: Raised if the target unit `to` is not known in the conversions dictionary.
            ValueError: Raised if the original unit `convert_from` or `da.attr['units']` is not known in the conversions dictionary.

        Returns:
            xr.DataArray: The array containing the unit-converted data out of `da` with a new `units` and potentially `original_units` attribute set.
        """
        if convert_from is None:
            try:
                from_ = da.attrs['units']
            except (AttributeError, KeyError):
                raise KeyError(
                    "The 'units' attribute of the DataArray must be set and of type str.")
        else:
            from_ = convert_from

        if to.lower() == from_.lower():
            # If unit is the same, do not convert
            return da

        try:
            divisor = self.conversions[from_.lower()]
        except KeyError:
            raise ValueError(
                f"Can't convert {self.quantity_name} from {from_!r}, only from: {self.targets}")

        try:
            dividend = self.conversions[to.lower()]
        except KeyError:
            raise ValueError(
                f"Can't convert {self.quantity_name} to {to!r}, only to: {self.targets}")

        with xr.set_options(keep_attrs=True):
            res: xr.DataArray = da * dividend / divisor

        res.attrs.update({'units': to})
        if 'original_units' not in res.attrs:
            # Set an indicator for the original units of this array
            res.attrs.update({'original_units': from_})

        return res


# Helper to convert energies
convert_energy = Converter(
    units.unit_dimensions.energy, units._energy_unit_scales)

# Helper to convert forces
convert_force = Converter(units.unit_dimensions.force,
                          units._force_unit_scales)

# Helper to convert dipole moments
convert_dipole = Converter(
    units.unit_dimensions.dipole, units._dipole_unit_scales)

# Helper to convert lengths and distances
convert_length = Converter(
    units.unit_dimensions.length, units._distance_unit_scales)

# Helper to convert time
convert_time = Converter(units.unit_dimensions.time, units._time_unit_scales)

# Helper to convert nacs
convert_nacs = Converter(units.unit_dimensions.nacs, units._nacs_unit_scales)

# Helper to convert socs
convert_socs = Converter(units.unit_dimensions.time, units._socs_unit_scale)
