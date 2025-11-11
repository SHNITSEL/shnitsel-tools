import logging
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
        self.case_sensitive_mapping = {k.lower(): k for k in conversions.keys()}

        self.targets = list(conversions.keys())

    def __call__(
        self, da: xr.DataArray, to: str, convert_from: str | None = None
    ) -> xr.DataArray:
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

        if to == "1":
            logging.warning(
                f"Target is {to} for {da.name}, which means we do not care about the target unit or do not have a standard."
            )
            return da

        if convert_from is None:
            try:
                from_ = da.attrs['units']
            except (AttributeError, KeyError):
                raise KeyError(
                    "The 'units' attribute of the DataArray must be set and of type str."
                )
        else:
            from_ = convert_from

        if to.lower() == from_.lower():
            # If unit is the same, do not convert
            return da

        try:
            dividend = self.conversions[from_.lower()]
        except KeyError:
            raise ValueError(
                f"Can't convert {self.quantity_name} from {from_!r}, only from: {self.targets}"
            )

        try:
            divisor = self.conversions[to.lower()]
        except KeyError:
            raise ValueError(
                f"Can't convert {self.quantity_name} to {to!r}, only to: {self.targets}"
            )

        print(from_, "->", to, ": /", divisor, "*", dividend)

        with xr.set_options(keep_attrs=True):
            res: xr.DataArray = da * dividend / divisor

        res.attrs.update({'units': to})
        if 'original_units' not in res.attrs:
            # Set an indicator for the original units of this array
            res.attrs.update({'original_units': from_})

        return res

    def convert_value(
        self,
        value: float,
        convert_from: str,
        to: str,
    ) -> float:
        if to == "1":
            logging.warning(
                f"Target is {to}, which means we do not care about the target unit or do not have a standard."
            )
            return value

        from_ = convert_from

        if to.lower() == from_.lower():
            # If unit is the same, do not convert
            return value

        try:
            dividend = self.conversions[from_.lower()]
        except KeyError:
            raise ValueError(
                f"Can't convert {self.quantity_name} from {from_!r}, only from: {self.targets}"
            )

        try:
            divisor = self.conversions[to.lower()]
        except KeyError:
            raise ValueError(
                f"Can't convert {self.quantity_name} to {to!r}, only to: {self.targets}"
            )

        return value * dividend / divisor


# Helper to convert energies
convert_energy = Converter(
    definitions.unit_dimensions.energy, definitions._energy_unit_scales
)

# Helper to convert forces
convert_force = Converter(
    definitions.unit_dimensions.force, definitions._force_unit_scales
)

# Helper to convert dipole moments
convert_dipole = Converter(
    definitions.unit_dimensions.dipole, definitions._dipole_unit_scales
)

# Helper to convert lengths and distances
convert_length = Converter(
    definitions.unit_dimensions.length, definitions._distance_unit_scales
)

# Helper to convert time
convert_time = Converter(
    definitions.unit_dimensions.time, definitions._time_unit_scales
)

# Helper to convert nacs
convert_nacs = Converter(
    definitions.unit_dimensions.nacs, definitions._nacs_unit_scales
)

# Helper to convert socs
convert_socs = Converter(definitions.unit_dimensions.time, definitions._socs_unit_scale)


def convert_all_units_to_shnitsel_defaults(data: xr.Dataset) -> xr.Dataset:
    new_vars = {}

    time_unit = data["time"].attrs["units"]

    with xr.set_options(keep_attrs=True):
        for var_name in data.variables:
            if 'unitdim' in data[var_name].attrs:
                new_vars[var_name] = convert_datarray_with_unitdim_to_shnitsel_defaults(
                    data[var_name]
                )

    if "delta_t" in data.attrs:
        data.attrs["delta_t"] = convert_time.convert_value(
            data.attrs["delta_t"],
            convert_from=time_unit,
            to=new_vars["time"].attrs["units"],
        )
    if "t_max" in data.attrs:
        data.attrs["t_max"] = convert_time.convert_value(
            data.attrs["t_max"],
            convert_from=time_unit,
            to=new_vars["time"].attrs["units"],
        )

    logging.debug("Converting: " + str(list(new_vars.keys())))
    return data.assign(new_vars)


_CONVERTERS: Dict[str, Callable[[xr.DataArray, str], xr.DataArray]] = {
    definitions.unit_dimensions.energy: convert_energy,
    definitions.unit_dimensions.force: convert_force,
    definitions.unit_dimensions.dipole: convert_dipole,
    definitions.unit_dimensions.length: convert_length,
    definitions.unit_dimensions.time: convert_time,
    definitions.unit_dimensions.nacs: convert_nacs,
    definitions.unit_dimensions.socs: convert_socs,
}


def convert_datarray_with_unitdim_to_shnitsel_defaults(
    data: xr.DataArray,
) -> xr.DataArray:
    if 'unitdim' in data.attrs:
        unit_dimension = data.attrs['unitdim']

        if (
            unit_dimension in definitions.standard_shnitsel_units
            and unit_dimension in _CONVERTERS
        ):
            return _CONVERTERS[unit_dimension](
                data, definitions.standard_shnitsel_units[unit_dimension]
            )

    return data
