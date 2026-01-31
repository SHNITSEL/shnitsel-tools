import logging
from typing import Callable, Dict, Mapping
from pandas import MultiIndex
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
            logging.info(
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

        # print("Before:", from_, " ->", to, "(", dividend, "/", divisor, ")" ":")
        # print(da)
        with xr.set_options(keep_attrs=True):
            res: xr.DataArray = da * dividend / divisor

        res.attrs.update({'units': to})
        if 'original_units' not in res.attrs:
            # Set an indicator for the original units of this array
            res.attrs.update({'original_units': from_})

        # print("After:", from_, " ->", to, "(", dividend, "/", divisor, ")" ":")
        # print(res)

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

# Helper to convert charges
convert_charge = Converter(
    definitions.unit_dimensions.charge, definitions._charge_unit_scales
)

# Helper to convert velocities
convert_velocity = Converter(
    definitions.unit_dimensions.velocity, definitions._velocity_unit_scales
)


def convert_all_units_to_shnitsel_defaults(data: xr.Dataset) -> xr.Dataset:
    """Helper function to convert all variables and coordinates with units in the dataset to
    shnitsel defaults.

    Parameters
    ----------
    data : xr.Dataset
        The dataset to convert to defaults.

    Returns
    -------
    xr.Dataset
        The converted dataset.

    Raises
    ------
    ValueError
        If conversion is performed on a multi-index, an error may occur if the index does not support the replacement of its levels.
    """
    new_vars = {}

    if "time" in data:
        assert "units" in data["time"].attrs, (
            "Dataset is missing `units` attribute on `time` coordinate"
        )
        time_unit = data["time"].attrs["units"]
    else:
        logging.warning(
            f"Missing `time` coordinate on input dataset. Available variables: {[str(x) for x in data.variables.keys()]}"
        )
        time_unit = None

    with xr.set_options(keep_attrs=True):
        for var_name in data.data_vars:
            if 'unitdim' in data[var_name].attrs:
                conv_res = convert_datarray_with_unitdim_to_shnitsel_defaults(
                    data[var_name]
                )

                logging.debug(
                    "Converting %s from unit %s to %s",
                    var_name,
                    data[var_name].attrs['units'],
                    conv_res.attrs['units'],
                )

                if var_name in data.indexes:
                    var_index = data.indexes[var_name]
                    if isinstance(var_index, MultiIndex):
                        raise ValueError(
                            "We do not support MultiIndices on non-coordinate variables."
                        )

                new_vars[var_name] = conv_res

    # logging.debug("Converting Data: " + str(list(new_vars.keys())))
    # NOTE: For some reason, sometimes, assigning multiple variables at once resulted in all of them being filled with NaN values.
    # NOTE: It may be an issue of setting the coordinate "time" before setting the variables. Split setting variables and coordinates
    tmp = data.assign(new_vars)

    new_coords = {}

    with xr.set_options(keep_attrs=True):
        for coord_name in data.coords:
            if 'unitdim' in data[coord_name].attrs:
                conv_res = convert_datarray_with_unitdim_to_shnitsel_defaults(
                    data[coord_name]
                )

                if coord_name in data.indexes:
                    coord_index = data.indexes[coord_name]
                    if isinstance(coord_index, MultiIndex):
                        from shnitsel.data.multi_indices import assign_levels

                        tmp = assign_levels(tmp, {str(coord_name): conv_res})
                        continue

                logging.debug(
                    "Converting coordinate %s from unit %s to %s",
                    coord_name,
                    data[coord_name].attrs['units'],
                    conv_res.attrs['units'],
                )

                new_coords[coord_name] = conv_res

    # logging.debug("Converting Coords: " + str(list(new_coords.keys())))
    # NOTE: Alignment screws us over if we convert the time before assigning the other variables.
    tmp = tmp.assign_coords(new_coords)

    if time_unit is not None:
        if "delta_t" in tmp.attrs:
            tmp.attrs["delta_t"] = convert_time.convert_value(
                tmp.attrs["delta_t"],
                convert_from=time_unit,
                to=tmp["time"].attrs["units"],
            )
            logging.debug(
                "Converting attribute %s from unit %s to %s",
                "delta_t",
                time_unit,
                tmp["time"].attrs["units"],
            )
        if "t_max" in tmp.attrs:
            tmp.attrs["t_max"] = convert_time.convert_value(
                tmp.attrs["t_max"],
                convert_from=time_unit,
                to=tmp["time"].attrs["units"],
            )
            logging.debug(
                "Converting attribute %s from unit %s to %s",
                "t_max",
                time_unit,
                tmp["time"].attrs["units"],
            )

    return tmp


def convert_to_target_units(
    data: xr.Dataset, unit_map: str | Mapping[str, str | None] | None
) -> xr.Dataset:
    """Helper function to convert units of variables and coordinates on a dataset to
    target units.

    If no target units are specified, will convert to the units set as default by shnitsel tools.

    Parameters
    ----------
    data : xr.Dataset
        The dataset to perform the conversion on.
    units : str | Mapping[str, str] | None
        Either a unit to convert all variables to or a mapping between variables that should be converted and their target unit.
        Set a target unit of `None` to convert to shnitsel default.

    Returns
    -------
    xr.Dataset
        The resulting dataset after conversion.

    Raises
    ------
    ValueError
        If conversion of a multi-index level fails.
    ValueError
        If a global target unit could not be applied to all arrays in the Dataset.
    """
    if unit_map is None:
        return convert_all_units_to_shnitsel_defaults(data)
    new_vars = {}

    if "time" in data:
        assert "units" in data["time"].attrs, (
            "Dataset is missing `units` attribute on `time` coordinate"
        )
        time_unit = data["time"].attrs["units"]
    else:
        logging.warning(
            f"Missing `time` coordinate on input dataset. Available variables: {[str(x) for x in data.variables.keys()]}"
        )
        time_unit = None

    global_conversion = isinstance(unit_map, str)
    global_unit = unit_map if global_conversion else None

    with xr.set_options(keep_attrs=True):
        for var_name in data.data_vars:
            if (global_conversion or var_name in unit_map) and 'unitdim' in data[
                var_name
            ].attrs:
                conv_res = convert_datarray_with_unitdim(
                    data[var_name],
                    global_unit if global_conversion else unit_map[str(var_name)],
                )

                logging.debug(
                    "Converting %s from unit %s to %s",
                    var_name,
                    data[var_name].attrs['units'],
                    conv_res.attrs['units'],
                )

                if var_name in data.indexes:
                    var_index = data.indexes[var_name]
                    if isinstance(var_index, MultiIndex):
                        raise ValueError(
                            "We do not support MultiIndices on non-coordinate variables."
                        )

                new_vars[var_name] = conv_res

    # logging.debug("Converting Data: " + str(list(new_vars.keys())))
    # NOTE: For some reason, sometimes, assigning multiple variables at once resulted in all of them being filled with NaN values.
    # NOTE: It may be an issue of setting the coordinate "time" before setting the variables. Split setting variables and coordinates
    tmp = data.assign(new_vars)

    new_coords = {}

    with xr.set_options(keep_attrs=True):
        for coord_name in data.coords:
            if (global_conversion or coord_name in unit_map) and 'unitdim' in data[
                coord_name
            ].attrs:
                conv_res = convert_datarray_with_unitdim(
                    data[coord_name],
                    global_unit if global_conversion else unit_map[str(coord_name)],
                )

                if coord_name in data.indexes:
                    coord_index = data.indexes[coord_name]
                    if isinstance(coord_index, MultiIndex):
                        from shnitsel.data.multi_indices import assign_levels

                        tmp = assign_levels(tmp, {str(coord_name): conv_res})
                        continue

                logging.debug(
                    "Converting coordinate %s from unit %s to %s",
                    coord_name,
                    data[coord_name].attrs['units'],
                    conv_res.attrs['units'],
                )

                new_coords[coord_name] = conv_res

    # logging.debug("Converting Coords: " + str(list(new_coords.keys())))
    # NOTE: Alignment screws us over if we convert the time before assigning the other variables.
    tmp = tmp.assign_coords(new_coords)

    if time_unit is not None:
        target_time_unit = str(
            (global_unit if global_conversion else unit_map["delta_t"])
            or tmp["time"].attrs["units"]
        )
        if "delta_t" in tmp.attrs and (
            global_conversion or "delta_t" in unit_map or 'time' in tmp
        ):
            tmp.attrs["delta_t"] = convert_time.convert_value(
                tmp.attrs["delta_t"],
                convert_from=time_unit,
                to=target_time_unit,
            )
            logging.debug(
                "Converting attribute %s from unit %s to %s",
                "delta_t",
                time_unit,
                target_time_unit,
            )
        if "t_max" in tmp.attrsand(
            global_conversion or "t_max" in unit_map or 'time' in tmp
        ):
            tmp.attrs["t_max"] = convert_time.convert_value(
                tmp.attrs["t_max"],
                convert_from=time_unit,
                to=target_time_unit,
            )
            logging.debug(
                "Converting attribute %s from unit %s to %s",
                "t_max",
                time_unit,
                target_time_unit,
            )

    return tmp


_CONVERTERS: Dict[str, Callable[[xr.DataArray, str], xr.DataArray]] = {
    definitions.unit_dimensions.energy: convert_energy,
    definitions.unit_dimensions.force: convert_force,
    definitions.unit_dimensions.dipole: convert_dipole,
    definitions.unit_dimensions.length: convert_length,
    definitions.unit_dimensions.time: convert_time,
    definitions.unit_dimensions.nacs: convert_nacs,
    definitions.unit_dimensions.socs: convert_socs,
    definitions.unit_dimensions.charge: convert_charge,
    definitions.unit_dimensions.velocity: convert_velocity,
}


def convert_datarray_with_unitdim_to_shnitsel_defaults(
    data: xr.DataArray,
) -> xr.DataArray:
    """Helper function to convert a data array to default values for
    shnitsel units.

    Parameters
    ----------
    data : xr.DataArray
        The array to convert

    Returns
    -------
    xr.DataArray
        The converted or untouched array if no conversion was necessary
    """
    return convert_datarray_with_unitdim(data, None)


def convert_datarray_with_unitdim(
    data: xr.DataArray, target_unit: str | None = None
) -> xr.DataArray:
    """Helper function to convert a dataarray with a unit dimension to
    either a target unit or a shnitsel default unit.

    Parameters
    ----------
    data : xr.DataArray
        The array to convert
    target_unit : str | None, optional
        The target unit, by default None which applies shnitsel default units

    Returns
    -------
    xr.DataArray
        The converted or untouched array if no conversion was necessary
    """
    if 'unitdim' in data.attrs:
        unit_dimension = data.attrs['unitdim']
        if unit_dimension in _CONVERTERS:
            if target_unit is None:
                if unit_dimension in definitions.standard_shnitsel_units:
                    return _CONVERTERS[unit_dimension](
                        data, definitions.standard_shnitsel_units[unit_dimension]
                    )
            else:
                return _CONVERTERS[unit_dimension](data, target_unit)
    return data
