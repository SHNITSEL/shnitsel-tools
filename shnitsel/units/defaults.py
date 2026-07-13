from typing import Dict, Literal
from ..io.shared.helpers import LoadingParameters
from .definitions import standard_units_of_formats, unit_dimensions
import xarray as xr
import numpy as np

DUMMY_INT_FILL_VALUE = -2e63 + 1


def get_fill_value(da: xr.DataArray | xr.Variable | np.ndarray):
    """Helper function to obtain either the set or the default fill value
    for a certain array dtype.

    Parameters
    ----------
    da : xr.DataArray | np.ndarray
        The array to get the fill value for. Either from the `fill_value` attribute or based on the dtype

    Returns
    ----------
    The value used as a fill value for the array.
    """
    tmp_res = None
    if isinstance(da, (xr.DataArray, xr.Variable)):
        tmp_res = da.attrs.get("fill_value", None)

    if tmp_res is None:
        global DUMMY_INT_FILL_VALUE
        dtype = da.dtype
        return DUMMY_INT_FILL_VALUE if np.issubdtype(dtype, np.integer) else np.nan
    return np.nan


def get_default_input_attributes(
    kind: Literal["sharc", "newtonx", "ase", "pyrai2md", "shnitsel"],
    loading_parameters: LoadingParameters | None = None,
) -> Dict[str, Dict[str, str]]:
    """Function to get the default attribute setup to read input from a certain file format.

    Used to set descriptions, long names, unit dimensions and unit names.

    Parameters
    ----------
    kind : Literal["sharc", "newtonx", "ase", "pyrai2md"]
        The kind of input format to get default settings for
    loading_parameters : LoadingParameters | None, optional
        User-provided overrides for default setup. Defaults to None.

    Returns
    -------
    dict[str, dict[str, str]]
        The resulting set of attributes for each individual supported observable in a dataset.
    """
    format_default_units = standard_units_of_formats[kind]

    def override_defaults(unit_dimension, variable_name):
        if (
            loading_parameters is not None
            and loading_parameters.input_units is not None
            and variable_name in loading_parameters.input_units
        ):
            return loading_parameters.input_units[variable_name]
        else:
            return format_default_units[unit_dimension]

    res = {
        "atXYZ": {
            "long_name": "Positions",
            "unitdim": unit_dimensions.length,
            "units": override_defaults(unit_dimensions.length, "atXYZ"),
        },
        "energy": {
            "long_name": "Absolute energy",
            "unitdim": unit_dimensions.energy,
            "units": override_defaults(unit_dimensions.energy, "energy"),
        },
        "e_kin": {
            "long_name": "Kinetic_energy",
            "unitdim": unit_dimensions.energy,
            "units": override_defaults(unit_dimensions.energy, "e_kin"),
        },
        "dip_all": {
            "long_name": "Complete dipoles",
            "unitdim": unit_dimensions.dipole,
            "units": override_defaults(unit_dimensions.dipole, "dip_all"),
        },
        "dip_perm": {
            "long_name": "Permanent dipoles",
            "unitdim": unit_dimensions.dipole,
            "units": override_defaults(unit_dimensions.dipole, "dip_perm"),
        },
        "dip_trans": {
            "long_name": "Transitional dipoles",
            "unitdim": unit_dimensions.dipole,
            "units": override_defaults(unit_dimensions.dipole, "dip_trans"),
        },
        "time": {
            "long_name": "Time in trajectory or timestep",
            "unitdim": unit_dimensions.time,
            "units": override_defaults(unit_dimensions.time, "time"),
        },
        "phases": {"long_name": "Phase vector"},
        "astate_diag": {
            "long_name": "Active state (diag)",
            "fill_value": -1,
        },
        "astate": {
            "long_name": "Active state (MCH)",
            "fill_value": -1,
        },
        "state": {"long_name": "Index of relevant states for indexing in MCH"},
        "state2": {"long_name": "The second state to build state combinations out of."},
        "state_diag": {"long_name": "State index in diagonal basis."},
        "state_names": {"long_name": "String representations of the states."},
        "state_types": {
            "long_name": "Multiplicity to indicate whether the respective state is singlet (1), doublet (2), or triplet(3)",
            "fill_value": -1,
        },
        "state_charges": {
            "long_name": "Charge of the various states.",
            "unitdim": unit_dimensions.charge,
            "units": override_defaults(unit_dimensions.charge, "state_charge"),
        },
        "statecomb": {
            "long_name": "Combination of two states used to index inter-state properties that don't depend on state order"
        },
        "frame": {
            "long_name": "An index enumerating all momentous frames in a set of combined trajectory data"
        },
        "atrajectory": {
            "long_name": "An index in a multi-trajectory dataset to specify, from which original trajectory this entry was merged. (active trajectory)"
        },
        "trajectory": {
            "long_name": "Index for multi-trajectory datasets for addressing per-trajectory properties"
        },
        "from": {
            "long_name": "An alias for the first state of a statecomb combination"
        },
        "to": {"long_name": "An alias for the second state of a statecomb combination"},
        "full_statecomb": {
            "long_name": "Combination of two states used to index inter-state properties that do depend on the order of states"
        },
        "full_statecomb_from": {
            "long_name": "An alias for the first state of a full_statecomb combination"
        },
        "full_statecomb_to": {
            "long_name": "An alias for the second state of a full_statecomb combination"
        },
        "atNames": {"long_name": "Names of atomic elements (short form)"},
        "atNums": {
            "long_name": "Periodic number of atomic elements",
            "fill_value": -1,
        },
        "forces": {
            "long_name": "Per-atom forces",
            "unitdim": unit_dimensions.force,
            "units": override_defaults(unit_dimensions.force, "forces"),
        },
        "nacs": {
            "long_name": "Nonadiabatic couplings",
            "unitdim": unit_dimensions.nacs,
            "units": override_defaults(unit_dimensions.nacs, "nacs"),
        },
        "socs": {
            "long_name": "Spin-orbit couplings",
            "unitdim": unit_dimensions.socs,
            "units": override_defaults(unit_dimensions.socs, "socs"),
        },
        "velocities": {
            "long_name": "Velocities of the atoms",
            "unitdim": unit_dimensions.velocity,
            "units": override_defaults(unit_dimensions.velocity, "velocities"),
        },
        "state_coefs_diag": {
            "long_name": "Coefficients of current active state in diagonal basis",
            "units": "1",
        },
        "prob_hop_diag": {
            "long_name": "Hopping probability to states in diagonal basis",
            "description": "Hopping probabilities denote the probability of hopping from the currently active diagonal state to the respective target diagonal state.",
            "units": "1",
        },
        "u_matrix": {
            "long_name": "Diag to MCH transformation matrix",
            "units": "1",
            "description": "The U matrix converts from the diag into the MCH matrix like $c_{MCH} = U c_{diag}$. "
            "As such, the 'state_diag' dimension should be considered as the index in the diagonal basis and the 'state' dimension as the index in the MCH basis.",
        },
        # history
        # A global attribute for an audit trail. This is a character array with a line for each invocation of a program that has modified the dataset. Well-behaved generic netCDF applications should append a line containing: date, time of day, user name, program name and command arguments.
    }

    return res
