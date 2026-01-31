from typing import Dict, Literal
from ..io.shared.helpers import LoadingParameters
from .definitions import standard_units_of_formats, unit_dimensions


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
        "sdiag": {"long_name": "Active state (diag)"},
        "astate": {"long_name": "Active state in dynamic trajectories (MCH)"},
        "state": {"long_name": "Index of relevant states for indexing"},
        "state2": {"long_name": "The second state to build state combinations out of"},
        "state_names": {"long_name": "String representations of the states."},
        "state_types": {
            "long_name": "Multiplicity to indicate whether the respective state is singlet (1), doublet (2), or triplet(3)"
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
        "trajid": {
            "long_name": "An index in a multi-trajectory dataset to specify, from which original trajectory this entry was merged."
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
        "atNums": {"long_name": "Periodic number of atomic elements"},
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
    }

    return res
