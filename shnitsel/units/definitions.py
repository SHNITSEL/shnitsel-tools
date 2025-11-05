from typing import Callable, Dict, Literal
import ase.units as si
import numpy as np
from shnitsel.io.helpers import LoadingParameters
import xarray as xr

# TODO: FIXME: Check all atomic units for correctness


class unit_dimensions:
    length = "length"
    energy = "energy"
    force = "force"
    dipole = "dipole"
    nacs = "nacs"
    time = "time"
    soc = "soc"


class time:
    pico_seconds = "ps"
    femto_seconds = "fs"
    nano_seconds = "ns"
    seconds = "s"
    ase_time_units = "ase_time"
    au = "au"


_time_unit_scales = {
    time.pico_seconds: si.fs * 1e3,
    time.femto_seconds: si.fs,
    time.nano_seconds: si.fs * 1e6,
    time.seconds: si.second,
    "seconds": si.second * 1e-12,
    # TODO: FIXME: set actual atomic time unit possibly:  si._aut
    time.au: si._aut,
    # ASE uses this time unit: Ang/sqrt(u/eV) which might differ slightly from 10fs
    time.ase_time_units: si.Angstrom / np.sqrt(si._amu / si.eV),
}


class energy:
    Hartree = "Hartree"
    Eh = "Eh"
    eV = "eV"
    keV = "keV"
    J = "J"
    kJ = "kJ"
    kcal = "kcal"
    au = "au"


_energy_unit_scales = {
    energy.Hartree: si.Hartree,
    energy.Eh: si.Hartree,
    energy.eV: si.eV,
    energy.keV: si.eV * 1e3,
    "joule": si.J,
    energy.J: si.J,
    energy.kJ: si.kJ,
    energy.kcal: si.kcal,
    energy.au: si.Hartree,
}


class force:
    Hartree_per_Bohr = "Hartree/Bohr"
    Eh_per_Bohr = "Eh/Bohr"
    Hartree_per_Angstrom = "Hartree/Angstrom"
    Eh_per_Angstrom = "Eh/Angstrom"
    eV_per_Bohr = "eV/Bohr"
    eV_per_Angstrom = "eV/Angstrom"
    Newton = "N"
    au = "au"


_force_unit_scales = {
    force.Hartree_per_Bohr: si.Hartree / si.Bohr,
    force.Eh_per_Bohr: si.Hartree / si.Bohr,
    force.eV_per_Bohr: si.eV / si.Bohr,
    force.Hartree_per_Angstrom: si.Hartree / si.Angstrom,
    force.Eh_per_Angstrom: si.Hartree / si.Angstrom,
    force.eV_per_Angstrom: si.eV / si.Angstrom,
    "Hartree/A": si.Hartree / si.Angstrom,
    "eV/A": si.eV / si.Angstrom,
    force.Newton: si.J / si.m,
    force.au: si.Hartree / si.Bohr,
}


class distance:
    Bohr = "Bohr"
    Angstrom = "Angstrom"
    meter = "meter"
    pico_meter = "pm"
    nano_meter = "nm"
    au = "au"


_distance_unit_scales = {
    distance.Bohr: si.Bohr,
    distance.Angstrom: si.Angstrom,
    "A": si.Angstrom,
    distance.meter: si.m,
    distance.pico_meter: si.nm * 1e-3,
    distance.nano_meter: si.nm,
    "1": 1.0,
    distance.au: si.Bohr,
}

# An alias
length = distance
_length_unit_scales = _distance_unit_scales


class dipole:
    Debye = "debye"
    au = "au"


_dipole_unit_scales = {
    "1": 1.0,
    dipole.Debye: si.Debye,
    # TODO: FIXME: set actual atomic dipole unit.
    dipole.au: np.nan,
}


class nacs:
    # TODO: FIXME: Figure out nacs units
    one_per_Bohr = "one_per_Bohr"
    au = "au"


_nacs_unit_scales = {
    dipole.one_per_Bohr: 1.0,
    dipole.au: 1.0,
}


class socs:
    # TODO: FIXME: Figure out nacs units
    au = "au"


_socs_unit_scale = {
    socs.au: 1.0,
}


class amount_units:
    mol = "mol"


_amount_unit_scales = {"mol": si.mol, "1": 1}

# TODO: FIXME: Deal with different codata versions?
__codata_version__ = si.__codata_version__

standard_shnitsel_units = {
    unit_dimensions.length: length.Bohr,
    unit_dimensions.energy: energy.Hartree,
    unit_dimensions.force: force.Hartree_per_Bohr,
    # TODO: FIXME: Check which default time unit is used in shnitsel
    unit_dimensions.time: time.femto_seconds,
    unit_dimensions.nacs: "1",
    unit_dimensions.dipole: dipole.Debye,
    # "dipole_trans": "1",
    unit_dimensions.soc: socs.au,
}

"""
Previously used settings for SHARC input:
    attrs = {
        'atXYZ': {'long_name': "positions", 'units': 'Bohr', 'unitdim': 'Length'},
        'energy': {'units': 'hartree', 'unitdim': 'Energy'},
        'e_kin': {'units': 'hartree', 'unitdim': 'Energy'},
        'dip_perm': {'long_name': "permanent dipoles", 'units': 'au'},
        'dip_trans': {'long_name': "transition dipoles", 'units': 'au'},
        'sdiag': {'long_name': 'active state (diag)'},
        'astate': {'long_name': 'active state (MCH)'},
        'forces': {'units': 'hartree/bohr', 'unitdim': 'Force'},
        'nacs': {'long_name': "nonadiabatic couplings", 'units': 'au'},
    }
"""


def get_default_input_attributes(
    kind: Literal["sharc", "newtonx", "ase", "pyrai2md"],
    loading_parameters: LoadingParameters = None,
) -> Dict[str, Dict[str, str]]:
    format_default_units = standard_units_of_formats[kind]

    def override_defaults(unit_dimension, variable_name):
        if (
            loading_parameters is not None
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
            "long_name": "Transition dipoles",
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
        "state_names": {"long_name": "String representations of the states."},
        "state_types": {
            "long_name": "Index to indicate whether the state is singlet (1), doublet (2), or triplet(3)"
        },
        "state2": {"long_name": "The second state to build state combinations out of"},
        "statecomb": {
            "long_name": "Combination of two states used to index inter-state properties"
        },
        "atNames": {"long_name": "Names of atomic elements (short form)"},
        "atNums": {"long_name": "Periodic number of atomic elements"},
        "forces": {
            "long_name": "Per-atom forces",
            "unitdim": unit_dimensions.force,
            "units": format_default_units[unit_dimensions.force],
        },
        "nacs": {
            "long_name": "nonadiabatic couplings",
            "unitdim": unit_dimensions.nacs,
            "units": format_default_units[unit_dimensions.nacs],
        },
    }

    return res


standard_units_of_formats = {
    "sharc": {
        unit_dimensions.length: length.Bohr,
        unit_dimensions.energy: energy.Hartree,
        unit_dimensions.force: force.Hartree_per_Bohr,
        unit_dimensions.time: time.femto_seconds,  # TODO: Confirm default time unit
        unit_dimensions.nacs: nacs.au,
        unit_dimensions.dipole: dipole.au,
        # "dipole_trans": dipole.au,
        unit_dimensions.soc: socs.au,
    },
    "ase": {
        unit_dimensions.length: length.Bohr,
        unit_dimensions.energy: energy.Hartree,
        unit_dimensions.force: force.Hartree_per_Bohr,
        # Ang/sqrt(u/eV) approx 10fs
        unit_dimensions.time: time.ase_time_units,
        unit_dimensions.nacs: nacs.au,
        unit_dimensions.dipole: dipole.au,
        # "dipole_trans": dipole.au,
        unit_dimensions.soc: socs.au,
    },
    "newtonx": {  # Generally uses atomic units
        unit_dimensions.length: length.Angstrom,  # TODO: FIXME: Until 1.3 it was AU
        unit_dimensions.energy: energy.eV,  # Hartree or eV, it depends || energy.Hartree,
        unit_dimensions.force: force.Hartree_per_Bohr,
        unit_dimensions.time: time.femto_seconds,
        unit_dimensions.nacs: nacs.au,
        unit_dimensions.dipole: dipole.Debye,
        # "dipole_trans": dipole.au,
        unit_dimensions.soc: socs.au,
    },
    "shnitsel": standard_shnitsel_units,
    "xyz": {
        unit_dimensions.length: length.Bohr,
    },
    # Below is tentative support for Pyrai2MD file reading
    "pyrai2md": {  # TODO: FIXME: Pyrai2MD parameters are quite uncertain based on their documentation
        unit_dimensions.length: length.Angstrom,
        unit_dimensions.energy: energy.Hartree,
        unit_dimensions.force: force.Hartree_per_Bohr,
        # Ang/sqrt(u/eV) approx 10fs
        unit_dimensions.time: time.au,
        unit_dimensions.nacs: nacs.au,
        unit_dimensions.dipole: dipole.au,
        # "dipole_trans": dipole.au,
        unit_dimensions.soc: socs.au,
    },
}
