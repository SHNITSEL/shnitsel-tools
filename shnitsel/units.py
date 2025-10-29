
import ase.units as si
import numpy as np

# TODO: FIXME: Check all atomic units for correctness


class time:
    pico_seconds = 'ps'
    femto_seconds = 'fs'
    nano_seconds = 'ns'
    seconds = 's'
    ase_time_units = 'ase_time'


time_unit_scales = {
    time.pico_seconds: si.fs*1e3,
    time.femto_seconds: si.fs,
    time.nano_seconds: si.fs*1e6,
    time.seconds: si.second,
    'seconds': si.second*1e-12,
    # 'au': si._aut,
    # ASE uses this time unit: Ang/sqrt(u/eV) which might differ slightly from 10fs
    time.ase_time_units: si.Angstrom/np.sqrt(si._amu/si.eV)
}


class energy:
    Hartree = 'Hartree'
    eV = 'eV'
    keV = 'keV'
    J = 'J'
    kJ = 'kJ'
    kcal = 'kcal'


_energy_unit_scales = {
    energy.Hartree: si.Hartree,
    energy.eV: si.eV,
    energy.keV: si.eV*1e3,
    'joule': si.J,
    energy.J: si.J,
    energy.kJ: si.kJ,
    energy.kcal: si.kcal
}


class force_units:
    Hartree_per_Bohr = 'Hartree/Bohr'
    Hartree_per_Angstrom = 'Hartree/Angstrom'
    eV_per_Bohr = 'eV/Bohr'
    eV_per_Angstrom = 'eV/Angstrom'
    Newton = 'N'


_force_unit_scales = {
    force_units.Hartree_per_Bohr: si.Hartree/si.Bohr,
    force_units.eV_per_Bohr: si.eV/si.Bohr,
    force_units.Hartree_per_Angstrom: si.Hartree/si.Angstrom,
    force_units.eV_per_Angstrom: si.eV/si.Angstrom,
    'Hartree/A': si.Hartree/si.Angstrom,
    'eV/A': si.eV/si.Angstrom,
    force_units.Newton: si.J/si.m,
    'au': si.Hartree/si.Bohr,
}


class distance_units:
    Bohr = 'Bohr'
    Angstrom = 'Angstrom'
    meter = 'meter'
    pico_meter = 'pm'
    nano_meter = 'nm'


_distance_unit_scales = {
    distance_units.Bohr: si.Bohr,
    distance_units.Angstrom: si.Angstrom,
    'A': si.Angstrom,
    distance_units.meter: si.m,
    distance_units.pico_meter: si.nm*1e-3,
    distance_units.nano_meter: si.nm,
    "1": 1.0,
    "au": si.Bohr,
}

# An alias
length_units = distance_units
_length_unit_scales = _distance_unit_scales


class dipole_units:
    Debye = 'Debye'


_dipole_unit_scales = {
    "1": 1.0,
    dipole_units.Debye: si.Debye,  # 1/0.3934303,
}


class amount_units:
    mol = 'mol'


_amount_unit_scales = {
    "mol": si.mol,
    "1": 1
}

# TODO: FIXME: Deal with different codata versions?
__codata_version__ = si.__codata_version__

standard_shnitsel_units = {
    "length": "Bohr",
    "energy": "Hartree",
    "force": "Hartree/Bohr",
    "time": "ps",  # TODO: FIXME: Check which default time unit is used in shnitsel
    "nacs": "1",
    "dipole": "debye",
    "dipole_trans": "1",
    "soc": "1"
}

standard_units_of_formats = {
    "sharc": {
        "length": "Bohr",
        "energy": "eV",
        "force": "eV/Angstrom",
        "time": "fs",  # TODO: FIXME: Check which default unit
        "nacs": "1",
        "dipole": "debye",
        "dipole_trans": "1",
        "soc": "1"
    },
    "xyz": {
        "length": "Bohr",
    },
    "ase": {
        "length": "Bohr",
        "energy": "Hartree",
        "force": "Hartree/Bohr",
        "time": "ase_time",  # Ang/sqrt(u/eV)
        "nacs": "1",
        "dipole": "debye",
        "dipole_trans": "1",
        "soc": "1"
    },
    "newtonx": {  # Generallz uses atomic units
        "length": "Angstrom",  # TODO: FIXME: Until 1.3 it was AU
        "energy": "eV",  # Hartree or eV, it depends
        "force": "Hartree/Bohr",
        "time": "fs",
        "nacs": "1",
        "dipole": "debye",
        "dipole_trans": "1",
        "soc": "1"
    },
    "shnitsel": standard_shnitsel_units
}
