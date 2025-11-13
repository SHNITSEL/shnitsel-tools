import ase.units as si
import numpy as np

# TODO: FIXME: Check all atomic units for correctness


class unit_dimensions:
    length = "length"
    energy = "energy"
    force = "force"
    dipole = "dipole"
    nacs = "nacs"
    time = "time"
    socs = "socs"
    charge = "charge"


class charge:
    e = "e"
    Coulomb = "C"


_charge_unit_scales = {
    charge.e: si._e,
    charge.Coulomb: si.C,
}


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
    "seconds": si.second,
    time.au: si.AUT,
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
    CoulombMeter = "Cm"
    Debye = "debye"
    au = "au"


_dipole_unit_scales = {
    dipole.CoulombMeter: si.C * si.m,
    dipole.Debye: si.Debye,  # si.C* si.m * 1.e-21 / si._c,
    # TODO: FIXME: set actual atomic dipole unit.
    dipole.au: si.C
    * si.m
    * si._e
    * (1.0e-10 * si.Bohr),  # TODO: Check: 1D = 0.393456a.u.
}


class nacs:
    # TODO: FIXME: Figure out nacs units
    one_per_Bohr = "one_per_Bohr"
    au = "au"


_nacs_unit_scales = {
    nacs.one_per_Bohr: 1.0,
    nacs.au: 1.0,
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
    unit_dimensions.nacs: "1",  # TODO: FIXME: NACS in molcas: 1/Bohr, SHARC liest 1/Bohr
    unit_dimensions.dipole: dipole.Debye,
    # "dipole_trans": "1",
    unit_dimensions.socs: socs.au,
    unit_dimensions.charge: charge.e,
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


standard_units_of_formats = {
    "sharc": {  # These units are wild and only apply to certain files...
        unit_dimensions.length: length.Bohr,
        # Units are in eV relative to the ezero setting.
        unit_dimensions.energy: energy.eV,
        # TODO: FIXME: The SHARC documentation does not state the unit, the output.dat file does not state the unit, but this seems most likely. output.lis lists eV/Ang, which would not match the other data in output.log at all.
        unit_dimensions.force: force.eV_per_Bohr,
        unit_dimensions.time: time.femto_seconds,
        unit_dimensions.nacs: nacs.one_per_Bohr,
        unit_dimensions.dipole: dipole.au,
        # "dipole_trans": dipole.au,
        unit_dimensions.socs: socs.au,
        unit_dimensions.charge: charge.e,
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
        unit_dimensions.socs: socs.au,
        unit_dimensions.charge: charge.e,
    },
    "newtonx": {  # Generally uses atomic units?
        unit_dimensions.length: length.Angstrom,  # TODO: FIXME: Until 1.3 it was AU
        unit_dimensions.energy: energy.eV,  # Hartree or eV, it depends || energy.Hartree,
        unit_dimensions.force: force.Hartree_per_Bohr,
        unit_dimensions.time: time.femto_seconds,
        unit_dimensions.nacs: nacs.au,
        unit_dimensions.dipole: dipole.Debye,
        # "dipole_trans": dipole.au,
        unit_dimensions.socs: socs.au,
        unit_dimensions.charge: charge.e,
    },
    "shnitsel": standard_shnitsel_units,
    "xyz": {
        unit_dimensions.length: length.Angstrom,
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
        unit_dimensions.socs: socs.au,
        unit_dimensions.charge: charge.e,
    },
}
