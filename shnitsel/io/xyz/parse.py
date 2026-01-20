from io import TextIOWrapper
from typing import List, Tuple
import numpy as np

from shnitsel.core._api_info import internal
from shnitsel.data.atom_helpers import get_atom_number_from_symbol


@internal()
def parse_xyz(f: TextIOWrapper) -> tuple[list[str], list[int], np.ndarray]:
    """Read the inputs from a text file stream into a tuple of atom names, atom numbers and positions.

    Parameters
    ----------
    f : TextIOWrapper
        The File wrapper providing the contents of a `.xyz` file

    Returns
    -------
    tuple[list[str], list[int], np.ndarray]
        The tuple of (atom_symbols, atom_numbers, atom_positions), where the latter has dimensions [timestep][atom][direction]
    """

    # TODO: f has to be an open file, I assume?
    natoms = int(next(f).strip())
    # ts = 0

    atXYZ = []  # np.full((nsteps, natoms, 3), np.nan)
    atNames = []  # np.full((natoms), '')
    atNums = []  # np.full((natoms), '')

    # Skip one line of the input file
    next(f)
    thisXYZ = np.full((natoms, 3), np.nan)
    for iatom in range(natoms):
        geometry_line = next(f).strip().split()
        # atNames[iatom] = geometry_line[0]
        # atXYZ[ts, iatom] = [float(n) for n in geometry_line[1:]]
        atNames.append(geometry_line[0])
        atNums.append(get_atom_number_from_symbol(geometry_line[0]))
        thisXYZ[iatom] = [float(n) for n in geometry_line[1:]]
    atXYZ.append(thisXYZ)

    for line in f:
        assert line.startswith(' '), f'Expected empty line but got content: {line!r}'
        # ts += 1
        line = next(f)
        assert line.startswith(' '), f'Expected empty line but got content: {line!r}'

        thisXYZ = np.full((natoms, 3), np.nan)
        for iatom, atName in enumerate(atNames):
            geometry_line = next(f).strip().split()
            assert geometry_line[0] == atName, "Inconsistent atom order"
            # atXYZ[ts, iatom] = [float(n) for n in geometry_line[1:]]
            thisXYZ[iatom] = [float(n) for n in geometry_line[1:]]
        atXYZ.append(thisXYZ)

    return (atNames, atNums, np.stack(atXYZ, axis=0))


@internal()
def get_dipoles_per_xyz(file: TextIOWrapper, n: int, m: int) -> np.ndarray:
    """Read full dipole matrix from an xyz file

    Parameters
    ----------
    file : TextIOWrapper
        Wrapper providing the contents of an `.xyz` file
    n : int
        First index length of the dipole matrix
    m : int
        Second index length of the dipole matrix

    Returns
    -------
    np.ndarray
        The matrix of the dipole contents
    """
    dip = np.zeros((n, m))
    for istate in range(n):
        linecont = next(file).strip().split()
        # delete every second element in list (imaginary values, all zero)
        dip[istate] = [float(i) for i in linecont[::2]]

    return dip
