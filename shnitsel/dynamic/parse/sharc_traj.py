import numpy as np
import xarray as xr
from itertools import combinations
import pandas as pd
import logging
import os
import re
import math

from .common import get_dipoles_per_xyz, dip_sep, get_triangular

def read_traj(traj_path):
    with open(os.path.join(traj_path, 'output.dat')) as f:
        single_traj = parse_trajout_dat(f)

    nsteps = single_traj.sizes['ts']

    with open(os.path.join(traj_path, 'output.xyz')) as f:
        atNames, atXYZ = parse_trajout_xyz(nsteps, f)

    with open(os.path.join(traj_path, 'input')) as f:
        settings = parse_input(f)

    single_traj.coords['atNames'] = 'atom', atNames

    single_traj['atXYZ'] = xr.DataArray(
        atXYZ,
        coords={k: single_traj.coords[k] for k in ['ts', 'atom', 'direction']},
        dims=['ts', 'atom', 'direction'],
    )
    single_traj.attrs['delta_t'] = float(settings['stepsize'])

    return single_traj


def parse_trajout_dat(f):
    settings = {}
    for line in f:
        if line.startswith('*'):
            break

        parsed = re.split(' +', line.strip())
        if len(parsed) == 2:
            settings[parsed[0]] = parsed[1]
        elif len(parsed) > 2:
            settings[parsed[0]] = parsed[1:]
        else:
            logging.warning("Key without value in settings of output.dat")

    nsteps = int(settings['nsteps']) + 1  # let's not forget ts=0
    logging.debug(f"nsteps = {nsteps}")
    natoms = int(settings['natom'])  # yes, really 'natom', not 'natoms'!
    logging.debug(f"natoms = {natoms}")
    ezero = float(settings['ezero'])
    logging.debug(f"ezero = {ezero}")
    state_settings = [int(s) for s in settings['nstates_m']]
    state_settings += [0] * (3 - len(state_settings))
    nsinglets, ndoublets, ntriplets = state_settings
    nstates = nsinglets + 2 * ndoublets + 3 * ntriplets
    logging.debug(f"nstates = {nstates}")

    # now we know the number of steps, we can initialize the data arrays:
    energy = np.full((nsteps, nstates), np.nan)
    e_kin = np.full((nsteps,), np.nan)
    dip_all = np.full((nsteps, nstates, nstates, 3), np.nan)
    phases = np.full((nsteps, nstates), np.nan)
    sdiag = np.full((nsteps), -1, dtype=int)
    astate = np.full((nsteps), -1, dtype=int)
    forces = np.full((nsteps, nstates, natoms, 3), np.nan)
    nacs = np.full((nsteps, math.comb(nstates, 2), natoms, 3), np.nan)

    max_ts = -1

    # skip through until initial step:
    for line in f:
        if line.startswith('! 0 Step'):
            ts = int(re.split(' +', next(f).strip())[-1])
            if ts != 0:
                logging.warning("Initial timestep's index is not 0")
            max_ts = max(max_ts, ts)
            break

    for index, line in enumerate(f):
        if line.startswith('! 0 Step'):
            # update `ts` to current timestep #
            new_ts = int(next(f).strip().split()[-1])
            if new_ts != (ts or 0) + 1:
                logging.warning(f"Non-consecutive timesteps: {ts} -> {new_ts}")
            ts = new_ts
            max_ts = max(max_ts, ts)
            logging.debug(f"timestep = {ts}")

        if line.startswith('! 1 Hamiltonian'):
            for istate in range(nstates):
                energy[ts, istate] = float(next(f).strip().split()[istate * 2]) + ezero

        if line.startswith('! 3 Dipole moments X'):
            x_dip = get_dipoles_per_xyz(file=f, n=nstates, m=nstates)
            dip_all[ts, :, :, 0] = x_dip

        if line.startswith('! 3 Dipole moments Y'):
            y_dip = get_dipoles_per_xyz(file=f, n=nstates, m=nstates)
            dip_all[ts, :, :, 1] = y_dip

        if line.startswith('! 3 Dipole moments Z'):
            z_dip = get_dipoles_per_xyz(file=f, n=nstates, m=nstates)
            dip_all[ts, :, :, 2] = z_dip

        if line.startswith('! 4 Overlap matrix'):
            found_overlap = False
            phasevector = np.ones((nstates))

            wvoverlap = np.zeros((nstates, nstates))
            for j in range(nstates):
                linecont = [float(n) for n in re.split(' +', next(f).strip())]
                # delete every second element in list (imaginary values, all zero)
                wvoverlap[j] = linecont[::2]

            for istate in range(nstates):
                if np.abs(wvoverlap[istate, istate]) >= 0.5:
                    found_overlap = True
                    if wvoverlap[istate, istate] >= 0.5:
                        phasevector[istate] = +1
                    else:
                        phasevector[istate] = -1

            if found_overlap:
                phases[ts] = phasevector

        if line.startswith('! 7 Ekin'):
            e_kin[ts] = float(next(f).strip())

        if line.startswith('! 8 states (diag, MCH)'):
            pair = re.split(' +', next(f).strip())
            sdiag[ts] = int(pair[0])
            astate[ts] = int(pair[1])

        if line.startswith('! 15 Gradients (MCH)'):
            state = int(re.split(' +', line.strip())[-1]) - 1

            for atom in range(natoms):
                forces[ts, state, atom] = [
                    float(n) for n in re.split(' +', next(f).strip())
                ]

        if line.startswith('! 16 NACdr matrix element'):
            linecont = re.split(' +', line.strip())
            si, sj = int(linecont[-2]) - 1, int(linecont[-1]) - 1

            if si == sj == 0:
                nacs_matrix = np.zeros((nstates, nstates, natoms, 3))

            for atom in range(natoms):
                nacs_matrix[si, sj, atom] = [
                    float(n) for n in re.split(' +', next(f).strip())
                ]

            # get upper triangular of nacs matrix
            nacs_tril = get_triangular(nacs_matrix)
            nacs[ts] = nacs_tril

    # post-processing
    dip_perm = np.full((nsteps, nstates, 3), np.nan)
    dip_trans = np.full((nsteps, math.comb(nstates, 2), 3), np.nan)
    has_forces = np.zeros((nsteps), dtype=bool)

    for ts in range(nsteps):
        p, t = dip_sep(dip_all[ts])
        dip_perm[ts] = p
        dip_trans[ts] = t

        if np.any(forces[ts]):
            has_forces[ts] = True

    if not max_ts + 1 <= nsteps:
        raise ValueError(
            f"The output.dat header declared {nsteps=} timesteps, but the "
            f"greatest timestep index was {max_ts + 1=}"
        )
    completed = max_ts + 1 == nsteps

    # Currently 1-based numbering corresponding to internal SHARC usage.
    # Ultimately aiming to replace numbers with labels ('S0', 'S1', ...),
    # but that has disadvantages in postprocessing.
    states = np.arange(1, nstates + 1)

    statecomb = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(states, 2), names=['from', 'to']),
        dim='statecomb',
    )

    coords = statecomb.merge(
        {
            'ts': np.arange(nsteps),
            'state': states,
            'state2': states,
            'atom': np.arange(natoms),
            'direction': ['x', 'y', 'z'],
        }
    )

    res = xr.Dataset(
        {
            'energy': (
                ['ts', 'state'],
                energy,
                {'units': 'hartree', 'unitdim': 'Energy'},
            ),
            'e_kin': (
                ['ts'],
                e_kin,
                {'units': 'hartree', 'unitdim': 'Energy'},
            ),
            # 'dip_all': (['ts', 'state', 'state2', 'direction'], dip_all),
            'dip_perm': (
                ['ts', 'state', 'direction'],
                dip_perm,
                {'long_name': "permanent dipoles", 'units': 'au'},
            ),
            'dip_trans': (
                ['ts', 'statecomb', 'direction'],
                dip_trans,
                {'long_name': "transition dipoles", 'units': 'au'},
            ),
            'sdiag': (['ts'], sdiag, {'long_name': 'active state (diag)'}),
            'astate': (['ts'], astate, {'long_name': 'active state (MCH)'}),
            'forces': (
                ['ts', 'state', 'atom', 'direction'],
                forces,
                {'units': 'hartree/bohr', 'unitdim': 'Force'},
            ),
            # 'has_forces': (['ts'], has_forces),
            'phases': (['ts', 'state'], phases),
            'nacs': (
                ['ts', 'statecomb', 'atom', 'direction'],
                nacs,
                {'long_name': "nonadiabatic couplings", 'units': 'au'},
            ),
        },
        coords=coords,
        attrs={'max_ts': max_ts, 'completed': completed},
    )

    if not completed:
        res = res.sel(ts=res.ts <= res.attrs['max_ts'])

    return res


def parse_trajout_xyz(nsteps, f):
    first = next(f)
    assert first.startswith(' ' * 6)
    natoms = int(first.strip())

    atNames = np.full((natoms), '')
    atXYZ = np.full((nsteps, natoms, 3), np.nan)

    ts = 0

    for index, line in enumerate(f):
        if 't=' in line:
            assert ts < nsteps, f"ts={ts}, nsteps={nsteps}"
            for atom in range(natoms):
                linecont = re.split(' +', next(f).strip())
                if ts == 0:
                    atNames[atom] = linecont[0]
                atXYZ[ts, atom] = [float(n) for n in linecont[1:]]
            ts += 1

    return (atNames, atXYZ)

def parse_input(f):
    settings = {}
    for line in f:
        parsed = line.strip().split()
        if len(parsed) == 2:
            settings[parsed[0]] = parsed[1]
        elif len(parsed) > 2:
            settings[parsed[0]] = parsed[1:]
        elif len(parsed) == 1:
            settings[parsed[0]] = True
    return settings