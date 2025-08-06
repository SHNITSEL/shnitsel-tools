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
    # Read some settings from input
    # In particular, if a trajectory is extended by increasing
    # tmax and resuming, the header of output.dat will give
    # only the original nsteps, leading to an ndarray IndexError
    input_path = os.path.join(traj_path, 'input')
    if os.path.isfile(input_path):
        with open(input_path) as f:
            settings = parse_input(f)
        delta_t = float(settings['stepsize'])
        tmax = float(settings['tmax'])
        nsteps = int(tmax/delta_t) + 1
    else:
        delta_t = None
        tmax= None
        nsteps = None

    with open(os.path.join(traj_path, 'output.dat')) as f:
        single_traj = parse_trajout_dat(f, nsteps=nsteps)
    if nsteps is None:
        nsteps = single_traj.sizes['ts']

    with open(os.path.join(traj_path, 'output.xyz')) as f:
        atNames, atXYZ = parse_trajout_xyz(nsteps, f)

    single_traj.coords['atNames'] = 'atom', atNames

    single_traj['atXYZ'] = xr.DataArray(
        atXYZ,
        coords={k: single_traj.coords[k] for k in ['ts', 'atom', 'direction']},
        dims=['ts', 'atom', 'direction'],
    )
    if delta_t is not None:
        single_traj.attrs['delta_t'] = delta_t

    return single_traj


def parse_trajout_dat(f, nsteps: int | None = None):
    settings = {}
    for line in f:
        if line.startswith('*'):
            break

        parsed = line.strip().split()
        if len(parsed) == 2:
            settings[parsed[0]] = parsed[1]
        elif len(parsed) > 2:
            settings[parsed[0]] = parsed[1:]
        else:
            logging.warning("Key without value in settings of output.dat")

    nsteps_output_dat = int(settings['nsteps']) + 1  # let's not forget ts=0
    if nsteps is None or nsteps < nsteps_output_dat:
        nsteps = nsteps_output_dat
        logging.debug(f"nsteps = {nsteps}")
    else:
        logging.debug(f"(From input file) nsteps = {nsteps}")
    natoms = int(settings['natom'])  # yes, really 'natom', not 'natoms'!
    logging.debug(f"natoms = {natoms}")
    ezero = float(settings['ezero'])
    logging.debug(f"ezero = {ezero}")
    state_settings = [int(s) for s in settings['nstates_m']]
    state_settings += [0] * (3 - len(state_settings))
    nsinglets, ndoublets, ntriplets = state_settings
    nstates = nsinglets + 2 * ndoublets + 3 * ntriplets
    logging.debug(f"nstates = {nstates}")

    idx_table_nacs = {
        (si, sj): idx
        for idx, (si, sj) in enumerate(combinations(range(1, nstates+1), 2))
    }

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
            ts = int(next(f).strip())
            if ts != 0:
                logging.warning("Initial timestep's index is not 0")
            max_ts = max(max_ts, ts)
            break

    for index, line in enumerate(f):
        if line[0] != '!':
            continue

        if line.startswith('! 0 Step'):
            # update `ts` to current timestep #
            new_ts = int(next(f).strip())
            if new_ts != (ts or 0) + 1:
                logging.warning(f"Non-consecutive timesteps: {ts} -> {new_ts}")
            ts = new_ts
            max_ts = max(max_ts, ts)
            logging.debug(f"timestep = {ts}")

        if line.startswith('! 1 Hamiltonian'):
            for istate in range(nstates):
                energy[ts, istate] = float(next(f).strip().split()[istate * 2]) + ezero

        if line.startswith('! 3 Dipole moments'):
            direction = {'X': 0, 'Y': 1, 'Z': 2}[line.strip().split()[4]]
            for istate in range(nstates):
                linecont = next(f).strip().split()
                # delete every second element in list (imaginary values, all zero)
                dip_all[ts, istate, :, direction] = [float(i) for i in linecont[::2]]

        if line.startswith('! 4 Overlap matrix'):
            found_overlap = False
            phasevector = np.ones((nstates))

            wvoverlap = np.zeros((nstates, nstates))
            for j in range(nstates):
                linecont = next(f).strip().split()
                # delete every second element in list (imaginary values, all zero)
                wvoverlap[j] = [float(n) for n in linecont[::2]]

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
            pair = next(f).strip().split()
            sdiag[ts] = int(pair[0])
            astate[ts] = int(pair[1])

        if line.startswith('! 15 Gradients (MCH)'):
            state = int(line.strip().split()[-1]) - 1

            for atom in range(natoms):
                forces[ts, state, atom] = [float(n) for n in next(f).strip().split()]

        if line.startswith('! 16 NACdr matrix element'):
            linecont = line.strip().split()
            si, sj = int(linecont[-2]), int(linecont[-1])

            if si < sj: # elements (si, si) are all zero; elements (sj, si) = -(si, sj)
                sc = idx_table_nacs[(si, sj)]  # statecomb index
                for atom in range(natoms):
                    nacs[ts, sc, atom, :] = [float(n) for n in next(f).strip().split()]
            else:  # we can skip the block
                for _ in range(natoms):
                    next(f)


    # post-processing
    # np.diagonal swaps state and direction, so we transpose them back
    dip_perm = np.diagonal(dip_all, axis1=1, axis2=2).transpose(0, 2, 1)
    idxs_dip_trans = (slice(None), *np.triu_indices(nstates, k=1))
    dip_trans = dip_all[idxs_dip_trans].transpose(0, 2, 1)
    # has_forces = forces.any(axis=(1, 2, 3))

    if not max_ts + 1 <= nsteps:
        raise ValueError(
            f"Metadata declared {nsteps=} timesteps, but the "
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
            # 'state2': states,
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