

import numpy as np
import xarray as xr
from itertools import combinations
import pandas as pd
import logging
import os
import re
import math

from ..xyz import parse_xyz


def parse_newtonx(traj_path) -> xr.Dataset:
    with open(os.path.join(traj_path, 'RESULTS', 'nx.log')) as f:
        single_traj = parse_nx_log(f)

    # nsteps = single_traj.sizes['ts']

    with open(os.path.join(traj_path, 'RESULTS', 'dyn.xyz')) as f:
        atNames, atNums, atXYZ = parse_xyz(f)

    # single_traj.attrs['atNames'] = atNames
    dims = ['ts', 'atom', 'direction']
    single_traj = single_traj.assign_coords(
        {'atNames': ('atom', atNames), 'atNums': ('atom', atNums)})

    if (
        not single_traj.attrs['completed']
        and atXYZ.shape[0] == single_traj.sizes['ts'] + 1
    ):
        logging.info("Geometry file contains ts after error. Truncating.")
        atXYZ = atXYZ[:-1]

    single_traj['atXYZ'] = xr.DataArray(
        atXYZ, coords={k: single_traj.coords[k] for k in dims}, dims=dims
    )

    return single_traj


def parse_nx_log(f) -> xr.Dataset:
    completed = True

    # hack to deal with restarts obscuring real tmax
    real_tmax = float(0)
    for line in f:
        stripline = line.strip()
        if stripline.startswith('FINISHING'):
            real_tmax = max(real_tmax, float(stripline.split()[4]))
        elif stripline.startswith('xx:'):
            splitline = stripline.split()
            if splitline[1] == '::ERROR::':
                real_tmax = max(real_tmax, float(splitline[7]))
                completed = False

    logging.debug(f'found real_tmax: {real_tmax}')
    f.seek(0)

    # skip to settings
    for line in f:
        if line.strip().startswith('Initial parameters:'):
            break

    settings = {}
    for line in f:
        # blank line marks end of settings
        if line.strip() == '':
            break

        key, val = re.split(' *= *', line.strip())
        if '.' in val:
            val = float(val)
        else:
            val = int(val)
        settings[key] = val

    delta_t = settings['dt']
    max_ts = int(real_tmax / delta_t)
    nsteps = max_ts + 1
    nstates = settings['nstat']
    natoms = settings['Nat']
    assert isinstance(nstates, int)
    assert isinstance(natoms, int)

    # TODO: so far, we haven't needed to parse data that
    #       had dipoles etc.
    energy = np.full((nsteps, nstates), np.nan)
    # dip_all = np.full((nsteps, nstates, nstates, 3), np.nan)
    # phases = np.full((nsteps, nstates), np.nan)
    # sdiag = np.full((nsteps), -1, dtype=int)
    astate = np.full((nsteps), -1, dtype=int)
    forces = np.full((nsteps, natoms, 3), np.nan)
    nacs = np.full((nsteps, math.comb(nstates, 2), natoms, 3), np.nan)

    ts: int = 0
    time: float = 0
    # parse actual data
    for line in f:
        stripline = line.strip()

        # TODO: unfortunately, Newton-X reports current step number, time and
        #       active state at the end of a timestep.
        #       Currently we use this to set the step number and time for the
        #       following time step. This is confusing and possibly vulnerable to
        #       strangely-formatted data -- is there a guarantee that time steps are
        #       always in order? There's almost certainly a better way to do this.
        #       For example, instead of assuming times follow expected order,
        #       lookahead to next FINISHING
        if stripline.startswith('FINISHING'):
            # Set active state for _current_ time step
            astate[ts] = stripline.split()[8]

            # Set step number and time for _the following_ time step
            time = float(stripline.split()[4]) + delta_t
            ts = int(time / delta_t)
            logging.debug(f'finished ts {ts - 1}')

        elif stripline.startswith('Gradient vectors'):
            for iatom in range(natoms):
                forces[ts, iatom] = [float(n) for n in next(f).strip().split()]

        elif stripline.startswith('Nonadiabatic coupling vectors'):
            for icomb in range(math.comb(nstates, 2)):
                nacs[ts, icomb] = [float(n) for n in next(f).strip().split()]

        elif stripline.startswith('Energy ='):
            for istate in range(nstates):
                energy[ts, istate] = float(next(f).strip())

    states = list(range(1, nstates + 1))

    statecomb = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(
            states, 2), names=['from', 'to']),
        dim='statecomb',
    )

    coords = statecomb.merge(
        {
            'ts': np.arange(nsteps),
            'state': states,
            'state2': states,
            'atom': np.arange(natoms),
            # 'statecomb': np.arange(math.comb(nstates, 2)),
            'direction': ['x', 'y', 'z'],
        }
    )

    # TODO: Are these units even correct?
    return xr.Dataset(
        {
            'energy': (
                ['ts', 'state'],
                energy,
                {'units': 'hartree', 'unitdim': 'Energy'},
            ),
            # 'dip_all': (['ts', 'state', 'state2', 'direction'], dip_all),
            # 'dip_perm': (['ts', 'state', 'direction'], dip_perm),
            # 'dip_trans': (['ts', 'statecomb', 'direction'], dip_trans),
            # 'sdiag': (['ts'], sdiag),
            'astate': (['ts'], astate, {'long_name': 'active state'}),
            'forces': (
                ['ts', 'atom', 'direction'],
                forces,
                {'units': 'hartree/bohr', 'unitdim': 'Force'},
            ),
            # 'has_forces': (['ts'], has_forces),
            # 'phases': (['ts', 'state'], phases),
            'nacs': (
                ['ts', 'statecomb', 'atom', 'direction'],
                nacs,
                {'long_name': "nonadiabatic couplings", 'units': "au"},
            ),
        },
        coords=coords,
        attrs={
            'max_ts': max_ts,
            'real_tmax': real_tmax,
            'delta_t': delta_t,
            'completed': completed,
        },
    )
