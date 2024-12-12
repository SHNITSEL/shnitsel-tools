from dboverview.parsecommon import *
import numpy as np
import xarray as xr
from itertools import combinations
import pandas as pd
import logging, os, re, math

##################################
# Functions for reading TRAJ files
# read_trajs parses everything in traj/
# read_traj to read a single TRAJ_0000n
##################################

def list_trajs(trajs_path=None):
    if not trajs_path:
        trajs_path = os.path.join(os.getcwd(), 'traj')

    names = sorted(
        [n
         for n in os.listdir(trajs_path)
         if n.startswith('TRAJ_0')
            and 'output.dat' in os.listdir(os.path.join(trajs_path, n))])
    
    return [
        (int(n[5:]), os.path.join(trajs_path, n))
        for n in names
    ]

# TODO: adapt to xarray!
# def read_trajs(nstates=3, natoms=12, trajs_path='', subset=None):  
#     data = {}
#     for number, traj_path in list_trajs(trajs_path):
#         if subset is None or number in subset:
#             print(traj_path)
#             data[number] = read_traj(nstates, natoms, traj_path)  

#     return Traj(data, nsinglets=nsinglets, ntriplets=ntriplets,
#                 natoms=natoms)

def read_traj(traj_path):
    with open(os.path.join(traj_path, 'output.dat')) as f:
        single_traj = parse_trajout_dat(f)
    
    nsteps = single_traj.sizes['ts']

    with open(os.path.join(traj_path, 'output.xyz')) as f:
        atNames, atXYZ = parse_trajout_xyz(nsteps, f)
    
    # single_traj.attrs['atNames'] = atNames
    single_traj = single_traj.assign_coords(
        {'atNames': ('atom', atNames)})
    single_traj['atXYZ'] = xr.DataArray(
        atXYZ,
        coords={k:single_traj.coords[k] for k in ['ts', 'atom', 'direction']},
        dims=['ts', 'atom', 'direction']
    )

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
            
    nsteps = int(settings['nsteps']) + 1 # let's not forget ts=0
    logging.debug(f"nsteps = {nsteps}")
    natoms = int(settings['natom']) # yes, really 'natom', not 'natoms'!
    logging.debug(f"natoms = {natoms}")
    ezero = float(settings['ezero'])
    logging.debug(f"ezero = {ezero}")
    state_settings = [int(s) for s in settings['nstates_m']]
    state_settings += [0] * (3 - len(state_settings))
    nsinglets, ndoublets, ntriplets = state_settings
    nstates = nsinglets + 2 * ndoublets + 3 * ntriplets
    print(nstates)
    logging.debug(f"nstates = {nstates}")

    # now we know the number of steps, we can initialize the data arrays:
    energies = np.full((nsteps, nstates), np.nan)
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
                energies[ts, istate] = \
                    float(next(f).strip().split()[istate * 2]) + ezero

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

            if found_overlap == True:
                phases[ts] = phasevector

        if line.startswith('! 8 states (diag, MCH)'):
            pair = re.split(' +', next(f).strip())
            sdiag[ts] = int(pair[0])
            astate[ts] = int(pair[1])

        if line.startswith('! 15 Gradients (MCH)'):
            # if astate[ts] == np.nan:
            #     raise ValueError(f"Gradients given before active state for timestep {ts}")

            state = int(re.split(' +', line.strip())[-1]) - 1

            # if gstate == astate[ts]:
            for atom in range(natoms):
                forces[ts, state, atom] = [float(n) for n in re.split(' +', next(f).strip())] 

        if line.startswith('! 16 NACdr matrix element'):
            linecont = re.split(' +', line.strip())
            si, sj = int(linecont[-2])-1, int(linecont[-1])-1

            if si == sj == 0:
                nacs_matrix = np.zeros((nstates, nstates, natoms, 3))

            for atom in range(natoms):
                nacs_matrix[si, sj, atom] = [float(n) for n
                                             in re.split(' +', next(f).strip())]

            # get upper triangular of nacs matrix
            nacs_tril = get_triangular(nacs_matrix)
            nacs[ts] = nacs_tril

    # post-processing
    dip_perm = np.full((nsteps, nstates, 3), np.nan)
    dip_trans = np.full((nsteps, math.comb(nstates, 2), 3), np.nan)
    has_forces = np.zeros((nsteps), dtype=bool)

    for ts in range(nsteps):
        a = dip_sep(dip_all[ts])
        p, t = a
        dip_perm[ts] = p
        dip_trans[ts] = t

        if np.any(forces[ts]):
            has_forces[ts] = True
    
    # Currently 1-based numbering corresponding to internal SHARC usage.
    # Ultimately aiming to replace numbers with labels ('S0', 'S1', ...),
    # but that has disadvantages in postprocessing.
    states = np.arange(1, nstates + 1)

    statecomb = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(
            combinations(states, 2),
            names=['from', 'to']
        ),
        dim='statecomb'
    )

    coords = statecomb.merge({
        'ts': np.arange(nsteps),
        'state': states,
        'state2': states,
        'atom': np.arange(natoms),
        # 'statecomb': np.arange(math.comb(nstates, 2)),
        'direction': ['x', 'y', 'z']
    })

    return xr.Dataset(
        {
            'energies': (['ts', 'state'], energies,
                         {'unit': 'hartree', 'unitdim': 'Energy'}),
            'dip_all': (['ts', 'state', 'state2', 'direction'], dip_all),
            'dip_perm': (['ts', 'state', 'direction'], dip_perm),
            'dip_trans': (['ts', 'statecomb', 'direction'], dip_trans),
            'sdiag': (['ts'], sdiag),
            'astate': (['ts'], astate),
            'forces': (['ts', 'state', 'atom', 'direction'], forces,
                       {'unit': 'hartree/bohr', 'unitdim': 'Force'}),
            'has_forces': (['ts'], has_forces),
            'phases': (['ts', 'state'], phases),
            'nacs': (['ts', 'statecomb', 'atom', 'direction'], nacs)
        },
        coords=coords,
        attrs={'max_ts': max_ts}
    )

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
                
            #atXYZ_align = align_mol(atXYZ=atXYZ, align_index=[0,1,2])
            #traj_data[number][timestep]['atXYZ_align'] = atXYZ_align
            ts += 1

    # atNames  = np.char.encode(atNames, 'utf8')

    return (atNames, atXYZ)
