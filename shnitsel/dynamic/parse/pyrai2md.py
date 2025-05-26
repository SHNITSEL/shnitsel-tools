import os
from itertools import combinations
from glob import glob
import xarray as xr
import pandas as pd
import numpy as np


def parse_md_energies(path):
    df = pd.read_csv(path, sep=r'\s+', header=None, skiprows=1).set_index(0)
    df.index.name = 'time'
    energy = df.loc[:, 4:]
    nstates = len(energy.columns)
    return (
        xr.Dataset.from_dataframe(energy)
        .to_array('state')
        .assign_coords(state=np.arange(1, nstates + 1))
    )


def parse_log(f, nsteps):
    ## First, read the settings.
    # *---------------------------------------------------*
    # |                                                   |
    # |          Nonadiabatic Molecular Dynamics          |
    # |                                                   |
    # *---------------------------------------------------*
    #
    #
    # State order:         1   2   3
    # Multiplicity:        1   1   1
    #
    # QMMM key:         None
    # QMMM xyz          Input
    # Active atoms:     45
    # Inactive atoms:   0
    # Link atoms:       0
    # Highlevel atoms:  45
    # Midlevel atoms:   0
    # Lowlevel atoms:   0

    while not next(f).startswith(' *---'):
        pass

    for _ in range(18):
        stripline = next(f).strip()
        if stripline.startswith('State order:'):
            states = np.asarray(stripline.split()[2:], dtype=int)
            nstates = len(states)
        if stripline.startswith('Multiplicity:'):
            multiplicity = stripline.split()[1:]
            for m in multiplicity:
                if m != '1':
                    raise ValueError("Only handling singlets for now.")
        if stripline.startswith('Active atoms:'):
            natoms = int(stripline.split()[2])
    del multiplicity, m

    # Set up numpy arrays
    astate = np.full((nsteps), -1, dtype=int)
    forces = np.full((nsteps, nstates, natoms, 3), np.nan)
    atXYZ = np.full((nsteps, natoms, 3), np.nan)
    atNames = np.full((natoms), '', dtype=str)
    got_atNames = False
    veloc = np.full((nsteps, natoms, 3), np.nan)
    dcmat = np.full((nsteps, nstates, nstates), np.nan)

    for line in f:
        ## The start of a timestep
        # Iter:        1  Ekin =           0.1291084223229551 au T =   300.00 K dt =         20 CI:   3
        # Root chosen for geometry opt   2
        if line.startswith('  Iter:'):
            ts = int(line.strip().split()[1]) - 1

            for _ in range(10):
                ## Get active state
                # A surface hopping is not allowed
                # **
                # At state:   2
                line = next(f)
                if line.startswith('  At state'):
                    astate[ts] = int(line.strip().split()[2])
                    break
                # A surface hopping event happened
                # **
                # From state:   2 to state:   3 *
                elif line.startswith('  From state'):
                    astate[ts] = int(line.strip().split()[5])
                    break
            else:
                raise ValueError(f"No state info found for Iter: {ts+1}")

        ## Positions:
        #   &coordinates in Angstrom
        # -------------------------------------------------------------------------------
        # C          0.5765950000000000     -0.8169010000000000     -0.0775610000000000
        # C          1.7325100000000000     -0.1032670000000000      0.1707480000000000
        # -------------------------------------------------------------------------------
        if line.startswith('  &coordinates'):
            assert next(f).startswith('---')
            if got_atNames:
                for iatom in range(natoms):
                    atXYZ[ts, iatom] = np.asarray(
                        next(f).strip().split()[1:], dtype=float
                    )
            else:
                for iatom in range(natoms):
                    content = next(f).strip().split()
                    atXYZ[ts, iatom] = np.asarray(content[1:], dtype=float)
                    atNames[iatom] = str(content[0])

            assert next(f).startswith('---')

        ## Velocities:
        #   &velocities in Bohr/au
        # -------------------------------------------------------------------------------
        # C          0.0003442000000000      0.0001534200000000     -0.0000597200000000
        # C         -0.0005580000000000      0.0003118300000000     -0.0000154900000000
        # -------------------------------------------------------------------------------
        if line.startswith('  &velocities'):
            assert next(f).startswith('---')
            for iatom in range(natoms):
                veloc[ts, iatom] = np.asarray(next(f).strip().split()[1:], dtype=float)
            assert next(f).startswith('---')

        ## Forces:
        #   &gradient state               1 in Eh/Bohr
        # -------------------------------------------------------------------------------
        # C         -0.0330978534152795      0.0073099255379017      0.0082666356536386
        # C          0.0313629524413876      0.0196036465968827      0.0060952442704520
        # -------------------------------------------------------------------------------
        if line.startswith('  &gradient'):
            istate = int(line.strip().split()[2]) - 1
            assert next(f).startswith('---')
            for iatom in range(natoms):
                forces[ts, istate, iatom] = np.asarray(
                    next(f).strip().split()[1:], dtype=float
                )
            assert next(f).startswith('---')

        ## Derivative coupling matrix:
        #  &derivative coupling matrix
        # -------------------------------------------------------------------------------
        #       0.0000000000000000       0.0000000000000004      -0.0000000000000001
        #      -0.0000000000000004       0.0000000000000000       0.0000000000000003
        #       0.0000000000000001      -0.0000000000000003       0.0000000000000000
        # -------------------------------------------------------------------------------
        if line.startswith('  &derivative coupling matrix'):
            assert next(f).startswith('---')
            for istate1 in range(nstates):
                dcmat[ts, istate1] = np.asarray(next(f).strip().split(), dtype=float)
            assert next(f).startswith('---')

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
            'atNames': ('atom', atNames),
            # 'statecomb': np.arange(math.comb(nstates, 2)),
            'direction': ['x', 'y', 'z'],
        }
    )

    return xr.Dataset(
        {
            # 'dip_all': (['ts', 'state', 'state2', 'direction'], dip_all),
            # 'dip_perm': (['ts', 'state', 'direction'], dip_perm),
            # 'dip_trans': (['ts', 'statecomb', 'direction'], dip_trans),
            # 'sdiag': (['ts'], sdiag),
            'astate': (['ts'], astate, {'long_name': 'active state'}),
            'forces': (
                ['ts', 'state', 'atom', 'direction'],
                forces,
                {'units': 'hartree/bohr', 'unitdim': 'Force'},
            ),
            # 'has_forces': (['ts'], has_forces),
            # 'phases': (['ts', 'state'], phases),
            # 'nacs': (
            #     ['ts', 'statecomb', 'atom', 'direction'],
            #     nacs,
            #     {'long_name': "nonadiabatic couplings", 'units': "au"},
            # ),
            'atXYZ': (['ts', 'atom', 'direction'], atXYZ),
            'dcmat': (['ts', 'state', 'state2'], dcmat),
        },
        coords=coords,
        attrs={
            # 'max_ts': max_ts,
            # 'real_tmax': real_tmax,
            # 'delta_t': delta_t,
            # 'completed': completed,
        },
    )


def read_traj(traj_path):
    md_energies_paths = glob(os.path.join(traj_path, '*.md.energies'))
    if (n := len(md_energies_paths)) != 1:
        raise FileNotFoundError(
            "Expected to find a single file ending with '.md.energies' "
            f"but found {n} files: {md_energies_paths}"
        )
    log_paths = glob(os.path.join(traj_path, '*.log'))
    if (n := len(md_energies_paths)) != 1:
        raise FileNotFoundError(
            "Expected to find a single file ending with '.log' "
            f"but found {n} files: {log_paths}"
        )

    energy = parse_md_energies(md_energies_paths[0])
    nsteps = energy.sizes['time']
    with open(os.path.join(traj_path, log_paths[0])) as f:
        single_traj = parse_log(f, nsteps)
    single_traj = single_traj.rename(ts='time').assign_coords(time=energy['time'])
    single_traj['energy'] = energy

    return single_traj