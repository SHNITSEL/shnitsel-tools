
import os
from itertools import combinations
from glob import glob
import pathlib
import xarray as xr
import pandas as pd
import numpy as np

from shnitsel.io.helpers import LoadingParameters, PathOptionsType


def parse_pyrai2md(traj_path: PathOptionsType,
                   loading_parameters: LoadingParameters | None = None) -> xr.Dataset:
    """Function to read a trajector of the PyrAI2md format.

    Args:
        pathlist (PathOptionsType): Path to the directory containing a PyrAI2md output file list
        loading_parameters (LoadingParameters | None, optional): Parameter settings for e.g. standard units or state names.

    Returns:
        xr.Dataset: The Dataset object containing all of the loaded data in default shnitsel units
    """
    # TODO: FIXME: use loading_parameters to configure units and state names
    md_energies_paths = glob(os.path.join(traj_path, '*.md.energies'))
    if (n := len(md_energies_paths)) != 1:
        raise FileNotFoundError(
            "Expected to find a single file ending with '.md.energies' "
            f"but found {n} files: {md_energies_paths}"
        )
    log_paths = glob(os.path.join(traj_path, '*.log'))
    if (n := len(log_paths)) != 1:
        raise FileNotFoundError(
            "Expected to find a single file ending with '.log' "
            f"but found {n} files: {log_paths}"
        )

    energy = parse_md_energies(md_energies_paths[0])
    with open(os.path.join(log_paths[0])) as f:
        single_traj = parse_log(f)

    # TODO: FIXME: conflicting dimension sizes "time". We need to deal with trajectory not finishing its full run.
    # One test trajectory did not finish its full number of steps as denoted in the log, so we need to check if the trajectory has finished
    # before sizing the output.
    single_traj = single_traj.rename(
        ts='time').assign_coords(time=energy['time'])
    single_traj['energy'] = energy
    single_traj['energy'].attrs['units'] = 'hartree'

    return single_traj


def parse_md_energies(path) -> xr.DataArray:
    df = pd.read_csv(path, sep=r'\s+', header=None, skiprows=1).set_index(0)
    df.index.name = 'time'
    df.index *= 0.5 / 20.67  # convert a.u. to fs
    energy = df.loc[:, 4:]
    nstates = len(energy.columns)
    return (
        xr.Dataset.from_dataframe(energy)
        .to_array('state')
        .assign_coords(state=np.arange(1, nstates + 1))
    )


def parse_log(f) -> xr.Dataset:
    # Read MD settings:
    #  &md
    # -------------------------------------------------------
    #  Initial state:              2
    #  Initialize random velocity  0
    #  Temperature (K):            300
    #  Step:                       4000
    #  Dt (au):                    20.67

    # Don't use .startswith(), because there's also "&md velocity control"
    while not next(f).strip() == '&md':
        pass

    hline = next(f)  # Assertions may be skipped, so must not have side-effects
    assert hline.startswith('---')

    nsteps: int | None = None
    delta_t: float | None = None
    for line in f:
        stripline = line.strip()
        if stripline.startswith("Step:"):
            nsteps = int(stripline.split()[-1])
        if stripline.startswith("Dt (au):"):
            delta_t = float(stripline.split()[-1]) * 0.5 / 20.67
        if stripline.startswith('---'):
            break

    if nsteps is None:
        raise ValueError("Could not read `nsteps` ('Step:' in '&md' block)")
    if delta_t is None:
        raise ValueError(
            "Could not read `delta_t` ('Dt (au):' in '&md' block)")

    # Read final settings:
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
    explicit_ts = np.full((nsteps,), -1, dtype=int)
    astate = np.full((nsteps), -1, dtype=int)
    forces = np.full((nsteps, nstates, natoms, 3), np.nan)
    atXYZ = np.full((nsteps, natoms, 3), np.nan)
    atNames = np.full((natoms), '', dtype=str)
    got_atNames = False
    veloc = np.full((nsteps, natoms, 3), np.nan)
    dcmat = np.full((nsteps, nstates, nstates), np.nan)

    ts_idx = -1
    end_msg_count = 0
    for line in f:
        # The start of a timestep
        # Iter:        1  Ekin =           0.1291084223229551 au T =   300.00 K dt =         20 CI:   3
        # Root chosen for geometry opt   2
        if line.startswith('  Iter:'):
            ts_idx += 1
            explicit_ts[ts_idx] = int(line.strip().split()[1])

            for _ in range(10):
                # Get active state
                # A surface hopping is not allowed
                # **
                # At state:   2
                line = next(f)
                if line.startswith('  At state'):
                    astate[ts_idx] = int(line.strip().split()[2])
                    break
                # A surface hopping event happened
                # **
                # From state:   2 to state:   3 *
                elif line.startswith('  From state'):
                    astate[ts_idx] = int(line.strip().split()[5])
                    break
            else:
                raise ValueError(f"No state info found for Iter: {ts_idx+1}")

        # Positions:
        #   &coordinates in Angstrom
        # -------------------------------------------------------------------------------
        # C          0.5765950000000000     -0.8169010000000000     -0.0775610000000000
        # C          1.7325100000000000     -0.1032670000000000      0.1707480000000000
        # -------------------------------------------------------------------------------
        if line.startswith('  &coordinates'):
            hline = next(f)
            assert hline.startswith('---')
            if got_atNames:
                for iatom in range(natoms):
                    atXYZ[ts_idx, iatom] = np.asarray(
                        next(f).strip().split()[1:], dtype=float
                    )
            else:
                for iatom in range(natoms):
                    content = next(f).strip().split()
                    atXYZ[ts_idx, iatom] = np.asarray(content[1:], dtype=float)
                    atNames[iatom] = str(content[0])
                    got_atNames = True

            hline = next(f)
            assert hline.startswith('---')

        # Velocities:
        #   &velocities in Bohr/au
        # -------------------------------------------------------------------------------
        # C          0.0003442000000000      0.0001534200000000     -0.0000597200000000
        # C         -0.0005580000000000      0.0003118300000000     -0.0000154900000000
        # -------------------------------------------------------------------------------
        if line.startswith('  &velocities'):
            hline = next(f)
            assert hline.startswith('---')
            for iatom in range(natoms):
                veloc[ts_idx, iatom] = np.asarray(
                    next(f).strip().split()[1:], dtype=float
                )
            hline = next(f)
            assert hline.startswith('---')

        # Forces:
        #   &gradient state               1 in Eh/Bohr
        # -------------------------------------------------------------------------------
        # C         -0.0330978534152795      0.0073099255379017      0.0082666356536386
        # C          0.0313629524413876      0.0196036465968827      0.0060952442704520
        # -------------------------------------------------------------------------------
        if line.startswith('  &gradient'):
            istate = int(line.strip().split()[2]) - 1
            hline = next(f)
            assert hline.startswith('---')
            for iatom in range(natoms):
                forces[ts_idx, istate, iatom] = np.asarray(
                    next(f).strip().split()[1:], dtype=float
                )
            hline = next(f)
            assert hline.startswith('---')

        # Derivative coupling matrix:
        #  &derivative coupling matrix
        # -------------------------------------------------------------------------------
        #       0.0000000000000000       0.0000000000000004      -0.0000000000000001
        #      -0.0000000000000004       0.0000000000000000       0.0000000000000003
        #       0.0000000000000001      -0.0000000000000003       0.0000000000000000
        # -------------------------------------------------------------------------------
        if line.startswith('  &derivative coupling matrix'):
            hline = next(f)
            assert hline.startswith('---')
            for istate1 in range(nstates):
                dcmat[ts_idx, istate1] = np.asarray(
                    next(f).strip().split(), dtype=float
                )
            hline = next(f)
            assert hline.startswith('---')

        # Surface hopping information at the end of each timestep:
        #  &surface hopping information
        # -------------------------------------------------------
        #
        #     Random number:             0.15725129
        #     Accumulated probability:   0.00000000
        #     state mult  level   probability
        #     1     1     1       0.00000000
        #     2     1     2       0.00000000
        #     3     1     3       0.00000000
        #
        #
        # -------------------------------------------------------
        if line.startswith('  &surface hopping information'):
            hline = next(f)
            assert hline.startswith('---')
            # We don't currently parse this
            while not next(f).startswith('---'):
                pass

        # Completion indicator:
        # Nonadiabatic Molecular Dynamics End:  2025-04-13 01:12:26 Total:     0 days    15 hours    59 minutes    20 seconds
        if line.startswith('Nonadiabatic Molecular Dynamics End:'):
            end_msg_count += 1

    if end_msg_count > 1:
        raise ValueError(
            'Completion message "Nonadiabatic Molecular Dynamics End:" appeared '
            f"{end_msg_count} times"
        )

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
            'max_ts': explicit_ts.max(),
            # 'real_tmax': real_tmax,
            'delta_t': delta_t,
            'completed': end_msg_count == 1,
        },
    )
