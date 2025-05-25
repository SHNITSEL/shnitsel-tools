import os
from glob import glob
from ast import literal_eval
import xarray as xr
import pandas as pd
import numpy as np


def parse_md_energies(path):
    df = pd.read_csv(path, sep=r'\s+', header=None, skiprows=1).set_index(0)
    df.index.name = 'time'
    energy = df.loc[:, 4:]
    state = len(energy.columns)
    return (
        xr.Dataset.from_dataframe(energy)
        .to_array('state')
        .assign_coords(state=[1, 2, 3])
    )


def parse_log(f, nsteps):
    ## First we read the settings.
    #   &molecule
    # -------------------------------------------------------
    #   States:                     [3]
    #   Spin:                       [0]
    #   Interstates:                [[1, 2], [1, 3], [2, 3]]
    #   QMMM keyfile:               None
    #   QMMM xyzfile:               Input
    ##  etc.
    for line in f:
        if line.startswith('  &molecule'):
            # Found settings section.
            next(f)
            break
    else:  # if no break occured
        raise ValueError("No '&molecule' section found in log file.")

    settings = {}
    for line in f:
        parts = line.strip().split(':')
        if len(parts) != 2:
            raise ValueError(
                "Expected exactly one colon on each line of the '&molecule block"
            )
        k, v = parts
        settings[k] = v

    states = literal_eval(settings['States'])
    spin = literal_eval(settings['Spin'])
    if len(states) != 1 or spin[0] != 0:
        raise ValueError("Only handling singlets for now.")
    nstates = states[0]
    interstates = literal_eval(settings['Interstates'])

    ## Get the number of atoms from a separate settings block.
    #   Active atoms:     45
    for line in f:
        if line.startswith('  Active atoms'):
            natoms = int(line.strip().split()[1])

    # Set up numpy arrays
    astate = np.full((nsteps), -1, dtype=int)
    forces = np.full((nsteps, natoms, 3), np.nan)
    # time = np.full((nsteps), np.nan)
    # nacs = np.full((nsteps, math.comb(nstates, 2), natoms, 3), np.nan)

    for line in f:
        ## The start of a timestep
        # Iter:        1  Ekin =           0.1291084223229551 au T =   300.00 K dt =         20 CI:   3
        # Root chosen for geometry opt   2
        if line.startswith('  Iter:'):
            ts = int(line.strip().split()[1])

            for line in f:
                ## Get active state
                # A surface hopping is not allowed
                # **
                # At state:   2
                if line.startswith('  At state'):
                    astate[ts] = int(line.strip().split()[2])
                # A surface hopping event happened
                # **
                # From state:   2 to state:   3 *
                elif line.startswith('  From state'):
                    astate[ts] = int(line.strip().split()[5])

    print(astate)
    raise KeyboardInterrupt


def read_traj(traj_path):
    # with open(os.path.join(traj_path, 'RESULTS', 'nx.log')) as f:
    #     single_traj = parse_nx_log(f)

    # nsteps = single_traj.sizes['ts']

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
    single_traj['energy'] = energy

    return single_traj