import argparse
import sys

from shnitsel.io.shared.trajectory_setup import create_initial_dataset
from shnitsel.data.tree import complete_shnitsel_tree

from types import ModuleType, FunctionType
from gc import get_referents

# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType


def getsize(obj):
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: ' + str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


def main():
    argument_parser = argparse.ArgumentParser(
        sys.argv[0],
        f"{sys.argv[0]} ",
        description="Script to test the memory consumption of xarray Shnitsel Datasets.",
    )

    args = argument_parser.parse_args()

    num_steps = [0, 10, 100, 10000]
    num_atoms = [0, 10, 100, 1000]
    num_states = [0, 2, 4, 6, 8]

    stats = {}
    stats_db = {}

    for nsteps in num_steps:
        for natoms in num_atoms:
            for statecount in num_states:
                print(nsteps, natoms, statecount)
                traj, default_param = create_initial_dataset(
                    nsteps, statecount, natoms, 'sharc', None
                )
                stats[(nsteps, natoms, statecount)] = getsize(traj)
                print(f"{stats[(nsteps, natoms, statecount)]=}")
                db_traj = complete_shnitsel_tree(traj)
                stats_db[(nsteps, natoms, statecount)] = getsize(db_traj)
                print(f"{stats_db[(nsteps, natoms, statecount)]=}")
                del traj
                del db_traj

    base = stats[(0, 0, 0)]

    print("done")
    # print(stats)
    # print(stats_db)
    with open("./stats_xarray.dat", "w") as stats_out:
        stats_out.write("# num_steps\t# num_atoms\t# num_states\t# size [bytes]\n")
        for (steps, atoms, states), size in stats.items():
            stats_out.write(f"{steps}\t{atoms}\t{states}\t{size}\n")
    with open("./stats_db.dat", "w") as stats_out:
        stats_out.write("# num_steps\t# num_atoms\t# num_states\t# size [bytes]\n")
        for (steps, atoms, states), size in stats_db.items():
            stats_out.write(f"{steps}\t{atoms}\t{states}\t{size}\n")

    def scale_size(param, base2, states_scale, states_scale2):
        # print(param)
        nonlocal base
        m_scale = 1
        (steps, atoms, states) = param
        return base + m_scale * steps * atoms * (
            base2 + states_scale * states + states_scale2 * states * states
        )

    from scipy.optimize import curve_fit
    import numpy as np

    x = np.array(list(zip(*list(stats.keys()))))
    y = np.array([float(x) for x in stats.values()])

    popt, pcov = curve_fit(scale_size, x, y, p0=(0, 1, 1))

    print(popt)
    perr = np.sqrt(np.diag(pcov))
    print(perr)
    m_scale = 1
    base2, states_scale, states_scale2 = popt

    print(f"{base=} {m_scale=} {base2=} {states_scale=} {states_scale2=}")

    import matplotlib.pyplot as plt

    plt.clf()

    sys.exit(0)


if __name__ == "__main__":
    main()
