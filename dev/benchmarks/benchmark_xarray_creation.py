import argparse
import sys
import timeit

from shnitsel.data.shnitsel_db_format import (
    build_shnitsel_db,
)
from shnitsel.io.shared.trajectory_setup import create_initial_dataset

from types import ModuleType, FunctionType
from gc import get_referents
from pprint import pprint

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

    num_steps = [0, 1, 10, 100, 10000]
    num_atoms = [0, 10, 100, 1000]
    num_states = [0, 2, 4, 6, 8]

    num_steps = [1]
    num_atoms = [12]
    num_states = [2]

    reps = 100

    stats = {}

    for nsteps in num_steps:
        for natoms in num_atoms:
            for statecount in num_states:
                print(nsteps, natoms, statecount)

                def create_set():
                    traj = create_initial_dataset(
                        nsteps, statecount, natoms, 'sharc', None
                    )
                    del traj

                duration_create = timeit.timeit(create_set, number=reps)
                stats[(nsteps, natoms, statecount)] = duration_create

                print(f"{stats[(nsteps, natoms, statecount)]=}")

    print(stats)

    sys.exit(0)


if __name__ == "__main__":
    main()
