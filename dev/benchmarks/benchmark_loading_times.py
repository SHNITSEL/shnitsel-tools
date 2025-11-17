import logging
import os
import timeit

from shnitsel.data.trajectory_format import Trajectory
from shnitsel.io import read

num_traj_newton_x = 0
num_frames_newton_x = 0
num_frames_sharc_icond = 0
num_traj_sharc_icond = 0
num_frames_sharc_traj = 0
num_traj_sharc_traj = 0
num_traj_pyraim2d = 0
num_frames_pyraim2d = 0


def benchmark_newton_x():
    global num_traj_newton_x, num_frames_newton_x
    # ...
    input_paths = [
        "./tutorials/test_data/newtonx/test_R02/",
    ]

    for path in input_paths:
        if os.path.exists(path):
            traj: list[Trajectory] = read(
                path, kind='newtonx', parallel=False, concat_method='list'
            )  # type: ignore

            for tr in traj:
                num_traj_newton_x += 1
                num_frames_newton_x += tr.sizes["time"]
        else:
            logging.warning(f"Skipping {path}")


def benchmark_sharc_icond():
    global num_traj_sharc_icond, num_frames_sharc_icond
    # ...
    input_paths = ["./tutorials/test_data/sharc/iconds_butene/"]

    for path in input_paths:
        if os.path.exists(path):
            traj: list[Trajectory] = read(
                path, kind='sharc', parallel=False, concat_method='list'
            )  # type: ignore

            for tr in traj:
                num_traj_sharc_icond += 1
                num_frames_sharc_icond += tr.sizes["time"]
        else:
            logging.warning(f"Skipping {path}")


def benchmark_sharc_traj():
    global num_traj_sharc_traj, num_frames_sharc_traj
    # ...
    input_paths = [
        "./tutorials/test_data/sharc/traj_butene/",
        "./tutorials/test_data/sharc/traj_I01/",
    ]

    for path in input_paths:
        if os.path.exists(path):
            traj: list[Trajectory] = read(
                path, kind='sharc', parallel=False, concat_method='list'
            )  # type: ignore

            for tr in traj:
                num_traj_sharc_traj += 1
                num_frames_sharc_traj += tr.sizes["time"]
        else:
            logging.warning(f"Skipping {path}")


def benchmark_pyrai2md():
    global num_traj_pyraim2d, num_frames_pyraim2d
    # ...
    input_paths = ["./tutorials/test_data/pyrai2md/"]

    for path in input_paths:
        if os.path.exists(path):
            traj: list[Trajectory] = read(
                path, kind='pyrai2md', parallel=False, concat_method='list'
            )  # type: ignore

            for tr in traj:
                num_traj_pyraim2d += 1
                num_frames_pyraim2d += tr.sizes["time"]
        else:
            logging.warning(f"Skipping {path}")


logging.basicConfig()

reps = 1000

# logging.getLogger().setLevel(logging._nameToLevel["DEBUG".upper()])

with open("read_benchmark.dat", "w") as out:
    duration_newton_x = timeit.timeit(benchmark_newton_x, number=reps)
    out.write(
        f"Newtonx:\t{num_traj_newton_x},\t{num_frames_newton_x},\t{duration_newton_x}\n"
    )
    out.flush()
    duration_sharc_icond = timeit.timeit(benchmark_sharc_icond, number=reps)
    out.write(
        f"SHARC(icond):\t{num_traj_sharc_icond},\t{num_frames_sharc_icond},\t{duration_sharc_icond}\n"
    )
    out.flush()
    duration_sharc_traj = timeit.timeit(benchmark_sharc_traj, number=reps)
    out.write(
        f"SHARC(traj):\t{num_traj_sharc_traj},\t{num_frames_sharc_traj},\t{duration_sharc_traj}\n"
    )
    out.flush()
    duration_pyrai2md = timeit.timeit(benchmark_pyrai2md, number=reps)
    out.write(
        f"PyrAI2md:\t{num_traj_pyraim2d},\t{num_frames_pyraim2d},\t{duration_pyrai2md}\n"
    )
    out.flush()
