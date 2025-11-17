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
    global num_traj_sharc_traj, num_frames_sharc_traj
    # ...
    input_paths = ["./tutorials/test_data/pyrai2md/traj/"]

    for path in input_paths:
        if os.path.exists(path):
            traj: list[Trajectory] = read(
                path, kind='pyrai2md', parallel=False, concat_method='list'
            )  # type: ignore

            for tr in traj:
                num_traj_sharc_traj += 1
                num_frames_sharc_traj += tr.sizes["time"]
        else:
            logging.warning(f"Skipping {path}")


logging.basicConfig()

logging.getLogger().setLevel(logging._nameToLevel["DEBUG".upper()])

duration_newton_x = timeit.timeit(benchmark_newton_x, number=100)
print(f"{num_traj_newton_x},{num_frames_newton_x},{duration_newton_x}")
duration_sharc_icond = timeit.timeit(benchmark_sharc_icond, number=100)
print(f"{num_traj_sharc_icond},{num_frames_sharc_traj},{duration_sharc_icond}")
duration_sharc_traj = timeit.timeit(benchmark_sharc_traj, number=100)
print(f"{num_traj_sharc_traj},{num_frames_sharc_traj},{duration_sharc_traj}")
duration_pyrai2md = timeit.timeit(benchmark_pyrai2md, number=100)
print(f"{num_traj_pyraim2d},{num_frames_pyraim2d},{duration_pyrai2md}")
