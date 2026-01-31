import logging
import os
import timeit

from shnitsel.data.dataset_containers import Frames, Trajectory
from shnitsel.io import read

if __name__ == '__main__':
    num_traj_newton_x = 0
    num_frames_newton_x = 0
    num_frames_sharc_icond = 0
    num_traj_sharc_icond = 0
    num_frames_sharc_traj = 0
    num_traj_sharc_traj = 0
    num_traj_pyraim2d = 0
    num_frames_pyraim2d = 0

    parallel_mode = False


    def benchmark_newton_x():
        global num_traj_newton_x, num_frames_newton_x
        # ...
        input_paths = [
            "./tutorials/test_data/newtonx/test_R02_a_v2.2/",
            "./tutorials/test_data/newtonx/test_R02_b/",
        ]

        for path in input_paths:
            if os.path.exists(path):
                traj = read(
                    path, kind='newtonx', parallel=parallel_mode,
                )  # type: ignore

                for tr in traj.collect_data():
                    num_traj_newton_x += 1
                    num_frames_newton_x += tr.sizes["time"]
            else:
                logging.warning(f"Skipping {path}")


    def benchmark_sharc_icond():
        global num_traj_sharc_icond, num_frames_sharc_icond
        # ...
        input_paths = [
            "./tutorials/test_data/sharc/iconds_butene/",
        ]

        for path in input_paths:
            if os.path.exists(path):
                traj = read(
                    path, kind='sharc', parallel=parallel_mode,
                )  # type: ignore

                for tr in traj.collect_data():
                    num_traj_sharc_icond += 1
                    num_frames_sharc_icond += tr.sizes["time"]
            else:
                logging.warning(f"Skipping {path}")


    def benchmark_sharc_traj():
        global num_traj_sharc_traj, num_frames_sharc_traj
        # ...
        input_paths = [
            "./tutorials/test_data/sharc/traj_butene_v2.1/",
            "./tutorials/test_data/sharc/traj_I01_v2.0/",
            "./tutorials/test_data/sharc/traj_I01_v3.0_triplets/",
            "./tutorials/test_data/sharc/traj_I01_v3.0_triplets_nacs_socs/",
            "./tutorials/test_data/sharc/traj_I01_v4.0/",
        ]

        for path in input_paths:
            if os.path.exists(path):
                traj = read(
                    path, kind='sharc', parallel=parallel_mode,
                )  # type: ignore

                for tr in traj.collect_data():
                    num_traj_sharc_traj += 1
                    num_frames_sharc_traj += tr.sizes["time"]
            else:
                logging.warning(f"Skipping {path}")


    def benchmark_pyrai2md():
        global num_traj_pyraim2d, num_frames_pyraim2d
        # ...
        input_paths = ["./tutorials/test_data/pyrai2md/traj_I01/"]

        for path in input_paths:
            if os.path.exists(path):
                traj = read(
                    path, kind='pyrai2md', parallel=parallel_mode,
                )  # type: ignore

                for tr in traj.collect_data():
                    num_traj_pyraim2d += 1
                    num_frames_pyraim2d += tr.sizes["time"]
            else:
                logging.warning(f"Skipping {path}")


    logging.basicConfig()
    reps = 10

    # logging.getLogger().setLevel(logging._nameToLevel["DEBUG".upper()])

    with open("read_benchmark.dat", "w") as out:
        duration_newton_x = timeit.timeit(benchmark_newton_x, number=reps)
        out.write(
            f"Newtonx:\t{num_traj_newton_x / reps},\t{num_frames_newton_x / reps},\t{duration_newton_x / reps},\t{duration_newton_x * 100 / num_frames_newton_x}\n"
        )
        out.flush()

        duration_sharc_icond = timeit.timeit(benchmark_sharc_icond, number=reps)
        out.write(
            f"SHARC(icond):\t{num_traj_sharc_icond / reps},\t{num_frames_sharc_icond / reps},\t{duration_sharc_icond / reps},\t{duration_sharc_icond * 100 / num_frames_sharc_icond}\n"
        )
        out.flush()

        duration_sharc_traj = timeit.timeit(benchmark_sharc_traj, number=reps)
        out.write(
            f"SHARC(traj):\t{num_traj_sharc_traj / reps},\t{num_frames_sharc_traj / reps},\t{duration_sharc_traj / reps},\t{duration_sharc_traj * 100 / num_frames_sharc_traj}\n"
        )
        out.flush()

        duration_pyrai2md = timeit.timeit(benchmark_pyrai2md, number=reps)
        out.write(
            f"PyrAI2md:\t{num_traj_pyraim2d / reps},\t{num_frames_pyraim2d / reps},\t{duration_pyrai2md / reps},\t{duration_pyrai2md * 100 / num_frames_pyraim2d}\n"
        )
        out.flush()

        parallel_mode = True
        num_traj_newton_x = 0
        num_frames_newton_x = 0
        num_frames_sharc_icond = 0
        num_traj_sharc_icond = 0
        num_frames_sharc_traj = 0
        num_traj_sharc_traj = 0
        num_traj_pyraim2d = 0
        num_frames_pyraim2d = 0
        
        duration_newton_x = timeit.timeit(benchmark_newton_x, number=reps)
        out.write(
            f"Newtonx (parallel):\t{num_traj_newton_x / reps},\t{num_frames_newton_x / reps},\t{duration_newton_x / reps},\t{duration_newton_x * 100 / num_frames_newton_x}\n"
        )
        out.flush()

        duration_sharc_icond = timeit.timeit(benchmark_sharc_icond, number=reps)
        out.write(
            f"SHARC(icond) (parallel):\t{num_traj_sharc_icond / reps},\t{num_frames_sharc_icond / reps},\t{duration_sharc_icond / reps},\t{duration_sharc_icond * 100 / num_frames_sharc_icond}\n"
        )
        out.flush()

        duration_sharc_traj = timeit.timeit(benchmark_sharc_traj, number=reps)
        out.write(
            f"SHARC(traj) (parallel):\t{num_traj_sharc_traj / reps},\t{num_frames_sharc_traj / reps},\t{duration_sharc_traj / reps},\t{duration_sharc_traj * 100 / num_frames_sharc_traj}\n"
        )
        out.flush()

        duration_pyrai2md = timeit.timeit(benchmark_pyrai2md, number=reps)
        out.write(
            f"PyrAI2md (parallel):\t{num_traj_pyraim2d / reps},\t{num_frames_pyraim2d / reps},\t{duration_pyrai2md / reps},\t{duration_pyrai2md * 100 / num_frames_pyraim2d}\n"
        )
        out.flush()
