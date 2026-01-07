import argparse
import logging
import pathlib
import sys

import shnitsel
from shnitsel.data.tree import ShnitselDB, complete_shnitsel_tree, tree_merge
from pprint import pprint


def main():
    argument_parser = argparse.ArgumentParser(
        sys.argv[0],
        f"{sys.argv[0]} <input_path1> <input_path2> [-o <output_path>] [OPTIONS]",
        description="Script to merge multiple individual Shnitsel files into one shared Shnitsel file.",
    )

    argument_parser.add_argument(
        "input_path",
        nargs='+',
        help="The paths to the input shnitsel files. If less than two are provided, nothing is done.",
    )

    argument_parser.add_argument(
        "-o",
        "--output_path",
        default=None,
        type=str,
        help="The path to put the converted shnitsel file to. if not provided will be the base name of the directory input_path is pointing to extended with `.nc` suffix. Should end on `.nc` or will be extended with `.nc`",
    )

    argument_parser.add_argument(
        "--loglevel",
        "-log",
        type=str,
        default="warn",
        help="The log level, `error`, `warn`, `info`, `debug`. ",
    )

    args = argument_parser.parse_args()

    input_paths = args.input_path

    output_path = args.output_path
    loglevel = args.loglevel

    if len(input_paths) < 2:
        logging.error("At least 2 Input paths need to be provided for merge.")

    logging.basicConfig()

    logging.getLogger().setLevel(logging._nameToLevel[loglevel.upper()])

    if output_path is None:
        path0 = pathlib.Path(input_paths[0])
        output_path = path0.parent / (path0.name + ".nc")
    else:
        output_path = pathlib.Path(output_path)
        if output_path.suffix != ".nc":
            output_path = output_path.parent / (output_path.name + ".nc")

    if output_path.exists():
        logging.error(
            f"Conversion would override {output_path}. For safety reasons, we will not proceed."
        )
        sys.exit(1)

    merged_tree = None

    tree_inputs: list[ShnitselDB] = []

    for input_path in input_paths:
        input_path = pathlib.Path(input_path)
        if not input_path.exists():
            logging.error(f"Input path {input_path} does not exist")
            sys.exit(1)

        trajectory = shnitsel.io.read(input_path, concat_method="db", kind='shnitsel')

        if trajectory is None:
            logging.error("Trajectory failed to load.")
            sys.exit(1)
        elif isinstance(trajectory, list):
            logging.error(
                "Trajectories failed to merge. Numbers of atoms or numbers of states differ. Please restrict your loading to a subset of trajectories with consistent parameters."
            )
            sys.exit(1)
        else:
            # pprint(trajectory)
            if not isinstance(trajectory, ShnitselDB):
                tree = complete_shnitsel_tree(trajectory)
            else:
                tree = trajectory

            tree_inputs.append(tree)

    merged_tree = tree_merge(*tree_inputs)

    if merged_tree is None:
        logging.error("Merge failed")
        sys.exit(1)

    num_compounds = len(merged_tree.children)
    list_compounds = [str(k) for k in merged_tree.children.keys()]
    num_trajectories = len(list(merged_tree.collect_data()))
    print(f"Number of compounds in trajectory: {num_compounds}")
    print(f"Present compounds: {list_compounds}")
    print(f"Number of Trajectories: {num_trajectories}")

    print("Resulting hierarchical structure:")
    shnitsel.io.write_shnitsel_file(merged_tree, output_path)
    pprint(merged_tree)
    sys.exit(0)


if __name__ == "__main__":
    main()
