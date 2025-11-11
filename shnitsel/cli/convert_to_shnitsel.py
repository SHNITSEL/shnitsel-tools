import argparse
import logging
import pathlib
import sys

import shnitsel
from shnitsel.data.shnitsel_db.db_compound_group import CompoundInfo
from shnitsel.data.shnitsel_db_format import ShnitselDB, build_shnitsel_db


def main():
    argument_parser = argparse.ArgumentParser(
        sys.argv[0],
        f"{sys.argv[0]} <input_path> [OPTIONS]",
        description="Script to read in an individual trajectory or a directory containing multiple "
        "sub-trajectories and convert them into a shnitsel-style file format.\n\n"
        "Currently supports reading of NewtonX, SHARC (ICOND and TRAJ) and PyrAI2md files.",
    )

    argument_parser.add_argument(
        "input_path",
        help="The path to the input directory to read. Can point to an individual trajectory or a parent directory of multiple trajectories.",
    )

    argument_parser.add_argument(
        "-p",
        "--pattern",
        help="A glob pattern to use to identify subdirectories from which trajectories should be read. E.g. `TRAJ_*`.",
    )

    argument_parser.add_argument(
        "-o",
        "--output_path",
        default=None,
        type=str,
        help="The path to put the converted shnitsel file to. if not provided will be the base name of the directory input_path is pointing to extended with `.nc` suffix. Should end on `.nc` or will be extended with `.nc`",
    )

    argument_parser.add_argument(
        "-c",
        "--compound_name",
        default=None,
        type=str,
        help="The name of the compound group to add the read trajectories to. E.g. `R02` or `I01`.",
    )

    argument_parser.add_argument(
        "-g",
        "--group_name",
        default=None,
        type=str,
        help="If set, all read trajectories will be added to a group of this name. This allows for differentiation between different trajectories within the same compound.",
    )

    argument_parser.add_argument(
        "--kind",
        "-k",
        required=False,
        type=str,
        default=None,
        help="Optionally an indication of the kind of trajectory you want to read, `shnitsel`, `sharc`, `newtonx`, `pyrai2md`. Will be guessed based on directory contents if not provided. If not set, the conversion may fail if ambiguous trajectory formats are found within the folder.",
    )

    argument_parser.add_argument(
        "--loglevel",
        "-log",
        type=str,
        default="warn",
        help="The log level, `error`, `warn`, `info`, `debug`. ",
    )

    args = argument_parser.parse_args()

    input_path = pathlib.Path(args.input_path)
    input_kind = args.kind
    input_path_pattern = args.pattern
    input_group = args.group_name
    input_compound = args.compound_name

    output_path = args.output_path
    loglevel = args.loglevel

    logging.basicConfig()

    logging.getLogger().setLevel(logging._nameToLevel[loglevel.upper()])

    if output_path is None:
        output_path = input_path / (input_path.name + ".nc")
    else:
        output_path = pathlib.Path(output_path)
        if output_path.suffix != ".nc":
            output_path = output_path.parent / (output_path.name + ".nc")

    if output_path.exists():
        logging.error(
            f"Conversion would override {output_path}. For safety reasons, we will not proceed."
        )
        sys.exit(1)

    if not input_path.exists():
        logging.error(f"Input path {input_path} does not exist")
        sys.exit(1)

    trajectory = shnitsel.io.read(
        input_path, sub_pattern=input_path_pattern, concat_method="db", kind=input_kind
    )

    from pprint import pprint

    if trajectory is None:
        logging.error("Trajectory failed to load.")
        sys.exit(1)
    elif isinstance(trajectory, list):
        logging.error(
            "Trajectories failed to merge. Numbers of atoms or numbers of states differ. Please restrict your loading to a subset of trajectories with consistent parameters."
        )
        sys.exit(1)
    else:
        if not isinstance(trajectory, ShnitselDB):
            trajectory = build_shnitsel_db(trajectory)

        compound_info = CompoundInfo()
        if input_compound:
            compound_info.compound_name = input_compound
            trajectory = trajectory.set_compound_info(compound_info=compound_info)

        if input_group:
            trajectory = trajectory.add_trajectory_group(input_group)

        num_compounds = len(trajectory.children)
        list_compounds = [str(k) for k in trajectory.children.keys()]
        print(f"Number of compounds in trajectory: {num_compounds}")
        print(f"Present compounds: {list_compounds}")

        print(f"Resulting trajectory:")
        pprint(trajectory)
        shnitsel.io.write_shnitsel_file(trajectory, output_path)
        sys.exit(0)


if __name__ == "__main__":
    main()
