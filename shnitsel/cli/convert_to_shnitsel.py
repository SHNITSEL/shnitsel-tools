import argparse
import logging
import pathlib
import sys

import shnitsel
from shnitsel.data.shnitsel_db.db_compound_group import CompoundInfo
from shnitsel.data.shnitsel_db_format import (
    MetaInformation,
    ShnitselDB,
    build_shnitsel_db,
)


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
        "--est_level",
        "-est",
        required=True,
        type=str,
        default=None,
        help="Level of applied Electronic Structure Theory.",
    )

    argument_parser.add_argument(
        "--basis_set",
        "-basis",
        required=True,
        type=str,
        default=None,
        help="The basis set used for for calculations.",
    )

    argument_parser.add_argument(
        "--loglevel",
        "-log",
        type=str,
        default="warn",
        help="The log level, `error`, `warn`, `info`, `debug`. ",
    )

    argument_parser.add_argument(
        "-f--force_write",
        action="store_true",
        help="A flag to make the script override existing files instead of halting if the output path already exists.",
    )

    argument_parser.add_argument(
        "--force_sequential",
        action="store_true",
        help="A flag to force sequential execution of trajectory conversion. Defaults to False to allow for parallel import and conversion.",
    )

    args = argument_parser.parse_args()

    input_path = pathlib.Path(args.input_path)
    input_kind = args.kind
    input_path_pattern = args.pattern
    input_group = args.group_name
    input_compound = args.compound_name

    input_est_level = args.est_level
    input_basis_set = args.basis_set

    output_path = args.output_path
    loglevel = args.loglevel

    force_sequential = args.force_sequential
    force_write = args.force_write

    found_file_at_beginning = False

    logging.basicConfig()

    logging.getLogger().setLevel(logging._nameToLevel[loglevel.upper()])

    if output_path is None:
        output_path = input_path / (input_path.name + ".nc")
    else:
        output_path = pathlib.Path(output_path)
        if output_path.suffix != ".nc":
            output_path = output_path.parent / (output_path.name + ".nc")

    if not input_path.exists():
        logging.error(f"Input path {input_path} does not exist")
        sys.exit(1)

    if output_path.exists():
        if force_write:
            logging.warning(
                f"Conversion will overwrite {output_path}. Will procede because of set `--force` flag."
            )
            found_file_at_beginning = True
        else:
            logging.error(
                f"Conversion would override {output_path}. For safety reasons, we will not proceed."
            )
            sys.exit(1)

    trajectory = shnitsel.io.read(
        input_path,
        sub_pattern=input_path_pattern,
        concat_method="db",
        kind=input_kind,
        parallel=not force_sequential,
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

        meta_info = MetaInformation(
            est_level=input_est_level, theory_basis_set=input_basis_set
        )
        trajectory.apply_trajectory_setup_properties(meta_info)

        num_compounds = len(trajectory.children)
        list_compounds = [str(k) for k in trajectory.children.keys()]
        print(f"Number of compounds in trajectory: {num_compounds}")
        print(f"Present compounds: {list_compounds}")
        num_trajectories = len(trajectory.leaves)
        print(f"Number of Trajectories: {num_trajectories}")

        if output_path.exists() and not found_file_at_beginning:
            logging.warning(
                f"File at {output_path} was written while conversion was running."
            )

            output_parent = output_path.parent
            output_file_name = output_path.stem
            output_suffix = output_path.suffix
            curr_index = 0

            while True:
                alternative_path = output_parent / (
                    output_file_name + f"_{curr_index}" + output_suffix
                )

                if not alternative_path.exists():
                    shnitsel.io.write_shnitsel_file(trajectory, alternative_path)
                    logging.warning(
                        f"To avoid loss of data, we swapped the output path to: {alternative_path}"
                    )
                    break
        else:
            shnitsel.io.write_shnitsel_file(trajectory, output_path)

        print("Wrote resulting trajectory collection:")
        pprint(trajectory)
        sys.exit(0)


if __name__ == "__main__":
    main()
