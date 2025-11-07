import argparse
import logging
import pathlib
import sys

import shnitsel


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
        "-o", "--output_path",
        default=None,
        type=str,
        help="The path to put the converted shnitsel file to. if not provided will be the base name of the directory input_path is pointing to extended with `.nc` suffix. Should end on `.nc` or will be extended with `.nc`",
    )

    argument_parser.add_argument(
        "--kind",
        "-k",
        required=False,
        type=str,
        default=None,
        help="Optionally an indication of the kind of trajectory you want to read, `shnitsel`, `sharc`, `newtonx`, `pyrai2md`. Wil be guessed based on directory contents if not provided.",
    )

    args = argument_parser.parse_args()

    input_path = pathlib.Path(args.input_path)
    input_kind = args.kind

    output_path = args.output_path

    if output_path is None:
        output_path = input_path / (input_path.name + ".nc")
    else:
        output_path = pathlib.Path(output_path)
        if output_path.suffix != ".nc":
            output_path = output_path.parent / (output_path.name + ".nc")

    if output_path.exists():
        logging.error(
            "Conversion would override {output_path}. For safety reasons, we will not proceed."
        )
        sys.exit(1)

    if not input_path.exists():
        logging.error(f"Input path {input_path} does not exist")
        sys.exit(1)

    trajectory = shnitsel.io.read(input_path, kind=input_kind)
    if trajectory is None:
        logging.error("Trajectory failed to load.")
        sys.exit(1)
    elif isinstance(trajectory, list):
        logging.error(
            "Trajectories failed to merge. Numbers of atoms or numbers of states differ. Please restrict your loading to a subset of trajectories with consistent parameters."
        )
        sys.exit(1)
    else:
        shnitsel.io.write_shnitsel_file(trajectory, output_path)
        sys.exit(0)


if __name__ == "__main__":
    main()
