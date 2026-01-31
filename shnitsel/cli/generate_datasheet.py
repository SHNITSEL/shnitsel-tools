import argparse
from dataclasses import asdict
import logging
import pathlib
import sys

import shnitsel
from shnitsel.data.dataset_containers.trajectory import (
    MetaInformation,
)
from shnitsel.data.tree import (
    complete_shnitsel_tree,
    ShnitselDB,
    CompoundInfo,
    GroupInfo,
)
from shnitsel.data.tree.compound import CompoundGroup
from shnitsel.data.tree.node import TreeNode
from shnitsel.vis.datasheet.datasheet import Datasheet


def main():
    argument_parser = argparse.ArgumentParser(
        sys.argv[0],
        f"{sys.argv[0]} <input_path> <output_path> [OPTIONS]",
        description="Script to read in an individual trajectory, a shnitsel-style file, an ASE database"
        "or a directory containing multiple sub-trajectories and generate a Shnitsel-Datasheet for it.\n\n"
        "Supports the reading of NewtonX, SHARC (ICOND and TRAJ) and PyrAI2md raw input files.",
    )

    argument_parser.add_argument(
        "input_path",
        help="The path to the input directory to read. Can point to an individual trajectory, an ASE db, a shnitsel `.nc` file or a parent directory of multiple trajectories.",
    )

    argument_parser.add_argument(
        "output_path",
        help="The path to place the fully generated datasheet into. Should have a `.pdf` extension or will be extended by `.pdf`",
    )

    argument_parser.add_argument(
        "-p",
        "--pattern",
        help="A glob pattern to use to identify subdirectories from which trajectories should be read. E.g. `TRAJ_*`.",
    )

    argument_parser.add_argument(
        "-c",
        "--compound_name",
        default=None,
        type=str,
        help="The name of the compound group for which to generate the datasheet for."
        "If using raw input from non-shnitsel formats, this instead specifies what name the input should have.",
    )

    argument_parser.add_argument(
        "-g",
        "--group_name",
        default=None,
        type=str,
        help="If the compound name is set and the input has a level of grouping underneath the chosen `compound_name`, this can pick a group within the data hierarchy to limit the scope of the datasheet",
    )

    argument_parser.add_argument(
        "--kind",
        "-k",
        required=False,
        type=str,
        default=None,
        help="Optionally an indication of the kind of input you are trying to read, `shnitsel`, `sharc`, `newtonx`, `pyrai2md`, `ase`. Will be guessed based on directory contents if not provided. If not set, the `read()` operation may fail if ambiguous trajectory formats are found within the folder.",
    )


    argument_parser.add_argument(
        "--charge",
        "-ch",
        required=False,
        type=int,
        default=None,
        help="Optional parameter to specify the charge of your imported molecule if it has not been set in the input. Must be an integer and is in units of electron charges (e)",
    )

    # argument_parser.add_argument(
    #     "--est_level",
    #     "-est",
    #     required=True,
    #     type=str,
    #     default=None,
    #     help="Level of applied Electronic Structure Theory.",
    # )

    # argument_parser.add_argument(
    #     "--basis_set",
    #     "-basis",
    #     required=True,
    #     type=str,
    #     default=None,
    #     help="The basis set used for for calculations.",
    # )

    argument_parser.add_argument(
        "--loglevel",
        "-log",
        type=str,
        default="warn",
        help="The log level, `error`, `warn`, `info`, `debug`. ",
    )

    argument_parser.add_argument(
        "-f",
        "--force_write",
        action="store_true",
        help="A flag to make the script override existing at the output position instead of appending a sequential number.",
    )

    # argument_parser.add_argument(
    #     "--force_sequential",
    #     action="store_true",
    #     help="A flag to force sequential execution of trajectory conversion. Defaults to False to allow for parallel import and conversion.",
    # )

    args = argument_parser.parse_args()

    input_path = pathlib.Path(args.input_path)
    output_path = pathlib.Path(args.output_path)
    input_kind = args.kind
    input_path_pattern = args.pattern
    input_group: str | None = args.group_name
    input_compound: str | None = args.compound_name

    charge = args.charge

    # input_est_level = args.est_level
    # input_basis_set = args.basis_set

    # output_path = args.output_path
    loglevel = args.loglevel

    # force_sequential = args.force_sequential
    force_write = args.force_write

    found_file_at_beginning = False

    logging.basicConfig()

    logging.getLogger().setLevel(logging._nameToLevel[loglevel.upper()])

    # if output_path is None:
    #     output_path = input_path / (input_path.name + ".pdf")
    # else:
    #     output_path = pathlib.Path(output_path)
    if output_path.suffix != ".pdf":
        output_path = output_path.parent / (output_path.name + ".pdf")

    if not input_path.exists():
        logging.error(f"Input path {input_path} does not exist")
        sys.exit(1)

    if output_path.exists():
        found_file_at_beginning = True
        if force_write:
            logging.warning(
                f"Conversion will overwrite {output_path}. Will procede because of set `--force` flag."
            )
        else:
            logging.warning(
                f"Conversion would override {output_path} and `--force` flag is not set."
            )
            # Try and find a collision-free name
            i = 1
            while True:
                tmp_filename = output_path.name
                tmp_new_path = tmp_filename[:-4] + f"_{i}.pdf"
                tmp_path = output_path.parent / tmp_new_path

                if not tmp_path.exists():
                    logging.warning(f"Changed output path to {tmp_path}.")
                    break
                else:
                    i += 1

    tree = shnitsel.io.read(
        input_path,
        sub_pattern=input_path_pattern,
        concat_method="db",
        kind=input_kind,
        parallel=True,  # not force_sequential,
    )

    if tree is None:
        logging.error("Input failed to load.")
        sys.exit(1)
    elif isinstance(tree, list):
        logging.error(
            "Imported trajectories failed to merge. Numbers of atoms or numbers of states differ. \n"
            "Please restrict your loading to a subset of trajectories with consistent parameters."
        )
        sys.exit(1)
    else:
        if charge is not None:
            tree =tree.set_charge(charge)
        if isinstance(tree, TreeNode):
            if not isinstance(tree, ShnitselDB):
                tree = complete_shnitsel_tree(tree)

            num_compounds = len(tree.children)
            print(f"Number of compounds in input: {num_compounds}")
            if input_compound:
                if "unknown" in tree.compounds:
                    target = tree.set_compound_info(input_compound)[input_compound]
                else:
                    if input_compound not in tree.compounds:
                        logging.error(
                            "Could not restrict analysis to compound %s. Only compounds in the input are: %s",
                            input_compound,
                            tree.compounds,
                        )
                        sys.exit(1)

                    target = tree[input_compound]
                # We now have a pointer to the compound_group level
                assert isinstance(target, CompoundGroup)
                if input_group:
                    if input_group not in target.subgroups:
                        logging.error(
                            "Could not restrict analysis to group %s. Only groups in the input under compound %s are: %s",
                            input_group,
                            input_compound,
                            tree.compounds,
                        )
                        sys.exit(1)
                    target = target[input_group]
            else:
                target = tree
            num_trajectories = len(list(target.collect_data()))
            print(f"Number of Trajectories in selection: {num_trajectories}")
        else:
            target = tree

        datasheet = Datasheet(target)
        print("Datasheet is being generated... (This may take a moment)")
        res = datasheet.plot(path=output_path)
        print(f"Datasheet has been written to {output_path}")
        sys.exit(0)


if __name__ == "__main__":
    main()
