#! /usr/bin/env python
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# cmdoret, 20181214
"""
Simple Hi-C pipeline for generating and manipulating contact matrices.

usage:
    hicstuff [-hv] <command> [<args>...]

options:
    -h, --help                  shows the help
    -v, --version               shows the version

The subcommands are:
    convert         Convert Hi-C data between different formats.
    digest          Digest genome into a list of fragments.
    cutsite         Preprocess fastq files by digesting reads at religation site.
    distancelaw     Analyse and plot distance law.
    filter          Filters Hi-C pairs to exclude spurious events.
    iteralign       Iteratively aligns reads to a reference genome.
    missview        Preview missing Hi-C bins in based on the genome and read length.
    pipeline        Hi-C pipeline to generate contact matrix from fastq files.
    rebin           Bin the matrix and regenerate files accordingly.
    subsample       Bootstrap subsampling of contacts from a Hi-C map.
    view            Visualize a Hi-C matrix.    
"""

from docopt import docopt
from docopt import DocoptExit
import hicstuff.commands as commands
from hicstuff.version import __version__


def main():
    args = docopt(__doc__, version=__version__, options_first=True)
    # Retrieve the command to execute.
    command_name = args.pop("<command>").capitalize()

    # Retrieve the command arguments.
    command_args = args.pop("<args>")
    if command_args is None:
        command_args = {}
    # After 'popping' '<command>' and '<args>', what is left in the
    # args dictionary are the global arguments.

    # Retrieve the class from the 'commands' module.
    try:
        command_class = getattr(commands, command_name)
    except AttributeError:
        print("Unknown command.")
        raise DocoptExit()
    # Create an instance of the command.
    command = command_class(command_args, args)
    # Execute the command.
    command.execute()


if __name__ == "__main__":
    main()
