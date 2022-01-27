#!/usr/bin/env python3
"""Main script for running the computational pipeline.

This script acts as an entrypoint into OCTOPUS. This does not replace any
components of the pipeline, but simply spares users the trouble of dealing with
the Makefile interface. All pipeline-level arguments (like "--mem") can be
documented and dealt with here. Ultimately, this script simply calls `make`.
Example usage of this script:

    /path/to/main.py

In practice, this script is used as an entrypoint to our Docker container

Users who want more control should simply bypass this script and continue to 
work with `make` directly.
"""
import argparse
from pathlib import Path

OCTOPUS_PATH = Path(__file__).parent.parent


def check_illumina_dir(seq_dir: Path):
    """Checks that seq_dir is valid output dir from an Illumina sequencer."""

    # Check if it exists
    #  hint to user that they may have messed up their docker mounting
    # Check that it is readable
    # Check for all of the requisite files
    #  perhaps print out a tree in the case of an issue?
    # Check that a SampleSheet is present
    # Warn on no FASTQ


def main(seq_dir, to):
    """Starts an OCTOPUS run based on the given sequence directory."""

    check_illumina_dir(seq_dir)


DESCRIPTION = """
🐙 OCTOPUS
A plasmid sequencing pipeline from Octant Bio (https://octant.bio).
"""

parser = argparse.ArgumentParser(
    prog="octantbio/octopus",
    description=DESCRIPTION,
    epilog="Visit https://github.com/octantbio/octopus for further help.",
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    "seq_dir",
    metavar="ILLUMINA_SEQ_DIR",
    help="Path to an output directory from an Illumina sequencer")

# This is easy to do with `make`, but assigning names to pipeline stages in the
# future could be nice.
parser.add_argument(
    "-t",
    "--to",
    metavar="FILE",
    help=
    "Run pipeline up to and including 'FILE' (default: 'aggregated-stats.tsv')",
    default="aggregated-stats.tsv")
# TODO: add a 'from' argument here

if __name__ == "__main__":
    args = parser.parse_args()
    main(**args)
