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
import tempfile
import os
import subprocess
from pathlib import Path

OCTOPUS_PATH = Path(__file__).parent.parent


def check_illumina_dir(seq_dir: Path):
    """Checks that seq_dir is valid output dir from an Illumina sequencer."""

    # TODO:
    # Check if it exists
    #  hint to user that they may have messed up their docker mounting
    #  Warn on relative paths (~, ..)
    # Check that it is readable (permissions)
    # Check for all of the requisite files
    #  perhaps print out a tree in the case of an issue?
    # Check that a SampleSheet is present
    # Warn on no input.fa


def setup_run_dir(run_dir: Path, seq_dir: Path, out_dir: Path) -> Path:
    """Creates a run directory in the format that the Makefile anticipates.

    To allow for more flexible workflows than what `make` allows, this create
    a 'mock' directory with symlinks imitidating the correct file structure.

    Args:
        run_dir: Path to existing directory to perform run in.
        seq_dir: Path to input directory
        out_dir: Path to output directory

    Raises:
        FileNotFoundError: run_dir, seq_dir, or out_dir do not exist
    """

    (run_dir / "data").mkdir(exist_ok=True)
    (run_dir / "data" / seq_dir.name).symlink_to(seq_dir.resolve())
    (run_dir / "pipeline").mkdir()
    (run_dir / "pipeline" / out_dir.name).symlink_to(out_dir.resolve())

    # HACK: Makefile uses relative paths, so we link to the src as a workaround
    (run_dir / "src").symlink_to(OCTOPUS_PATH / "src")


def start_run(seq_dir: Path, out_dir: Path):
    """Starts an OCTOPUS run based with given sequence directory as input."""
    # Warn if running as root?

    check_illumina_dir(seq_dir)

    run_dir = Path(tempfile.mkdtemp())
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_run_dir(run_dir, seq_dir, out_dir)
    os.chdir(run_dir)

    # yapf: disable
    args = [
        "make",
        # point to Makefile in OCTOPUS directory
        "-f", f"{OCTOPUS_PATH}/Makefile",
        # `make` recipe
        # TODO: make this stop at the appropriate file
        f"pipeline/{out_dir.name}/aggregated-stats.tsv"
    ]
    # yapf: enable

    # TODO: tee this output into a file in the output folder
    # pylint: disable=subprocess-run-check
    pipeline_proc = subprocess.run(args, stderr=subprocess.STDOUT)

    return pipeline_proc.returncode


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
    help="Path to an output directory from an Illumina sequencer",
    type=Path)

parser.add_argument("-o",
                    "--out-dir",
                    metavar="DIR",
                    help="Output directory for resulting files "
                    "(default: 'pipeline/[ILLUMINA_SEQ_DIR]')")

# It's easier to refer to intermediate files with `make`, but assigning names to
# pipeline stages in the future could be nice.
# TODO: add 'start-from' and 'stop-at'

if __name__ == "__main__":
    args = parser.parse_args()
    if args.out_dir is None:
        args.out_dir = Path(f"pipeline/{args.seq_dir}")
    returncode = start_run(**vars(args))
    exit(returncode)
