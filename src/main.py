#!/usr/bin/env python3
"""Main script for running the computational pipeline.

This script acts as an entrypoint into the OCTOPUS docker image. This does not
replace any components of the pipeline, but simply spares users the trouble of
dealing with the Makefile interface. All pipeline-level arguments (like "--mem")
can be documented and dealt with here. Ultimately, this script simply calls
`make`.
Example usage of this script:

    /path/to/main.py

Users who want more control should simply bypass this script and continue to
work with `make` directly.
"""
import argparse
import tempfile
import os
import subprocess
from pathlib import Path

OCTOPUS_PATH = Path(__file__).parent.parent


def check_illumina_dir(seq_dir: Path) -> bool:
    """Checks that seq_dir is valid output dir from an Illumina sequencer."""

    try:
        seq_dir_exists = seq_dir.exists()
    except PermissionError as pe:
        print(f"I can't find '{seq_dir}' because I don't have permission to "
              f"access '{pe.filename}'.")
        return False

    if not seq_dir_exists:
        print(f"I couldn't find a folder / directory with name '{seq_dir}'.\n")

        if ".." in str(seq_dir):
            print("Warning: relative paths with '../' don't work well with "
                  "Docker.\n"
                  "You may need to change to an absolute path.")
        else:
            print("If you're using Docker, be sure to mount the appropriate "
                  "volumes with '-v'.\nFor example:\n"
                  f'    docker run -v "$(pwd):/data" octant/octopus {seq_dir}')
        return False
    if not seq_dir.is_dir():
        print(
            f"I was expecting a folder / directory, but input '{seq_dir}' is "
            "not a directory.\n")
        return False

    # smoke test to see if directory can be accessed
    try:
        next(seq_dir.iterdir())
    except StopIteration:
        pass
    except PermissionError:
        print(f"'{seq_dir}' exists but I can't seem to access it.\n"
              "Please ensure your permissions are appropriate.")
        return False

    try:
        if not (seq_dir / "SampleSheet.csv").exists():
            # could not find SampleSheet.csv inside illumina sequencing directory.
            print("I couldn't find a 'SampleSheet.csv' inside directory "
                  f"'{seq_dir}'.\n"
                  "You can find an example SampleSheet here:\n"
                  "https://github.com/octantbio/octopus/blob/master/test/"
                  "pOK_barcode_test/SampleSheet.csv")
            return False

        fastq_location = seq_dir / "Data" / "Intensities" / "BaseCalls"
        files = list(fastq_location.glob("*.fastq.gz"))
        if not files:
            print(
                f"I couldn't find any '.fastq.gz' files in '{fastq_location}'\n"
                "Please make sure your data have been demultiplexed and placed "
                "in that location.")
            return False

        if not (seq_dir / "input.fasta").exists():
            print(
                "Warning: I couldn't find a reference 'input.fasta' inside input "
                f"directory '{seq_dir}'.\n"
                "Pipeline will only perform de-novo alignment.\n")

    except PermissionError as perm_err:
        print(
            f"File '{perm_err.filename}' exists, but I don't have permission "
            "to access it.")

    return True


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

    print(
        "Checking that the provided sequencing directory has all the expected "
        "files...\n")
    if not check_illumina_dir(seq_dir):
        print("\nAborting pipeline due to the above error.")
        exit(1)

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
        # TODO: optionally make this stop at a user-defined file
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

parser.add_argument("seq_dir",
                    metavar="ILLUMINA_SEQ_DIR",
                    help="Path to an directory from an Illumina sequencer",
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
        print("No output directory provided, I'll put pipeline results in "
              f"'{args.out_dir}'.")
    returncode = start_run(**vars(args))
    exit(returncode)
