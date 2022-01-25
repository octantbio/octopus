"""High-level integration tests covering the entire pipeline."""
import os
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np


def all_equal_or_NA(df: pd.DataFrame, field1: str, field2: str):
    """Compares two fields of a data frame with 'NA' values considered equal."""
    return ((df[field1] == df[field2])
            | df[field1].isnull() & df[field2].isnull()).all()


def compare_aggregated_stats(expected_fn: str, actual_fn: str):
    """Checks the relevant fields of two 'aggregated-stats.tsv' files.

    Compares the reference identified by the DeNovo, the barcode, the barcode
    contamination status are identical and that the depth at each base pair 
    (LT_10 and LT_3) are relatively close.

    Args:
        expected_fn: Path to "aggregated-stats.tsv" with expected values.
        actual_fn: Path to an actual "aggregated-stats.tsv" from the pipeline. 

    Raises:
        AssertionError: one or more fields does not match.
    """

    expected_df = pd.read_csv(expected_fn, sep="\t")
    actual_df = pd.read_csv(actual_fn, sep="\t")
    test_df = expected_df.merge(actual_df,
                                on=["Plate", "Well", "Plate_Well"],
                                how="outer",
                                suffixes=("_expected", "_actual"))

    assert all_equal_or_NA(test_df, "DeNovo_Ref_expected", "DeNovo_Ref_actual")
    assert all_equal_or_NA(test_df, "bc1_expected", "bc1_actual")
    assert all_equal_or_NA(test_df, "BC_Contam_expected", "BC_Contam_actual")

    # Exact LT_10 / LT_3 values fluctuate due to mapping differences
    assert np.allclose(test_df["LT_10_actual"],
                       test_df["LT_10_expected"],
                       rtol=0.01,
                       equal_nan=True)
    assert np.allclose(test_df["LT_3_actual"],
                       test_df["LT_3_expected"],
                       rtol=0.01,
                       equal_nan=True)


def test_pOK_barcode_full(tmp_path):
    """Tests the whole pipeline using the pOK_barcode dataset.

    Assumes that the current directory is the root of the OCTOPUS repository. 
    Intermediate files are created in a temporary directory.

    Args:
        tmp_path: Path to a new temporary directory.

    Raises:
        AssertionError: if the pipeline crashes or "aggregated-stats.tsv" does 
            not match the expected output.
    """
    project_dir = Path(os.getcwd())
    make_path = project_dir / "Makefile"

    # Create the data folder and symlink to it
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "pOK_barcode_test").symlink_to(
        project_dir / "test" / "pOK_barcode_test", target_is_directory=True)

    # HACK: Makefile uses relative paths, so we link to the src as a workaround
    (tmp_path / "src").symlink_to(project_dir / "src")

    output_file = "pipeline/pOK_barcode_test/aggregated-stats.tsv"
    os.chdir(tmp_path)
    args = ["make", "-f", str(make_path), str(output_file)]
    # pylint: disable=subprocess-run-check
    cmd = subprocess.run(args,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    assert cmd.returncode == 0
    compare_aggregated_stats(
        project_dir / "test" / "pOK_barcode-aggregated-stats.tsv",
        tmp_path / output_file)
