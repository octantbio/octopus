"""High-level integration tests covering the entire pipeline."""
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
