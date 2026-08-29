import gc
import os
import time
import traceback
from pathlib import Path
import bioframe as bf
import numpy as np
import pandas as pd
from ancIBD.IO.batch_run import to_ibd_df_batches


def find_short_ibd_overlapping_long(
    df_short,
    df_long,
    normalize_pair_order=True,
    return_matches=False,
    return_clean=False,
):
    """
    Find short IBD segments that overlap long IBD segments
    using physical base-pair coordinates.

    An overlap is identified only when the short and long segments
    belong to:
        1. the same individual pair;
        2. the same chromosome;
        3. overlapping StartBP-EndBP intervals.

    Bioframe uses half-open interval semantics:
        [StartBP, EndBP)

    Therefore, intervals that only touch at one endpoint are not
    considered overlapping.

    Parameters
    ----------
    df_short : pandas.DataFrame
        Short-IBD segment-level table.

        Required columns:
            iid1, iid2, ch, StartBP, EndBP

    df_long : pandas.DataFrame
        Long-IBD segment-level table.

        Required columns:
            iid1, iid2, ch, StartBP, EndBP

    normalize_pair_order : bool, default=True
        If True, treat A-B and B-A as the same individual pair.

    return_matches : bool, default=False
        If True, also return the complete short-long overlap table.

    return_clean : bool, default=False
        If True, also return the short-IBD table after removing
        all short segments overlapping any long IBD segment.

    Returns
    -------
    Default:
        df_short_overlap

    return_matches=True:
        df_short_overlap, matches

    return_clean=True:
        df_short_clean, df_short_overlap

    return_clean=True and return_matches=True:
        df_short_clean, df_short_overlap, matches

    Notes
    -----
    df_short_overlap
        Unique short IBD segments that overlap at least one
        long IBD segment.

    df_short_clean
        Short IBD segments remaining after all overlapping
        short segments have been removed.

    matches
        One row per short-long overlap relationship. One short
        segment may appear more than once when it overlaps more
        than one long segment.
    """

    # ------------------------------------------------------------
    # 1. Check required columns
    # ------------------------------------------------------------
    required = {
        "iid1",
        "iid2",
        "ch",
        "StartBP",
        "EndBP",
    }

    missing_short = required.difference(df_short.columns)
    missing_long = required.difference(df_long.columns)

    if missing_short:
        raise ValueError(
            "Short IBD table is missing required columns: "
            f"{sorted(missing_short)}"
        )

    if missing_long:
        raise ValueError(
            "Long IBD table is missing required columns: "
            f"{sorted(missing_long)}"
        )

    # ------------------------------------------------------------
    # 2. Handle empty input tables
    # ------------------------------------------------------------
    if df_short.empty:
        raise ValueError(
            "The short IBD table is empty. No segments to process. "
            "Please provide a non-empty df_short."
        )

    if df_long.empty:
        raise ValueError(
            "The long IBD table is empty. Cannot perform overlap analysis "
            "without anchor segments. Please provide a non-empty df_long."
        )

    # ------------------------------------------------------------
    # 3. Create working copies
    # ------------------------------------------------------------
    short = df_short.reset_index(drop=True).copy()
    long = df_long.reset_index(drop=True).copy()

    short["_short_id"] = np.arange(
        len(short),
        dtype=np.int64,
    )

    long["_long_id"] = np.arange(
        len(long),
        dtype=np.int64,
    )

    # ------------------------------------------------------------
    # 4. Standardize chromosome and BP coordinate types
    # ------------------------------------------------------------
    short["ch"] = short["ch"].astype(str)
    long["ch"] = long["ch"].astype(str)

    for table_name, table in [
        ("short", short),
        ("long", long),
    ]:
        # Convert BP coordinates to numeric.
        table["StartBP"] = pd.to_numeric(
            table["StartBP"],
            errors="raise",
        )

        table["EndBP"] = pd.to_numeric(
            table["EndBP"],
            errors="raise",
        )

        # Check missing coordinates.
        invalid_coordinates = (
            table["StartBP"].isna()
            | table["EndBP"].isna()
        )

        if invalid_coordinates.any():
            raise ValueError(
                f"{table_name.capitalize()} IBD table contains "
                f"{int(invalid_coordinates.sum())} rows with "
                "missing StartBP or EndBP."
            )

        # BP coordinates must represent integer positions.
        non_integer_coordinates = (
            (table["StartBP"] % 1 != 0)
            | (table["EndBP"] % 1 != 0)
        )

        if non_integer_coordinates.any():
            raise ValueError(
                f"{table_name.capitalize()} IBD table contains "
                f"{int(non_integer_coordinates.sum())} rows with "
                "non-integer StartBP or EndBP coordinates."
            )

        # Convert to numpy int64, as required by bioframe.
        table["StartBP"] = table["StartBP"].astype(np.int64)
        table["EndBP"] = table["EndBP"].astype(np.int64)

        # Check interval direction.
        invalid_intervals = (
            table["EndBP"] < table["StartBP"]
        )

        if invalid_intervals.any():
            raise ValueError(
                f"{table_name.capitalize()} IBD table contains "
                f"{int(invalid_intervals.sum())} rows with "
                "EndBP < StartBP."
            )

        # Bioframe intervals should have positive width.
        zero_width_intervals = (
            table["EndBP"] == table["StartBP"]
        )

        if zero_width_intervals.any():
            raise ValueError(
                f"{table_name.capitalize()} IBD table contains "
                f"{int(zero_width_intervals.sum())} zero-width "
                "intervals where StartBP == EndBP."
            )

    # ------------------------------------------------------------
    # 5. Normalize individual-pair orientation
    # ------------------------------------------------------------
    if normalize_pair_order:
        short_iid1 = short["iid1"].astype(str).to_numpy()
        short_iid2 = short["iid2"].astype(str).to_numpy()

        long_iid1 = long["iid1"].astype(str).to_numpy()
        long_iid2 = long["iid2"].astype(str).to_numpy()

        short["pair1"] = np.minimum(
            short_iid1,
            short_iid2,
        )

        short["pair2"] = np.maximum(
            short_iid1,
            short_iid2,
        )

        long["pair1"] = np.minimum(
            long_iid1,
            long_iid2,
        )

        long["pair2"] = np.maximum(
            long_iid1,
            long_iid2,
        )

    else:
        short["pair1"] = short["iid1"].astype(str)
        short["pair2"] = short["iid2"].astype(str)

        long["pair1"] = long["iid1"].astype(str)
        long["pair2"] = long["iid2"].astype(str)

    # ------------------------------------------------------------
    # 6. Find all short-long BP-coordinate overlaps
    # ------------------------------------------------------------
    matches = bf.overlap(
        short,
        long,
        how="inner",
        cols1=("ch", "StartBP", "EndBP"),
        cols2=("ch", "StartBP", "EndBP"),
        on=["pair1", "pair2"],
        suffixes=("_short", "_long"),
        return_input=True,
        return_index=True,
        return_overlap=True,
        ensure_int=True,
    )

    # ------------------------------------------------------------
    # 7. Handle the situation in which no overlaps are found
    # ------------------------------------------------------------
    if matches.empty:
        df_short_overlap = df_short.iloc[0:0].copy()
        df_short_clean = df_short.copy().reset_index(drop=True)

        if return_clean and return_matches:
            return df_short_clean, df_short_overlap, matches

        if return_clean:
            return df_short_clean, df_short_overlap

        if return_matches:
            return df_short_overlap, matches

        return df_short_overlap

    # ------------------------------------------------------------
    # 8. Identify the short and long index columns
    # ------------------------------------------------------------
    # Most bioframe versions return index_short/index_long when
    # suffixes=("_short", "_long") is used.
    if "index_short" in matches.columns:
        short_index_col = "index_short"
    elif "index_1" in matches.columns:
        short_index_col = "index_1"
    elif "index" in matches.columns:
        short_index_col = "index"
    else:
        raise KeyError(
            "Could not identify the short-table index column in "
            "the bioframe overlap output. Available columns: "
            f"{matches.columns.tolist()}"
        )

    if "index_long" in matches.columns:
        long_index_col = "index_long"
    elif "index_2" in matches.columns:
        long_index_col = "index_2"
    else:
        long_index_col = None

    # ------------------------------------------------------------
    # 9. Calculate overlap length in base pairs
    # ------------------------------------------------------------
    # Different bioframe versions may use different overlap
    # coordinate names, so calculate the intersection directly
    # when necessary.
    if (
        "overlap_start" in matches.columns
        and "overlap_end" in matches.columns
    ):
        matches["overlap_BP"] = (
            matches["overlap_end"]
            - matches["overlap_start"]
        )

    elif (
        "overlap_StartBP" in matches.columns
        and "overlap_EndBP" in matches.columns
    ):
        matches["overlap_BP"] = (
            matches["overlap_EndBP"]
            - matches["overlap_StartBP"]
        )

    else:
        required_match_columns = {
            "StartBP_short",
            "EndBP_short",
            "StartBP_long",
            "EndBP_long",
        }

        missing_match_columns = (
            required_match_columns.difference(matches.columns)
        )

        if missing_match_columns:
            raise KeyError(
                "Could not calculate BP overlap length. "
                "Missing columns: "
                f"{sorted(missing_match_columns)}. "
                "Available columns: "
                f"{matches.columns.tolist()}"
            )

        matches["overlap_BP"] = (
            np.minimum(
                matches["EndBP_short"],
                matches["EndBP_long"],
            )
            - np.maximum(
                matches["StartBP_short"],
                matches["StartBP_long"],
            )
        )

    matches["overlap_BP"] = (
        matches["overlap_BP"].astype(np.int64)
    )

    # ------------------------------------------------------------
    # 10. Calculate short-segment overlap fraction
    # ------------------------------------------------------------
    if (
        "StartBP_short" in matches.columns
        and "EndBP_short" in matches.columns
    ):
        short_length_bp = (
            matches["EndBP_short"]
            - matches["StartBP_short"]
        )

    else:
        short_indices = (
            matches[short_index_col]
            .astype(np.int64)
            .to_numpy()
        )

        short_length_bp = (
            short.loc[short_indices, "EndBP"].to_numpy()
            - short.loc[short_indices, "StartBP"].to_numpy()
        )

    matches["short_length_BP"] = short_length_bp

    matches["short_overlap_fraction"] = np.where(
        short_length_bp > 0,
        matches["overlap_BP"] / short_length_bp,
        np.nan,
    )

    # ------------------------------------------------------------
    # 11. Determine whether short IBD is fully inside long IBD
    # ------------------------------------------------------------
    containment_columns = {
        "StartBP_short",
        "EndBP_short",
        "StartBP_long",
        "EndBP_long",
    }

    if containment_columns.issubset(matches.columns):
        matches["short_contained_in_long"] = (
            (
                matches["StartBP_long"]
                <= matches["StartBP_short"]
            )
            & (
                matches["EndBP_long"]
                >= matches["EndBP_short"]
            )
        )

    # ------------------------------------------------------------
    # 12. Obtain unique IDs of overlapping short IBD segments
    # ------------------------------------------------------------
    overlapping_short_ids = (
        matches[short_index_col]
        .dropna()
        .astype(np.int64)
        .unique()
    )

    # Sort IDs so output remains in the original short-table order.
    overlapping_short_ids = np.sort(
        overlapping_short_ids
    )

    # ------------------------------------------------------------
    # 13. Extract overlapping short IBD segments
    # ------------------------------------------------------------
    original_short = df_short.reset_index(drop=True)

    df_short_overlap = (
        original_short
        .iloc[overlapping_short_ids]
        .copy()
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------
    # 14. Extract short IBD segments remaining after removal
    # ------------------------------------------------------------
    keep_mask = np.ones(
        len(original_short),
        dtype=bool,
    )

    keep_mask[overlapping_short_ids] = False

    df_short_clean = (
        original_short
        .loc[keep_mask]
        .copy()
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------
    # 15. Internal validation
    # ------------------------------------------------------------
    if (
        len(df_short_clean)
        + len(df_short_overlap)
        != len(original_short)
    ):
        raise RuntimeError(
            "Internal consistency check failed: "
            "clean short IBD count + overlapping short IBD count "
            "does not equal the original short IBD count."
        )

    if matches["overlap_BP"].le(0).any():
        raise RuntimeError(
            "Bioframe returned one or more matches with "
            "non-positive BP overlap length."
        )

    # ------------------------------------------------------------
    # 16. Return requested tables
    # ------------------------------------------------------------
    if return_clean and return_matches:
        return (
            df_short_clean,
            df_short_overlap,
            matches,
        )

    if return_clean:
        return (
            df_short_clean,
            df_short_overlap,
        )

    if return_matches:
        return (
            df_short_overlap,
            matches,
        )

    return df_short_overlap

def summarize_filtered_ibd_vectorized(
    df_ibd,
    min_cms=(2, 4, 6, 8),
    normalize_pair_order=False,
    output=False,
    ind_savepath="",
):
    """
    Create a pair-level IBD summary from a filtered segment-level
    IBD dataframe using vectorized pandas groupby.

    Parameters
    ----------
    df_ibd : pandas.DataFrame
        Filtered segment-level IBD table.

        Required columns:
            iid1, iid2, lengthM

    min_cms : iterable
        IBD length cutoffs in cM.

        For every cutoff, the function calculates:
            sum_IBD>{cutoff}
            n_IBD>{cutoff}

    normalize_pair_order : bool, default=False
        If True, treat A-B and B-A as the same pair before grouping.

    output : bool
        Print summary information.

    ind_savepath : str
        Optional path for saving the pair-level summary table.
        No file is written when this is an empty string.

    Returns
    -------
    pandas.DataFrame
        One row per individual pair.
    """

    # ------------------------------------------------------------
    # 1. Validate input
    # ------------------------------------------------------------
    required = {
        "iid1",
        "iid2",
        "lengthM",
    }

    missing = required.difference(
        df_ibd.columns
    )

    if missing:
        raise ValueError(
            "IBD segment table is missing required columns: "
            f"{sorted(missing)}"
        )

    min_cms = list(
        min_cms
    )

    if not min_cms:
        raise ValueError(
            "min_cms must contain at least one cutoff."
        )

    pair_cols = [
        "iid1",
        "iid2",
    ]

    output_cols = (
        pair_cols
        + ["max_IBD"]
    )

    for cutoff in min_cms:
        output_cols.extend(
            [
                f"sum_IBD>{cutoff}",
                f"n_IBD>{cutoff}",
            ]
        )

    # ------------------------------------------------------------
    # 2. Handle an empty segment table
    # ------------------------------------------------------------
    if df_ibd.empty:
        df_ind = pd.DataFrame(
            columns=output_cols
        )

        if ind_savepath:
            parent = os.path.dirname(
                ind_savepath
            )

            if parent:
                os.makedirs(
                    parent,
                    exist_ok=True,
                )

            df_ind.to_csv(
                ind_savepath,
                sep="\t",
                index=False,
            )

        return df_ind

    # ------------------------------------------------------------
    # 3. Retain columns required for summary
    # ------------------------------------------------------------
    work = df_ibd[
        [
            "iid1",
            "iid2",
            "lengthM",
        ]
    ].copy()

    # Ensure lengthM is numeric.
    work["lengthM"] = pd.to_numeric(
        work["lengthM"],
        errors="raise",
    )

    if work["lengthM"].isna().any():
        raise ValueError(
            "The IBD segment table contains missing lengthM values."
        )

    if (work["lengthM"] < 0).any():
        raise ValueError(
            "The IBD segment table contains negative lengthM values."
        )

    # ------------------------------------------------------------
    # 4. Optionally normalize pair orientation
    # ------------------------------------------------------------
    if normalize_pair_order:
        iid1 = (
            work["iid1"]
            .astype(str)
            .to_numpy()
        )

        iid2 = (
            work["iid2"]
            .astype(str)
            .to_numpy()
        )

        work["iid1"] = np.minimum(
            iid1,
            iid2,
        )

        work["iid2"] = np.maximum(
            iid1,
            iid2,
        )

    # ------------------------------------------------------------
    # 5. Convert Morgan to cM once
    # ------------------------------------------------------------
    work["_length_cm"] = (
        work["lengthM"]
        .to_numpy(dtype=np.float64)
        * 100.0
    )

    # ------------------------------------------------------------
    # 6. Build vectorized aggregation columns
    # ------------------------------------------------------------
    aggregation = {
        "max_IBD": (
            "_length_cm",
            "max",
        ),
    }

    for cutoff in min_cms:
        sum_col = (
            f"_sum_gt_{cutoff}"
        )

        n_col = (
            f"_n_gt_{cutoff}"
        )

        above_cutoff = (
            work["_length_cm"] > cutoff
        )

        work[sum_col] = (
            work["_length_cm"]
            .where(
                above_cutoff,
                0.0,
            )
        )

        work[n_col] = (
            above_cutoff
            .astype(np.int8)
        )

        aggregation[
            f"sum_IBD>{cutoff}"
        ] = (
            sum_col,
            "sum",
        )

        aggregation[
            f"n_IBD>{cutoff}"
        ] = (
            n_col,
            "sum",
        )

    # ------------------------------------------------------------
    # 7. Generate pair-level summary
    # ------------------------------------------------------------
    df_ind = (
        work
        .groupby(
            pair_cols,
            as_index=False,
            sort=False,
            observed=True,
        )
        .agg(
            **aggregation
        )
    )

    # Count columns should be integers.
    for cutoff in min_cms:
        df_ind[
            f"n_IBD>{cutoff}"
        ] = (
            df_ind[
                f"n_IBD>{cutoff}"
            ]
            .astype(np.int64)
        )

    # Sort by the first requested cutoff.
    df_ind = (
        df_ind
        .sort_values(
            by=f"sum_IBD>{min_cms[0]}",
            ascending=False,
            ignore_index=True,
        )
    )

    # ------------------------------------------------------------
    # 8. Optionally save pair-level summary
    # ------------------------------------------------------------
    if ind_savepath:
        parent = os.path.dirname(
            ind_savepath
        )

        if parent:
            os.makedirs(
                parent,
                exist_ok=True,
            )

        df_ind.to_csv(
            ind_savepath,
            sep="\t",
            index=False,
        )

        if output:
            print(
                f"Saved individual-pair summary to:\n"
                f"{ind_savepath}"
            )

    if output:
        print(
            f"Generated summary for "
            f"{len(df_ind):,} individual pairs "
            f"from {len(df_ibd):,} IBD segments."
        )

    return df_ind

def process_one_batch_pair(
    b1,
    b2,
    folder_short,
    folder_long,
    output_root,
    chs=range(1, 23),
    short_min_cms=(2, 4, 6, 8),
    long_min_cms=(8, 12, 16, 20),
    snp_cm=220,
    short_min_cm=2,
    long_min_cm=8,
    normalize_pair_order=True,
    overwrite=False,
    verbose=True,
):
    """
    Process one short/long IBD batch pair and save all intermediate
    and final tables.

    Workflow
    --------
    1. Read and filter the original short-IBD segment table.
    2. Read and filter the original long-IBD segment table.
    3. Save both original filtered segment tables.
    4. Find short IBD segments overlapping long IBD.
    5. Remove overlapping short IBD segments.
    6. Summarize the cleaned short-IBD table.
    7. Save all output tables.

    Output directory
    ----------------
    output_root/batchX_Y/

    Output files
    ------------
    df_ibd_short_original.tsv
        Original filtered short-IBD segment table before removing
        overlap with long IBD.

    df_ibd_long_original.tsv
        Original filtered long-IBD segment table used as the
        reference for overlap filtering.

    df_short_clean.tsv
        Short-IBD segments remaining after overlap removal.

    df_short_overlap.tsv
        Unique short-IBD segments removed because they overlap
        at least one long-IBD segment.

    df_matches.tsv
        Complete short-long overlap relationships.

    df_ind_short_clean.tsv
        Pair-level summary calculated from df_short_clean.

    batch_summary.tsv
        Batch-level segment counts and runtime.

    _SUCCESS
        Marker indicating successful completion.

    _ERROR.txt
        Error details, generated only if the batch fails.
    """

    b1 = int(b1)
    b2 = int(b2)

    batch_name = f"batch{b1}_{b2}"
    batch_output_dir = Path(output_root) / batch_name

    batch_output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    success_file = batch_output_dir / "_SUCCESS"
    error_file = batch_output_dir / "_ERROR.txt"

    # ------------------------------------------------------------
    # 1. Skip completed batches unless overwrite=True
    # ------------------------------------------------------------
    if success_file.exists() and not overwrite:
        if verbose:
            print(
                f"[SKIP] {batch_name}: already completed."
            )

        return {
            "batch": batch_name,
            "b1": b1,
            "b2": b2,
            "status": "skipped",
            "runtime_seconds": 0.0,
        }

    # Remove old completion/error markers before rerunning.
    if success_file.exists():
        success_file.unlink()

    if error_file.exists():
        error_file.unlink()

    start_time = time.perf_counter()

    current_batch = np.array(
        [[b1, b2]],
        dtype=int,
    )

    if verbose:
        print()
        print("=" * 72)
        print(f"Processing {batch_name}")
        print("=" * 72)

    try:
        # --------------------------------------------------------
        # 2. Read and filter original short-IBD segments
        # --------------------------------------------------------
        if verbose:
            print("[1/6] Reading original short-IBD segments")

        df_ibd_short = to_ibd_df_batches(
            batches=current_batch,
            chs=chs,
            folder_out=folder_short,
            min_cms=list(short_min_cms),
            snp_cm=snp_cm,
            min_cm=short_min_cm,
            output=False,
            savepath="",
        )

        if df_ibd_short is None:
            df_ibd_short = pd.DataFrame()

        df_ibd_short = (
            df_ibd_short
            .reset_index(drop=True)
        )

        n_short_original = len(df_ibd_short)

        if verbose:
            print(
                f"      Original short IBD: "
                f"{n_short_original:,}"
            )

        # --------------------------------------------------------
        # 3. Read and filter original long-IBD segments
        # --------------------------------------------------------
        if verbose:
            print("[2/6] Reading original long-IBD segments")

        df_ibd_long = to_ibd_df_batches(
            batches=current_batch,
            chs=chs,
            folder_out=folder_long,
            min_cms=list(long_min_cms),
            snp_cm=snp_cm,
            min_cm=long_min_cm,
            output=False,
            savepath="",
        )

        if df_ibd_long is None:
            df_ibd_long = pd.DataFrame()

        df_ibd_long = (
            df_ibd_long
            .reset_index(drop=True)
        )

        n_long_original = len(df_ibd_long)

        if verbose:
            print(
                f"      Original long IBD: "
                f"{n_long_original:,}"
            )

        # --------------------------------------------------------
        # 4. Save original short and long segment tables
        # --------------------------------------------------------
        if verbose:
            print("[3/6] Saving original IBD segment tables")

        short_original_path = (
            batch_output_dir
            / "df_ibd_short_original.tsv"
        )

        long_original_path = (
            batch_output_dir
            / "df_ibd_long_original.tsv"
        )

        df_ibd_short.to_csv(
            short_original_path,
            sep="\t",
            index=False,
        )

        df_ibd_long.to_csv(
            long_original_path,
            sep="\t",
            index=False,
        )

        # --------------------------------------------------------
        # 5. Find and remove short-long overlap
        # --------------------------------------------------------
        if verbose:
            print("[4/6] Identifying short-long IBD overlaps")

        if df_ibd_short.empty:
            # No short IBD: all output tables are empty.
            df_short_clean = (
                df_ibd_short.copy()
            )

            df_short_overlap = (
                df_ibd_short.copy()
            )

            df_matches = pd.DataFrame()

        elif df_ibd_long.empty:
            # No long IBD: retain all short IBD segments.
            df_short_clean = (
                df_ibd_short
                .copy()
                .reset_index(drop=True)
            )

            df_short_overlap = (
                df_ibd_short
                .iloc[0:0]
                .copy()
            )

            df_matches = pd.DataFrame()

        else:
            (
                df_short_clean,
                df_short_overlap,
                df_matches,
            ) = find_short_ibd_overlapping_long(
                df_short=df_ibd_short,
                df_long=df_ibd_long,
                normalize_pair_order=(
                    normalize_pair_order
                ),
                return_clean=True,
                return_matches=True,
            )

        n_short_clean = len(df_short_clean)
        n_short_overlap = len(df_short_overlap)
        n_matches = len(df_matches)

        # Each original short segment must be either retained
        # or removed.
        if (
            n_short_clean
            + n_short_overlap
            != n_short_original
        ):
            raise RuntimeError(
                f"{batch_name}: short-IBD count check failed. "
                f"Original={n_short_original:,}; "
                f"clean={n_short_clean:,}; "
                f"overlap={n_short_overlap:,}."
            )

        if verbose:
            print(
                f"      Clean short IBD:   "
                f"{n_short_clean:,}"
            )
            print(
                f"      Removed short IBD: "
                f"{n_short_overlap:,}"
            )
            print(
                f"      Overlap matches:   "
                f"{n_matches:,}"
            )

        # --------------------------------------------------------
        # 6. Summarize cleaned short-IBD segments
        # --------------------------------------------------------
        if verbose:
            print("[5/6] Summarizing cleaned short IBD")

        if df_short_clean.empty:
            summary_columns = [
                "iid1",
                "iid2",
                "max_IBD",
            ]

            for cutoff in short_min_cms:
                summary_columns.extend(
                    [
                        f"sum_IBD>{cutoff}",
                        f"n_IBD>{cutoff}",
                    ]
                )

            df_ind_short_clean = pd.DataFrame(
                columns=summary_columns
            )

        else:
            df_ind_short_clean = (
                summarize_filtered_ibd_vectorized(
                    df_ibd=df_short_clean,
                    min_cms=short_min_cms,
                    normalize_pair_order=(
                        normalize_pair_order
                    ),
                    output=False,
                )
            )

        n_clean_pairs = len(df_ind_short_clean)

        # --------------------------------------------------------
        # 7. Save cleaned, removed, matched and summary tables
        # --------------------------------------------------------
        if verbose:
            print("[6/6] Saving final output tables")

        df_short_clean.to_csv(
            batch_output_dir
            / "df_short_clean.tsv",
            sep="\t",
            index=False,
        )

        df_short_overlap.to_csv(
            batch_output_dir
            / "df_short_overlap.tsv",
            sep="\t",
            index=False,
        )

        df_matches.to_csv(
            batch_output_dir
            / "df_matches.tsv",
            sep="\t",
            index=False,
        )

        df_ind_short_clean.to_csv(
            batch_output_dir
            / "df_ind_short_clean.tsv",
            sep="\t",
            index=False,
        )

        runtime_seconds = (
            time.perf_counter()
            - start_time
        )

        removed_fraction = (
            n_short_overlap / n_short_original
            if n_short_original > 0
            else np.nan
        )

        batch_summary = pd.DataFrame(
            [
                {
                    "batch": batch_name,
                    "b1": b1,
                    "b2": b2,
                    "n_short_original": (
                        n_short_original
                    ),
                    "n_long_original": (
                        n_long_original
                    ),
                    "n_short_clean": (
                        n_short_clean
                    ),
                    "n_short_overlap": (
                        n_short_overlap
                    ),
                    "short_removed_fraction": (
                        removed_fraction
                    ),
                    "n_overlap_matches": (
                        n_matches
                    ),
                    "n_clean_pairs": (
                        n_clean_pairs
                    ),
                    "runtime_seconds": (
                        runtime_seconds
                    ),
                    "status": "completed",
                }
            ]
        )

        batch_summary.to_csv(
            batch_output_dir
            / "batch_summary.tsv",
            sep="\t",
            index=False,
        )

        # Write success marker only after all output files
        # have been successfully saved.
        success_file.write_text(
            (
                f"batch\t{batch_name}\n"
                f"status\tcompleted\n"
                f"runtime_seconds\t"
                f"{runtime_seconds:.2f}\n"
                f"n_short_original\t"
                f"{n_short_original}\n"
                f"n_long_original\t"
                f"{n_long_original}\n"
                f"n_short_clean\t"
                f"{n_short_clean}\n"
                f"n_short_overlap\t"
                f"{n_short_overlap}\n"
            )
        )

        if verbose:
            print(
                f"[DONE] {batch_name}: "
                f"{runtime_seconds:.2f} seconds"
            )

        result = (
            batch_summary
            .iloc[0]
            .to_dict()
        )

        # Original tables have already been written to disk,
        # so they can now be safely removed from memory.
        del df_ibd_short
        del df_ibd_long
        del df_short_clean
        del df_short_overlap
        del df_matches
        del df_ind_short_clean
        del batch_summary

        gc.collect()

        return result

    except Exception as error:
        runtime_seconds = (
            time.perf_counter()
            - start_time
        )

        error_message = (
            f"Batch: {batch_name}\n"
            f"Runtime: {runtime_seconds:.2f} seconds\n"
            f"Error: {repr(error)}\n\n"
            f"{traceback.format_exc()}"
        )

        error_file.write_text(
            error_message
        )

        print(
            f"[ERROR] {batch_name}: "
            f"{repr(error)}"
        )

        # Remove any dataframes that may exist before continuing.
        for variable_name in [
            "df_ibd_short",
            "df_ibd_long",
            "df_short_clean",
            "df_short_overlap",
            "df_matches",
            "df_ind_short_clean",
            "batch_summary",
        ]:
            if variable_name in locals():
                del locals()[variable_name]

        gc.collect()

        return {
            "batch": batch_name,
            "b1": b1,
            "b2": b2,
            "status": "failed",
            "runtime_seconds": runtime_seconds,
            "error": repr(error),
        }
