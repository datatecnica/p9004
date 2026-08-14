"""Step 4 — left-join every block onto the clinical scaffold.

The scaffold is the base and its row count is preserved exactly: 19,450 in, 19,450 out.
Proteomic blocks join on `key` (PATNO_EVENT_ID); genetics blocks join on `PATNO` and
propagate to all of a participant's visits.

Join artifacts (`PATNO_p312inflam`, `EVENT_ID_p314csf`, …) are retained and later
documented as `Technical` stubs, matching how Projects 282/277/312/314 were merged
previously.

Emits:
  build_intermediates/merged-<ts>.tab
  build_intermediates/merge_sparsity-<ts>.tab
"""

from __future__ import annotations

import glob
import os
import sys

import pandas as pd

from build_common import TIMESTAMP, assert_no_fatal, log, make_key, make_patno
from build_step_2_format_blocks import BLOCKS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} — run the prior step first")
    return hits[-1]


def main() -> None:
    log(f"=== Step 4: merge ===  timestamp {TIMESTAMP}")

    scaffold = pd.read_csv(latest("scaffold-*.tab"), sep="\t", low_memory=False,
                           dtype={"key": str, "MERGE_INDEX": str, "PATNO": str,
                                  "EVENT_ID": str})
    n_rows = len(scaffold)
    log(f"scaffold: {n_rows:,} rows x {scaffold.shape[1]} columns")

    key_index = pd.Index(scaffold["key"])
    assert key_index.is_unique, "scaffold key is not unique"

    pieces: list[pd.DataFrame] = [scaffold]
    sparsity = []

    # --- proteomic blocks --------------------------------------------------
    for b in BLOCKS:
        name, path, suffix = b["name"], b["file"], b["suffix"]
        df = pd.read_csv(path, sep="\t", low_memory=False,
                         dtype={"MERGE_INDEX": str, "PATNO": str, "EVENT_ID": str})

        has_ids = "PATNO" in df.columns and "EVENT_ID" in df.columns
        if has_ids:
            key, f = make_key(df, source=name)
        else:
            parts = df["MERGE_INDEX"].astype(str).str.rsplit("_", n=1, expand=True)
            key, f = make_key(pd.DataFrame({"PATNO": parts[0], "EVENT_ID": parts[1]}),
                              source=name)
        assert_no_fatal(f)

        df = df.drop(columns=["MERGE_INDEX"])
        if has_ids:
            # Retain as documented Technical stubs, under the previous build's names.
            df = df.rename(columns={"PATNO": f"PATNO_{suffix}",
                                    "EVENT_ID": f"EVENT_ID_{suffix}"})

        df.index = key
        # Keep only samples the scaffold can host; count the rest as unmatched.
        matched = int(key.isin(key_index).sum())
        df = df[~df.index.duplicated(keep="first")]
        df = df.reindex(key_index)
        df.index = scaffold.index

        analyte_cols = [c for c in df.columns
                        if not c.startswith(("PATNO_", "EVENT_ID_")) and "PLATE_ID" not in c]
        rows_with_data = int(df[analyte_cols].notna().any(axis=1).sum())
        pct_missing = float(df[analyte_cols].isna().sum().sum()) / (n_rows * len(analyte_cols)) * 100

        log(f"  {name:20s} +{df.shape[1]:5,} cols   matched {matched:5,}   "
            f"rows with data {rows_with_data:5,}   {pct_missing:5.2f}% missing")

        sparsity.append(dict(block=name, columns=df.shape[1], analyte_columns=len(analyte_cols),
                             vendor_samples=len(key), matched=matched,
                             rows_with_data=rows_with_data, pct_missing=round(pct_missing, 2)))
        pieces.append(df)

    # --- genetics blocks ---------------------------------------------------
    patno_index = pd.Index(scaffold["PATNO"])
    for label, pattern in (("GP2", "genetics_gp2-*.tab"), ("p9001", "genetics_p9001-*.tab")):
        g = pd.read_csv(latest(pattern), sep="\t", low_memory=False, dtype={"PATNO": str})
        p, f = make_patno(g, source=label)
        assert_no_fatal(f)
        g.index = p
        g = g.drop(columns=["PATNO"])
        g = g.reindex(patno_index)
        g.index = scaffold.index

        rows_with_data = int(g.notna().any(axis=1).sum())
        pct_missing = float(g.isna().sum().sum()) / (n_rows * g.shape[1]) * 100
        log(f"  {label:20s} +{g.shape[1]:5,} cols   "
            f"rows with data {rows_with_data:5,}   {pct_missing:5.2f}% missing")

        sparsity.append(dict(block=label, columns=g.shape[1], analyte_columns=g.shape[1],
                             vendor_samples=len(p), matched=int(p.isin(patno_index).sum()),
                             rows_with_data=rows_with_data, pct_missing=round(pct_missing, 2)))
        pieces.append(g)

    # --- assemble ----------------------------------------------------------
    log("\nconcatenating...")
    merged = pd.concat(pieces, axis=1)
    del pieces

    assert len(merged) == n_rows, f"row count changed: {n_rows:,} -> {len(merged):,}"
    dupe_cols = merged.columns[merged.columns.duplicated()].tolist()
    assert not dupe_cols, f"duplicate column names after merge: {dupe_cols[:10]}"
    log(f"merged: {merged.shape[0]:,} rows x {merged.shape[1]:,} columns")
    log("row count preserved; no duplicate column names")

    f_out = os.path.join(OUT, f"merged-{TIMESTAMP}.tab")
    log(f"\nwriting {os.path.basename(f_out)} ...")
    merged.to_csv(f_out, sep="\t", index=False)
    pd.DataFrame(sparsity).to_csv(
        os.path.join(OUT, f"merge_sparsity-{TIMESTAMP}.tab"), sep="\t", index=False)
    log(f"wrote {f_out}  ({os.path.getsize(f_out) / 1e9:.2f} GB)")


if __name__ == "__main__":
    sys.exit(main())
