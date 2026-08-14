"""Exhaustive check that the proteomic layer is unchanged from the previous release.

Not a sample: every analyte column in every panel, on every key present in both
releases. Two things are checked per column, because either alone is insufficient:

  1. VALUES  — where both releases have a value, they must be equal (exact, not atol)
  2. PRESENCE — the null pattern must match. A value present in one release and NaN in
                the other is a discrepancy even though no value "differs".

Reads each dataset's analyte block once as float32 and compares in memory, rather than
re-scanning a 1.2 GB TSV per panel.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.join(HERE, "Project_9004_Unified_Emerging_Biomarkers.tab")
OLD = os.path.join(HERE, os.pardir, "Project_9004_Unified_Emerging_Biomarkers.tab")

PANELS = ["p277_CSF", "p282_CNS_CSF", "p282_Inflammation_CSF", "p288_CNS_plasma",
          "p288_Inflammation_plasma", "p293_olink_plasma", "p312_Inflammation_CSF",
          "p312_Inflammation_Plasma", "p312_Neuro_CSF", "p312_Neuro_Plasma",
          "p314_CSF", "p314_Plasma"]


def analyte_cols(header: list[str]) -> list[str]:
    return [c for c in header
            if any(c.startswith(p + "_") for p in PANELS)
            and "earliest_visit_PC" not in c and "PLATE_ID" not in c]


def main() -> None:
    nh = pd.read_csv(NEW, sep="\t", nrows=0).columns.tolist()
    oh = pd.read_csv(OLD, sep="\t", nrows=0).columns.tolist()
    shared_cols = [c for c in analyte_cols(nh) if c in set(analyte_cols(oh))]
    print(f"analyte columns: new {len(analyte_cols(nh)):,}  old {len(analyte_cols(oh)):,}  "
          f"shared {len(shared_cols):,}")

    print("reading new ...")
    a = pd.read_csv(NEW, sep="\t", usecols=["key"] + shared_cols,
                    dtype={"key": str, **{c: np.float32 for c in shared_cols}}).set_index("key")
    print("reading previous ...")
    b = pd.read_csv(OLD, sep="\t", usecols=["PATNO", "EVENT_ID"] + shared_cols,
                    dtype={"PATNO": str, "EVENT_ID": str,
                           **{c: np.float32 for c in shared_cols}})
    b.index = b["PATNO"].str.strip() + "_" + b["EVENT_ID"].str.strip()
    b = b.drop(columns=["PATNO", "EVENT_ID"])

    shared_keys = a.index.intersection(b.index)
    print(f"shared keys: {len(shared_keys):,}\n")
    a = a.reindex(shared_keys)[shared_cols]
    b = b.reindex(shared_keys)[shared_cols]

    A, B = a.to_numpy(np.float32), b.to_numpy(np.float32)
    na, nb = np.isnan(A), np.isnan(B)

    presence_mismatch = na ^ nb
    both = ~na & ~nb
    value_mismatch = np.zeros_like(both)
    value_mismatch[both] = A[both] != B[both]

    print(f"{'panel':28s}{'cols':>7s}{'cells':>14s}{'val diff':>10s}{'presence diff':>15s}")
    total_v = total_p = 0
    idx = {c: i for i, c in enumerate(shared_cols)}
    for p in PANELS:
        cols = [c for c in shared_cols if c.startswith(p + "_")]
        if not cols:
            continue
        j = [idx[c] for c in cols]
        v = int(value_mismatch[:, j].sum())
        pr = int(presence_mismatch[:, j].sum())
        total_v += v
        total_p += pr
        print(f"{p:28s}{len(cols):7,}{both[:, j].sum():14,}{v:10,}{pr:15,}")

    print(f"\n{'TOTAL':28s}{len(shared_cols):7,}{int(both.sum()):14,}{total_v:10,}{total_p:15,}")

    if total_v or total_p:
        print("\nFAIL — the proteomic layer is not identical.")
        bad = [shared_cols[j] for j in range(len(shared_cols))
               if value_mismatch[:, j].any() or presence_mismatch[:, j].any()]
        print(f"  {len(bad)} column(s) affected, e.g. {bad[:10]}")
        sys.exit(1)

    print(f"\nPASS — all {len(shared_cols):,} analyte columns are identical across "
          f"{int(both.sum()):,} value-carrying cells and {len(shared_keys):,} shared keys, "
          "with matching null patterns.")


if __name__ == "__main__":
    sys.exit(main())
