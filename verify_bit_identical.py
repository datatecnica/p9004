"""Exhaustive check that the proteomic layer is unchanged from the previous release.

Not a sample: every analyte column in every panel, on every key present in both
releases. Three things are checked, because any one alone is insufficient:

  0. COLUMN SET — the two releases must carry the SAME analyte columns. Checks 1 and 2
                run on the intersection, so on their own they cannot tell "identical"
                from "identical wherever they happen to overlap": a release that
                dropped a whole panel would still pass on whatever remained. Fatal
                unless --allow-column-diff.
  1. VALUES  — where both releases have a value, they must be equal (exact, not atol)
  2. PRESENCE — the null pattern must match. A value present in one release and NaN in
                the other is a discrepancy even though no value "differs".

Key coverage is reported the same way (rows only in one release are named, not
silently skipped) but is not fatal, since releases legitimately gain participants.

Reads each dataset's analyte block once as float32 and compares in memory, rather than
re-scanning a 1.2 GB TSV per panel.

Usage:
    python3 verify_bit_identical.py                      # defaults below
    python3 verify_bit_identical.py --new A.tab --old B.tab
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

from build_common import BASELINE_DATASET_STEM, DATASET_STEM, find_build

HERE = os.path.dirname(os.path.abspath(__file__))
# find_build, not require_build: this is only the DEFAULT for --new, and resolving it
# strictly at import would exit before argparse could see an explicit --new/--old pair.
NEW = find_build(DATASET_STEM)
OLD = find_build(BASELINE_DATASET_STEM)

PANELS = ["p277_CSF", "p282_CNS_CSF", "p282_Inflammation_CSF", "p288_CNS_plasma",
          "p288_Inflammation_plasma", "p293_olink_plasma", "p312_Inflammation_CSF",
          "p312_Inflammation_Plasma", "p312_Neuro_CSF", "p312_Neuro_Plasma",
          "p314_CSF", "p314_Plasma"]


def analyte_cols(header: list[str]) -> list[str]:
    return [c for c in header
            if any(c.startswith(p + "_") for p in PANELS)
            and "earliest_visit_PC" not in c and "PLATE_ID" not in c]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Verify the proteomic layer is unchanged between two releases")
    ap.add_argument("--new", default=NEW, help="the release being checked")
    ap.add_argument("--old", default=OLD, help="the baseline release")
    ap.add_argument("--allow-column-diff", action="store_true",
                    help="report asymmetric analyte columns instead of failing on them")
    args = ap.parse_args()

    for label, path in (("new", args.new), ("old", args.old)):
        if path is None:
            sys.exit(f"ERROR: no --{label} given and none could be resolved "
                     f"(expected {DATASET_STEM}_<YYYYMMDD>.tab)")
        if not os.path.exists(path):
            sys.exit(f"ERROR: --{label} not found: {path}")
    print(f"new: {args.new}\nold: {args.old}\n")

    nh = pd.read_csv(args.new, sep="\t", nrows=0).columns.tolist()
    oh = pd.read_csv(args.old, sep="\t", nrows=0).columns.tolist()
    an, ao = analyte_cols(nh), analyte_cols(oh)
    set_n, set_o = set(an), set(ao)
    shared_cols = [c for c in an if c in set_o]

    # Comparing only the intersection cannot distinguish "the releases are identical"
    # from "the releases overlap on these columns and differ everywhere else" — a
    # release that silently dropped a panel would still PASS on whatever remained.
    # Raised by the stats core; identity now requires shared == both totals.
    only_new = sorted(set_n - set_o)
    only_old = sorted(set_o - set_n)
    print(f"analyte columns: new {len(an):,}  old {len(ao):,}  shared {len(shared_cols):,}")
    if only_new or only_old:
        print(f"  only in new: {len(only_new):,}" + (f"  e.g. {only_new[:5]}" if only_new else ""))
        print(f"  only in old: {len(only_old):,}" + (f"  e.g. {only_old[:5]}" if only_old else ""))
        if not args.allow_column_diff:
            print("\nFAIL — the two releases do not carry the same analyte columns, so "
                  "'identical' cannot be established on the intersection alone.\n"
                  "       Pass --allow-column-diff to compare the overlap anyway.")
            sys.exit(1)
    else:
        print(f"  column sets are identical ({len(shared_cols):,} = {len(an):,} = {len(ao):,})")

    print("\nreading new ...")
    a = pd.read_csv(args.new, sep="\t", usecols=["key"] + shared_cols,
                    dtype={"key": str, **{c: np.float32 for c in shared_cols}}).set_index("key")
    print("reading previous ...")
    b = pd.read_csv(args.old, sep="\t", usecols=["PATNO", "EVENT_ID"] + shared_cols,
                    dtype={"PATNO": str, "EVENT_ID": str,
                           **{c: np.float32 for c in shared_cols}})
    b.index = b["PATNO"].str.strip() + "_" + b["EVENT_ID"].str.strip()
    b = b.drop(columns=["PATNO", "EVENT_ID"])

    shared_keys = a.index.intersection(b.index)
    # Same blind spot one axis over: rows present in only one release are never compared.
    # Not fatal — releases legitimately gain participants — but it must be visible, or
    # "all shared keys match" reads as "every row matches".
    print(f"shared keys: {len(shared_keys):,} "
          f"(new {len(a.index):,}, old {len(b.index):,}; "
          f"only-new {len(a.index.difference(b.index)):,}, "
          f"only-old {len(b.index.difference(a.index)):,})\n")
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

    scope = ("all" if not (only_new or only_old) else "the shared")
    print(f"\nPASS — {scope} {len(shared_cols):,} analyte columns are identical across "
          f"{int(both.sum()):,} value-carrying cells and {len(shared_keys):,} shared keys, "
          "with matching null patterns.")
    if only_new or only_old:
        print(f"       NOTE: {len(only_new):,} column(s) exist only in new and "
              f"{len(only_old):,} only in old; those were not compared.")


if __name__ == "__main__":
    sys.exit(main())
