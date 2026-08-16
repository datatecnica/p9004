#!/usr/bin/env python3
"""build_composition_table.py — Table 2: composition of the unified dataset.

Regenerates the release's composition table from the shipped dataset rather than from
the build logs, so the numbers are a property of the artifact a reader actually holds.

Per column group it reports:
  columns            how many columns the group contributes
  missingness        share of visit rows carrying no value anywhere in the group
                     (genetics is static per participant, so it is ALSO quoted as the
                     share of participants lacking the block — the row-level figure
                     alone would misdescribe a participant-level layer)
  participants       unique PATNOs with at least one value in the group
  visits             visit rows with at least one value in the group
  visits/participant mean (SD) and median, over the participants who have the group
                     at all. Both are reported because the distribution is skewed:
                     a block assayed once for most participants and repeatedly for a
                     longitudinal subset has a mean pulled above a median of 1, and
                     the mean alone would suggest longitudinal depth the block does
                     not have.

The last three are the columns the missingness percentage cannot carry: a block can be
"90% missing" because it was assayed once per participant or because it covers few
participants densely, and only the counts distinguish those.

Group membership is read from the build intermediates (scaffold, derived, PCs,
harmonized) rather than pattern-matched, so a column cannot land in two groups or be
silently dropped — the script asserts the groups partition the dataset exactly.

Emits:
  dataset_composition.md
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

from build_common import DATASET_STEM, require_build

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")
DATA = require_build(DATASET_STEM)
REPORT = os.path.join(HERE, "dataset_composition.md")

CHUNK = 1000  # rows per pass; 35.6k columns x 1k rows is a comfortable working set


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} in build_intermediates — run the build first")
    return hits[-1]


def header_of(path: str) -> list[str]:
    return pd.read_csv(path, sep="\t", nrows=0).columns.tolist()


def main() -> None:
    print(f"dataset: {os.path.basename(DATA)}")
    cols = header_of(DATA)
    print(f"columns: {len(cols):,}")

    # ---- group membership, from the intermediates that produced each block ----
    scaffold = set(header_of(latest("scaffold-*.tab")))
    derived = set(header_of(latest("derived-*.tab")))
    harmon = set(header_of(latest("harmonized-*.tab")))

    ID_COLS = {"key", "PATNO", "EVENT_ID", "MERGE_INDEX"}

    def group_of(c: str) -> str:
        if c in ID_COLS:
            return "Technical"
        # Per-block join artifacts: each merged block retains its own PATNO_<block> /
        # EVENT_ID_<block> so a bad join stays auditable after the fact. Named with a
        # suffix, so an exact-match ID test misses them.
        if c.startswith(("PATNO_", "EVENT_ID_")):
            return "Technical"
        if "earliest_visit_PC" in c:
            return "Derived + baseline-PC"
        if "PLATE_ID" in c or c.endswith("_src"):
            return "Technical"
        if c == "collection_era":
            return "Technical"
        if c.startswith("harmonized_"):
            return "Harmonized proteomic blocks (6)"
        for p, label in (("p277_", "Project 277 (Olink CSF)"),
                         ("p282_", "Project 282 (NULISA CSF)"),
                         ("p288_", "Project 288 (NULISA plasma)"),
                         ("p293_", "Project 293 (Olink plasma)"),
                         ("p312_", "Project 312 (NULISA CSF + plasma)"),
                         ("p314_", "Project 314 (Olink CSF + plasma)")):
            if c.startswith(p):
                return label
        if c.startswith("GP2_") or c.startswith("p9005_"):
            return "Genetics (GP2 + Project 9005)"
        if c in derived:
            return "Derived + baseline-PC"
        if c in scaffold:
            return "Clinical / demographic"
        if c in harmon:
            return "Technical"
        return "UNASSIGNED"

    groups: dict[str, list[str]] = {}
    for c in cols:
        groups.setdefault(group_of(c), []).append(c)

    if "UNASSIGNED" in groups:
        sys.exit(f"FATAL: {len(groups['UNASSIGNED'])} column(s) matched no group, "
                 f"e.g. {groups['UNASSIGNED'][:10]}")
    assert sum(len(v) for v in groups.values()) == len(cols), "groups do not partition"

    order = ["Clinical / demographic", "Genetics (GP2 + Project 9005)",
             "Project 282 (NULISA CSF)", "Project 288 (NULISA plasma)",
             "Project 277 (Olink CSF)", "Project 293 (Olink plasma)",
             "Project 312 (NULISA CSF + plasma)", "Project 314 (Olink CSF + plasma)",
             "Harmonized proteomic blocks (6)", "Derived + baseline-PC", "Technical"]
    order = [g for g in order if g in groups] + [g for g in groups if g not in order]

    # ---- one pass over the dataset, accumulating presence per group ----
    print(f"scanning in {CHUNK}-row chunks ...", flush=True)
    n_rows = 0
    patnos_all: set[str] = set()
    # per group: visit rows with data, and per-PATNO visit-with-data counts
    visits = {g: 0 for g in groups}
    per_patno = {g: {} for g in groups}

    reader = pd.read_csv(DATA, sep="\t", chunksize=CHUNK, low_memory=False,
                         dtype={"key": str, "PATNO": str, "EVENT_ID": str})
    for i, chunk in enumerate(reader):
        n_rows += len(chunk)
        pat = chunk["PATNO"].astype(str)
        patnos_all.update(pat)
        for g, gcols in groups.items():
            present = chunk[gcols].notna().any(axis=1)
            visits[g] += int(present.sum())
            if present.any():
                vc = pat[present].value_counts()
                d = per_patno[g]
                for p, n in vc.items():
                    d[p] = d.get(p, 0) + int(n)
        print(f"  chunk {i + 1}: {n_rows:,} rows", flush=True)

    n_pat = len(patnos_all)
    print(f"\n{n_rows:,} rows, {n_pat:,} participants\n")

    # ---- assemble ----
    rows = []
    for g in order:
        gcols = groups[g]
        counts = np.array(list(per_patno[g].values()), dtype=float)
        n_p = len(counts)
        mean = counts.mean() if n_p else 0.0
        sd = counts.std(ddof=1) if n_p > 1 else 0.0
        med = float(np.median(counts)) if n_p else 0.0
        vmin = int(counts.min()) if n_p else 0
        vmax = int(counts.max()) if n_p else 0
        row_miss = 100 * (1 - visits[g] / n_rows) if n_rows else 0.0
        pat_miss = 100 * (1 - n_p / n_pat) if n_pat else 0.0
        rows.append(dict(group=g, columns=len(gcols),
                         row_missing=row_miss, pat_missing=pat_miss,
                         participants=n_p, visits=visits[g],
                         mean_visits=mean, sd_visits=sd, median_visits=med,
                         min_visits=vmin, max_visits=vmax))
    tab = pd.DataFrame(rows)

    def fmt(r) -> str:
        miss = (f"{r.row_missing:.1f}%" if r.group != "Technical" else "—")
        return (f"| {r.group} | {r.columns:,} | {miss} | {r.pat_missing:.1f}% | "
                f"{r.participants:,} | {r.visits:,} | "
                f"{r.mean_visits:.2f} ({r.sd_visits:.2f}) | "
                f"{r.median_visits:.0f} | {r.min_visits}–{r.max_visits} |")

    lines = [
        "# Table 2 — composition of the unified dataset",
        "",
        f"`{os.path.basename(DATA)}` — **{n_rows:,} participant-visits × "
        f"{len(cols):,} columns**, **{n_pat:,} unique participants**.",
        "",
        "| Column group | Columns | Rows missing | Participants missing | "
        "Participants | Visits | Visits/participant, mean (SD) | Median | Range |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    lines += [fmt(r) for r in tab.itertuples()]
    lines.append(f"| **Total** | **{len(cols):,}** | — | — | **{n_pat:,}** | "
                 f"**{n_rows:,}** | — | — | — |")
    lines.append("")
    tab.to_csv(os.path.join(HERE, "dataset_composition.tab"), sep="\t", index=False)
    with open(REPORT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {os.path.basename(REPORT)} and dataset_composition.tab")


if __name__ == "__main__":
    main()
