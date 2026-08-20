#!/usr/bin/env python3
"""Analyte concordance on overlapping samples, for the AD/PD abstract.

Three comparisons, all restricted to participant-VISITS where both members of the pair
were actually measured (the overlapping sample), because a concordance computed on a
union of two disjoint cohorts is not a concordance at all:

  A. SAME ASSAY, CSF vs PLASMA     -- does the biofluid carry the same information?
                                      Olink CSF vs Olink plasma; NULISA CSF vs plasma.
  B. BETWEEN ASSAYS, WITHIN CSF    -- do two technologies agree on the same protein
  C. BETWEEN ASSAYS, WITHIN PLASMA    in the same fluid? NULISA (NPQ) vs Olink (NPX).

R-squared is the headline: it is the shared-variance number a general audience reads
directly as "how much of one measurement the other explains", and it is invariant to
the NPQ/NPX unit difference the way a raw slope would not be. Spearman is carried
alongside as a monotone-agreement check, since a single influential batch can prop up
Pearson on assay data.

Harmonized columns only, so each side is already pooled across its project pair.
"""
from __future__ import annotations

import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CWD = os.path.dirname(os.path.abspath(__file__))
INPUT_GLOB = "Project_9004_Unified_Emerging_Biomarkers_*.tab"
BLOCKS = ["olink_plasma", "olink_csf", "nulisa_cns_plasma", "nulisa_cns_csf",
          "nulisa_inf_plasma", "nulisa_inf_csf"]

# Same floor the pipeline applies to any reportable cell.
MIN_PAIRS = 20

# (label, block_a, block_b) -- the three comparisons above, NULISA split by panel.
COMPARISONS = [
    ("A. Same assay: CSF vs plasma (Olink)",   "olink_csf",         "olink_plasma"),
    ("A. Same assay: CSF vs plasma (NULISA)",  "nulisa_cns_csf",    "nulisa_cns_plasma"),
    ("A. Same assay: CSF vs plasma (NULISA)",  "nulisa_inf_csf",    "nulisa_inf_plasma"),
    ("B. Between assays, within CSF",          "nulisa_cns_csf",    "olink_csf"),
    ("B. Between assays, within CSF",          "nulisa_inf_csf",    "olink_csf"),
    ("C. Between assays, within plasma",       "nulisa_cns_plasma", "olink_plasma"),
    ("C. Between assays, within plasma",       "nulisa_inf_plasma", "olink_plasma"),
]


def main() -> None:
    import regressions as R

    path = R.resolve_input_path(INPUT_GLOB, CWD)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")

    pat = re.compile(r"^harmonized_(" + "|".join(BLOCKS) + r")_(.+)_(NPX|NPQ)$")
    by_block: dict[str, dict[str, str]] = defaultdict(dict)
    for c in header:
        m = pat.match(c)
        if m:
            by_block[m.group(1)][m.group(2)] = c

    shared = [(label, a, b, by_block[a][ana], by_block[b][ana], ana)
              for label, a, b in COMPARISONS
              for ana in sorted(set(by_block[a]) & set(by_block[b]))]

    # PATNO + EVENT_ID come along so the pairing unit is DEMONSTRABLY one participant at
    # one timepoint, not merely implied by the row. A row of this release is exactly one
    # participant-visit, so requiring both analytes non-null on the same row already
    # matches same-person/same-visit -- but the abstract claims it, so it gets counted
    # and reported rather than asserted.
    KEYS = ["PATNO", "EVENT_ID"]
    need = sorted({c for row in shared for c in row[3:5]} | set(KEYS))
    print(f"{len(shared):,} analyte pairs; reading {len(need):,} of {len(header):,} columns ...",
          flush=True)
    df = pd.read_csv(path, sep="\t", usecols=need, low_memory=False)
    print(f"  {df.shape[0]:,} rows x {df.shape[1]:,} cols", flush=True)
    dup = int(df.duplicated(KEYS).sum())
    print(f"  rows per (PATNO, EVENT_ID): {'1:1' if dup == 0 else f'NOT UNIQUE ({dup} dups)'}\n",
          flush=True)

    rows = []
    for label, a, b, ca, cb, ana in shared:
        s = df[KEYS + [ca, cb]].dropna(subset=[ca, cb])   # <- the overlapping sample
        if len(s) < MIN_PAIRS:
            continue
        x, y = s[ca].to_numpy(float), s[cb].to_numpy(float)
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        r, p = pearsonr(x, y)
        rho, _ = spearmanr(x, y)
        if not np.isfinite(r):
            continue
        rows.append({"comparison": label, "analyte": ana, "block_a": a, "block_b": b,
                     "n_overlap": len(s), "n_participants": int(s.PATNO.nunique()),
                     "n_visits": int(s.EVENT_ID.nunique()),
                     "r": float(r), "r2": float(r) ** 2,
                     "rho": float(rho), "P": float(p)})

    res = pd.DataFrame(rows)

    # An R-squared whose correlation is not significant is noise around zero, and
    # pooling those into a median describes the null rather than the agreement. A pair
    # at P > 0.05 is therefore reported as NS and carries no number at all.
    #
    # NOTE THE SAMPLE SIZE when reading the nominal column. These overlaps run 2,000-3,600
    # participant-visits, and at n=3,000 the nominal 0.05 threshold is cleared by
    # r ~ 0.036, i.e. R2 ~ 0.001 -- a correlation that is "significant" and still
    # explains a tenth of a percent of the variance. Bonferroni across the analytes in
    # each comparison is carried alongside for that reason; where the two disagree it is
    # the nominal column that is admitting noise, not the corrected one rejecting signal.
    #
    # Both the RATE and the surviving R-squared are reported, because either alone
    # misleads: the survivor median is upward-biased by the selection that produced it,
    # while the rate is what says how much of the panel agrees at all.
    # NS is the VERDICT, not an erasure: the r2 column keeps its value on every row so
    # nothing is lost from the record, and `verdict` is what a table or slide prints.
    res["ns"] = res["P"] > 0.05
    res["alpha_bonf"] = res.groupby("comparison")["P"].transform(lambda p: 0.05 / len(p))
    res["sig_bonf"] = res["P"] < res["alpha_bonf"]
    res["verdict"] = np.where(res["ns"], "NS", res["r2"].round(3).astype(str))

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "concordance.csv")
    res.sort_values(["comparison", "r2"], ascending=[True, False]).to_csv(out, index=False)

    print("=" * 78)
    for label, g in res.groupby("comparison"):
        sig, bonf, ns = g[~g.ns], g[g.sig_bonf], g[g.ns]
        print(f"\n{label}")
        print(f"   {len(g):>5,} analytes tested   median overlap {int(g.n_overlap.median()):,} "
              f"participant-visits, {int(g.n_participants.median()):,} participants across "
              f"{int(g.n_visits.median())} timepoints")
        print(f"   NS (P>0.05): {len(ns):,}/{len(g):,} ({100*len(ns)/len(g):.1f}%)"
              f"   —   note their R2 anyway: median {ns.r2.median():.4f}"
              if len(ns) else f"   NS (P>0.05): 0/{len(g):,}")
        if sig.empty:
            print("   nothing significant — no R2 to report")
            continue
        print(f"   P<0.05      : {len(sig):,} analytes   median R2 {sig.r2.median():.3f}   "
              f"IQR {sig.r2.quantile(.25):.3f}-{sig.r2.quantile(.75):.3f}")
        print(f"   Bonferroni  : {len(bonf):,} analytes   "
              + (f"median R2 {bonf.r2.median():.3f}   "
                 f"IQR {bonf.r2.quantile(.25):.3f}-{bonf.r2.quantile(.75):.3f}"
                 if not bonf.empty else "none survive")
              + f"   (α={g.alpha_bonf.iloc[0]:.2e})")
        top = sig.nlargest(min(6, len(sig)), "r2")
        print("   strongest: " + ", ".join(f"{r.analyte} {r.r2:.2f}" for r in top.itertuples()))

    # Roll the two NULISA panels together so the abstract can quote one number per row.
    print("\n" + "=" * 78)
    print("\nHEADLINE (panels pooled) — median R2 among analytes that reach each threshold:")
    print(f"   {'comparison':32s} {'tested':>7s} {'NS':>7s} {'P<0.05':>16s} {'Bonferroni':>16s}")
    for key, pref in [("Same assay, CSF vs plasma", "A."),
                      ("Between assays, within CSF", "B."),
                      ("Between assays, within plasma", "C.")]:
        g = res[res.comparison.str.startswith(pref)]
        sig, bonf = g[~g.ns], g[g.sig_bonf]
        f_sig = f"{len(sig):,} @ {sig.r2.median():.3f}" if not sig.empty else "NS"
        f_bon = f"{len(bonf):,} @ {bonf.r2.median():.3f}" if not bonf.empty else "NS"
        print(f"   {key:32s} {len(g):>7,} {len(g[g.ns]):>7,} {f_sig:>16s} {f_bon:>16s}")

    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
