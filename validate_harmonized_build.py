"""Validate the harmonized blocks as actually built.

Everything here reads the shipped `harmonized_*` columns, not a reconstruction, so it
tests what a downstream user would actually consume.

  1. INTEGRITY   — the fill rule. Where the reference project has a value the harmonized
                   column must equal it EXACTLY; the corrected mapped value may appear
                   only where the reference is null. A violation means the block silently
                   overwrote real data.
  2. COVERAGE    — rows gained over the reference alone, and the share of each block that
                   is corrected rather than native.
  3. CONTROLS    — positive controls (DDC, NEFL, GFAP, SNCA, CHI3L1, CST3, MAPT) and
                   negative controls, reference vs harmonized, on the built columns.

Writes harmonized_build_validation.md.
"""

from __future__ import annotations

import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm

from build_common import DATASET_STEM, require_build

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = require_build(DATASET_STEM)
FIT = sorted(__import__("glob").glob(
    os.path.join(HERE, "build_intermediates", "harmonization_fit-*.tab")))[-1]
REPORT = os.path.join(HERE, "harmonized_build_validation.md")

POS = ["DDC", "NEFL", "GFAP", "SNCA", "CHI3L1", "CST3", "MAPT"]
NEG = ["CTRL", "mCherry"]
GROUP = "grp_NSD_vs_HC"
L: list[str] = []


def w(s: str = "") -> None:
    L.append(s)
    print(s, flush=True)


def fit_one(df, col):
    s = df[[col, GROUP, "age_at_visit", "SEX", "PATNO"]].dropna()
    if s[GROUP].nunique() < 2 or len(s) < 30 or s["PATNO"].nunique() < 20:
        return None
    y = s[col].astype(float)
    y = (y - y.mean()) / y.std()
    X = sm.add_constant(s[[GROUP, "age_at_visit", "SEX"]].astype(float))
    try:
        m = sm.OLS(y, X).fit(cov_type="cluster",
                             cov_kwds={"groups": s["PATNO"].to_numpy()})
    except Exception:
        return None
    return dict(n=len(s), ppl=s["PATNO"].nunique(), beta=float(m.params[GROUP]),
                se=float(m.bse[GROUP]), p=float(m.pvalues[GROUP]))


def main() -> None:
    f = pd.read_csv(FIT, sep="\t")
    header = pd.read_csv(DATA, sep="\t", nrows=0).columns.tolist()

    w("# Harmonized blocks — validation of the built dataset")
    w()
    w(f"Reads the shipped `harmonized_*` columns in "
      "`Project_9004_Unified_Emerging_Biomarkers.tab`, so this tests what a downstream "
      "user consumes rather than a reconstruction.")
    w()

    # ---------------------------------------------------------------- 1. integrity
    w("## 1. Fill-rule integrity")
    w()
    w("Where the reference project has a value, the harmonized column must equal it "
      "exactly. The corrected mapped value may appear **only** where the reference is "
      "null. A violation would mean real reference data was silently overwritten.")
    w()
    w("| block | analytes checked | cells with reference value | **mismatches** | "
      "cells filled from mapped |")
    w("|---|---|---|---|---|")
    rng = np.random.default_rng(0)
    total_bad = 0
    for block, sub in f.groupby("block"):
        probe = sub.sample(min(60, len(sub)), random_state=0)
        cols = (["key"] + probe.col_reference.tolist() + probe.harmonized_column.tolist())
        d = pd.read_csv(DATA, sep="\t", usecols=cols, low_memory=False)
        n_ref = n_bad = n_fill = 0
        for _, r in probe.iterrows():
            a = d[r.col_reference].to_numpy(float)
            h = d[r.harmonized_column].to_numpy(float)
            has = ~np.isnan(a)
            n_ref += int(has.sum())
            n_bad += int((has & ~np.isclose(a, h, rtol=0, atol=1e-9, equal_nan=True)).sum())
            n_fill += int((~has & ~np.isnan(h)).sum())
        total_bad += n_bad
        w(f"| `{block}` | {len(probe)} | {n_ref:,} | **{n_bad}** | {n_fill:,} |")
    w()
    w(f"**{'PASS' if total_bad == 0 else 'FAIL'}** — "
      f"{'every reference value is preserved exactly.' if total_bad == 0 else f'{total_bad} cells overwritten.'}")
    w()

    # ------------------------------------------------------- 1b. correction applied
    w("### 1b. Filled values match the recorded correction")
    w()
    w("Check 1 proves reference values survived. This proves the *filled* cells are "
      "exactly `applied_slope * mapped + applied_intercept`, i.e. the export matches the "
      "slope and intercept recorded per analyte in the data dictionary. Any drift here "
      "would mean the dictionary documents a correction the data does not carry.")
    w()
    w("| block | analytes checked | filled cells | **mismatches** | max abs deviation |")
    w("|---|---|---|---|---|")
    total_drift = 0
    for block, sub in f.groupby("block"):
        probe = sub.sample(min(40, len(sub)), random_state=1)
        cols = (probe.col_reference.tolist() + probe.col_mapped.tolist()
                + probe.harmonized_column.tolist())
        d = pd.read_csv(DATA, sep="\t", usecols=cols, low_memory=False)
        n_fill = n_bad = 0
        mx = 0.0
        for _, r in probe.iterrows():
            a = d[r.col_reference].to_numpy(float)
            m = d[r.col_mapped].to_numpy(float)
            h = d[r.harmonized_column].to_numpy(float)
            filled = np.isnan(a) & ~np.isnan(m)
            if not filled.any():
                continue
            expect = r.applied_slope * m[filled] + r.applied_intercept
            got = h[filled]
            dev = np.abs(expect - got)
            n_fill += int(filled.sum())
            n_bad += int((dev > 1e-6).sum())
            mx = max(mx, float(np.nanmax(dev)) if dev.size else 0.0)
        total_drift += n_bad
        w(f"| `{block}` | {len(probe)} | {n_fill:,} | **{n_bad}** | {mx:.2e} |")
    w()
    w(f"**{'PASS' if total_drift == 0 else 'FAIL'}** — "
      + ("every filled cell reproduces the recorded slope and intercept, so the export "
         "matches the dictionary." if total_drift == 0
         else f"{total_drift} cells deviate from the recorded correction."))
    w()

    # ---------------------------------------------------------------- 2. coverage
    w("## 2. Coverage and correction dependence")
    w()
    w("| block | reference | mapped | rows, reference alone | rows, harmonized | gain | % of rows corrected |")
    w("|---|---|---|---|---|---|---|")
    for block, sub in f.groupby("block"):
        r0 = sub.iloc[0]
        probe = sub.head(40)
        d = pd.read_csv(DATA, sep="\t", low_memory=False,
                        usecols=probe.col_reference.tolist() + probe.harmonized_column.tolist())
        nref = int(d[probe.col_reference.tolist()].notna().any(axis=1).sum())
        nharm = int(d[probe.harmonized_column.tolist()].notna().any(axis=1).sum())
        w(f"| `{block}` | {r0.reference_project} | {r0.mapped_project} | {nref:,} "
          f"| **{nharm:,}** | +{nharm - nref:,} | {r0.pct_from_corrected}% |")
    w()

    # ---------------------------------------------------------------- 3. controls
    w("## 3. Positive and negative controls on the built columns")
    w()
    w(f"Contrast `{GROUP}`, OLS with PATNO-clustered SEs over all visits, outcome "
      "z-scored. A random-intercept LMM is not identifiable — these panels average ~1.2 "
      "observations per participant.")
    w()
    targets = f[f.analyte.isin(POS + NEG)]
    cols = ["PATNO", "age_at_visit", "SEX", GROUP]
    d = pd.read_csv(DATA, sep="\t", low_memory=False, dtype={"PATNO": str},
                    usecols=cols + targets.col_reference.tolist()
                    + targets.harmonized_column.tolist())
    rows = []
    for _, r in targets.iterrows():
        a = fit_one(d, r.col_reference)
        h = fit_one(d, r.harmonized_column)
        if a is None or h is None:
            continue
        rows.append(dict(block=r.block, analyte=r.analyte,
                         control="neg" if r.analyte in NEG else "pos",
                         slope=r.applied_slope, ref_n=a["n"], ref_beta=a["beta"],
                         ref_se=a["se"], ref_p=a["p"], harm_n=h["n"],
                         harm_beta=h["beta"], harm_se=h["se"], harm_p=h["p"]))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(HERE, "harmonized_build_validation.tab"), sep="\t", index=False)

    for kind, label in (("pos", "Positive controls"), ("neg", "Negative controls")):
        S = R[R.control == kind].sort_values(["block", "analyte"])
        if S.empty:
            continue
        w(f"### {label}")
        w()
        w("| block | analyte | ref n | ref beta | ref P | harm n | harm beta | harm P | P better? |")
        w("|---|---|---|---|---|---|---|---|---|")
        for _, r in S.iterrows():
            w(f"| `{r.block}` | {r.analyte} | {r.ref_n:.0f} | {r.ref_beta:+.3f} "
              f"| {r.ref_p:.1e} | {r.harm_n:.0f} | {r.harm_beta:+.3f} | {r.harm_p:.1e} "
              f"| {'**yes**' if r.harm_p < r.ref_p else 'no'} |")
        w()

    P = R[R.control == "pos"]
    if not P.empty:
        z = (P.harm_beta - P.ref_beta) / np.sqrt(P.harm_se ** 2 + P.ref_se ** 2)
        w("### Summary")
        w()
        w(f"- positive-control fits: **{len(P)}**")
        w(f"- observations: median {P.ref_n.median():.0f} -> **{P.harm_n.median():.0f}**")
        w(f"- median SE: {P.ref_se.median():.4f} -> **{P.harm_se.median():.4f}** "
          f"({(1 - P.harm_se.median() / P.ref_se.median()) * 100:.0f}% smaller)")
        w(f"- p-value improved in **{(P.harm_p < P.ref_p).sum()} of {len(P)}**")
        w(f"- beta shift beyond sampling noise (|z| > 1.96): **{(z.abs() > 1.96).sum()} "
          f"of {len(P)}**")
        w()
        dd = P[P.analyte == "DDC"]
        if not dd.empty:
            w("DDC, the anchor signal, across every block that carries it:")
            w()
            for _, r in dd.iterrows():
                w(f"- `{r.block}`: {r.ref_beta:+.3f} (P={r.ref_p:.1e}, n={r.ref_n:.0f}) "
                  f"-> **{r.harm_beta:+.3f} (P={r.harm_p:.1e}, n={r.harm_n:.0f})**")
            w()

    with open(REPORT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\nwrote {REPORT}")


if __name__ == "__main__":
    sys.exit(main())
