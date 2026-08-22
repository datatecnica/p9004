#!/usr/bin/env python3
"""validate_abstract_numbers.py — check every number in `adpd_abstract.md` against the
most recent meta-analyses in `meta/`, and confirm the app's snapshot agrees with them.

The abstract is written from the phase 1 wrap meta-analysis output (newest META_*.csv per
run). The browser app ships a snapshot of that same output, so the two must agree; this
script is the thing that says so. Failures print the expected and observed value.

Three sources are consulted:

  meta/META_*.csv     newest per run -- the authority for every result-level number
  meta_analysis_*.log matching timestamp -- the per run x term genomic-inflation cells
  app data/           optional cross-check; the snapshot must reproduce the robust set

The app filters `thin_within` rows out of its browsable table (LMM within-subject terms
with too few repeat observations to interpret). They are counted here so the row totals
reconcile rather than looking like a discrepancy: meta rows - thin_within = app rows.

Usage:
  python3 validate_abstract_numbers.py [--meta-dir meta] [--app-dir ../../proteomics_data_mine_app/data]
Exit status is 1 if any check fails.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

from abstract_results_buckets import (BUCKETS, FOUR, Q1, Q2, Q3, Q4, QUESTION_TEXT,
                                      hit_chips, load, newest_per_run, prs_summary)

CHECKS: list[tuple[str, object, object]] = []


def check(name: str, expected, observed, *, tol: float = 0.0) -> None:
    if isinstance(expected, float) and isinstance(observed, (int, float)):
        ok = abs(float(observed) - expected) <= tol
    else:
        ok = expected == observed
    CHECKS.append((name, expected, observed if ok else f"{observed}  <-- MISMATCH"))
    if not ok:
        CHECKS[-1] = (name, expected, f"{observed}  <-- MISMATCH")


def lambda_cells(meta_dir: str) -> np.ndarray:
    """λ_FE per run x reported term, parsed from the log of the run that made meta/."""
    ts = sorted({re.search(r"-(\d{8}_\d{6})\.csv$", p).group(1)
                 for p in newest_per_run(meta_dir).values()})[-1]
    logs = glob.glob(os.path.join(os.path.dirname(meta_dir) or ".",
                                  f"meta_analysis_{ts}.log"))
    if not logs:
        return np.array([])
    txt = open(logs[0], encoding="utf-8", errors="replace").read()
    return np.array([float(m) for m in re.findall(r"\[\w+\]\s+λ_FE = ([0-9.]+)", txt)])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta-dir", default="meta")
    ap.add_argument("--app-dir", default="../../proteomics_data_mine_app/data")
    a = ap.parse_args()

    df = load(a.meta_dir)
    rob = df[df["robust_FE"]].copy()
    four = rob[rob["bucket"].isin(FOUR)]
    fourp = four[four["is_proteomic"]]

    def q(bucket: str) -> pd.DataFrame:
        return rob[rob["bucket"] == bucket]

    def hit(run: str, protein: str, term: str = "effect") -> pd.DataFrame:
        return rob[(rob["run"] == run) & (rob["protein"] == protein) & (rob["term"] == term)]

    # ---------------------------------------------------------------- scale of the batch
    thin = int(df["thin_within"].sum())
    check("associations tested (meta rows - thin_within)", 711_927, len(df) - thin)
    check("analyses run", 52, df["run"].nunique())
    check("robust results", 296, len(rob))
    check("robust answering the four questions", 280, len(four))
    check("analyses with a robust answer to the four", 28, four["run"].nunique())
    check("distinct proteins across the four", 136, fourp["protein"].nunique())

    # ---------------------------------------------------------------- Q1
    q1 = q(Q1)
    check("Q1 robust results", 88, len(q1))
    check("Q1 distinct proteins", 53, q1[q1["is_proteomic"]]["protein"].nunique())
    ddc = hit("NSD_vs_HC", "DDC").sort_values("P_FE")
    check("Q1 DDC NSD_vs_HC raised (beta > 0)", True, float(ddc["beta_FE"].iloc[0]) > 0)
    check("Q1 DDC NSD_vs_HC quoted P (2e-40)", "2e-40", f'{float(ddc["P_FE"].iloc[0]):.0e}')
    check("Q1 DDC NSD_vs_HC fluids", {"CSF", "plasma"}, set(ddc["fluid"]))
    check("Q1 DDC NSD_vs_HC platforms", {"Olink", "NULISA"}, set(ddc["platform"]))
    stage = q1[q1["run"].str.startswith("NSD_stage")]
    check("Q1 DDC rises with stage (all beta > 0)", True,
          bool((stage[stage["protein"] == "DDC"]["beta_FE"] > 0).all())
          and len(stage[stage["protein"] == "DDC"]) > 0)
    cd276 = stage[stage["protein"] == "CD276"]
    check("Q1 CD276 stage analyses", 3, cd276["run"].nunique())
    check("Q1 CD276 stage betas all positive", True, bool((cd276["beta_FE"] > 0).all()))
    gba = q1[q1["run"] == "NSD_vs_notNSD_GBA"]
    for p in ("SPAG7", "LIMD1"):
        g = gba[gba["protein"] == p]
        check(f"Q1 {p} robust in GBA NSD contrast, CSF, positive", True,
              len(g) == 1 and g["fluid"].iloc[0] == "CSF" and float(g["beta_FE"].iloc[0]) > 0)

    # ---------------------------------------------------------------- Q2
    q2 = q(Q2)
    check("Q2 robust results", 11, len(q2))
    check("Q2 distinct proteins", 3, q2[q2["is_proteomic"]]["protein"].nunique())
    eno2 = hit("sPD_vs_LRRK2", "ENO2")
    check("Q2 ENO2 sPD_vs_LRRK2 lower (beta < 0)", True, float(eno2["beta_FE"].iloc[0]) < 0)
    check("Q2 ENO2 quoted P (3e-05)", "3e-05", f'{float(eno2["P_FE"].iloc[0]):.0e}')
    eno2_adj = hit("sPD_vs_LRRK2_SAAadj", "ENO2")
    check("Q2 ENO2 survives NSD-status adjustment (same sign)", True,
          len(eno2_adj) == 1 and float(eno2_adj["beta_FE"].iloc[0]) < 0)

    # ---------------------------------------------------------------- Q3
    q3 = q(Q3)
    check("Q3 robust results", 35, len(q3))
    check("Q3 NPTX2 protective on MoCA slope (beta > 0)", True,
          float(hit("slope_moca", "NPTX2")["beta_FE"].iloc[0]) > 0)
    for run in ("cox_moca_lt26", "cox_pm_cog_any"):
        h = hit(run, "NPTX2")
        check(f"Q3 low NPTX2 raises risk in {run} (beta < 0)", True,
              len(h) == 1 and float(h["beta_FE"].iloc[0]) < 0)
    ptau_slope = q3[(q3["run"] == "slope_moca") & q3["protein"].str.contains("Tau", case=False)]
    check("Q3 phospho-tau predicts faster MoCA decline (all beta < 0)", True,
          len(ptau_slope) > 0 and bool((ptau_slope["beta_FE"] < 0).all()))
    ptau_event = q3[q3["run"].isin(["cox_moca_lt26", "cox_pm_cog_any"])
                    & q3["protein"].str.contains("Tau", case=False)]
    check("Q3 phospho-tau raises event risk (all beta > 0)", True,
          len(ptau_event) > 0 and bool((ptau_event["beta_FE"] > 0).all()))
    ddc3 = q3[q3["protein"] == "DDC"]
    check("Q3 DDC predicts milestones and stage D", {"cox_pm_any", "cox_pm_cog_any", "cox_stage_d"},
          set(ddc3["run"]))
    check("Q3 DDC progression betas all positive", True, bool((ddc3["beta_FE"] > 0).all()))
    check("Q3 OFF-motor slope robust results", 0,
          len(q3[q3["run"].isin(["slope_updrs3_off", "slope_updrs3_off_LEDD_adj"])]))
    check("Q3 DaT slope robust results", 1, len(q3[q3["run"] == "slope_lowput_ratio"]))

    # ---------------------------------------------------------------- Q4
    q4 = q(Q4)
    within = q4[q4["term"] == "within"]
    check("Q4 robust results", 146, len(q4))
    check("Q4 contrasts with a robust result", 9, q4["run"].nunique())
    check("Q4 within-subject rate differences", 75, len(within))
    check("Q4 within-subject proteins", 32, within["protein"].nunique())
    check("Q4 DDC within-subject contrasts", 6, within[within["protein"] == "DDC"]["run"].nunique())
    check("Q4 DDC is the most recurrent within-subject protein", "DDC",
          within["protein"].value_counts().index[0])
    nxt = within[within["protein"] != "DDC"]["protein"].value_counts()
    check("Q4 next four within-subject proteins", ["FLT1", "AOC3", "NEFL", "VEGFD"],
          sorted(nxt.head(4).index, key=lambda p: (-nxt[p], p)))

    # ---------------------------------------------------------------- PRS, in one clause
    prs = prs_summary(rob)
    check("PRS157 robust results", 8, prs["total"])
    check("PRS157 dropped as carrier contrasts", 2, prs["carrier"])
    check("PRS157 reportable results", 6, prs["clean"])
    check("PRS157 NSD_vs_HC rises (beta > 0)", True, prs["beta_nsd_vs_hc"] > 0)
    check("PRS157 NSD_vs_HC quoted P (6e-13)", "6e-13", f'{prs["p_nsd_vs_hc"]:.0e}')
    check("PRS157 rises with stage and predicts stage D", True,
          all(r in prs["runs"] for r in
              ["NSD_stage_2A_vs_3", "NSD_stage_2B_vs_3", "NSD_stage_early_vs_late",
               "cox_stage_d"]))
    check("PRS interaction results excluded from every count", 15,
          len(rob[rob["bucket"] == "Outside: PRS157 interaction runs"]))

    # ---------------------------------------------------------------- conclusions
    ddc = rob[rob["protein"] == "DDC"]
    check("DDC robust in all four questions", 4, ddc["bucket"].nunique())
    check("DDC robust results (figure panel C)", 51, len(ddc))
    check("DDC per question (figure panel C)", [16, 1, 3, 31],
          [len(ddc[ddc["bucket"] == qkey]) for qkey in FOUR])
    check("Q1 DDC in six of the eight analyses that returned anything", (6, 8),
          (ddc[ddc["bucket"] == Q1]["run"].nunique(), q(Q1)["run"].nunique()))
    check("Q2 proteomic robust results (of 11)", 4,
          int(q(Q2)["is_proteomic"].sum()))
    lam = lambda_cells(a.meta_dir)
    check("genomic-inflation cells", 62, len(lam))
    if len(lam):
        check("median lambda", 1.10, round(float(np.median(lam)), 2), tol=0.005)
        check("lambda IQR", (1.02, 1.19),
              (round(float(np.percentile(lam, 25)), 2), round(float(np.percentile(lam, 75)), 2)))
        check("lambda cells within 1.2", 76, int(round(100 * float((lam <= 1.2).mean()))))

    # ---------------------------------------------------------------- figure panel B
    # The proteins named on the slide are also listed in the document's Panel B table, so
    # the two are checked against one selection rather than against each other by eye.
    md = "adpd_abstract.md"
    if os.path.exists(md):
        rows = re.findall(r"^\| ([1-4]) \| .+?\|\s*([^|]+?)\s*\|$",
                          open(md, encoding="utf-8").read(), flags=re.M)
        listed = {int(i): [c.strip().split(" ")[0] for c in cell.split("·")]
                  for i, cell in rows}
        for i, qkey in enumerate(FOUR, start=1):
            sub = rob[(rob["bucket"] == qkey) & rob["is_proteomic"]]
            _, _, k, term = QUESTION_TEXT[qkey]
            named = sub if term is None else sub[sub["term"] == term]
            check(f"Panel B table matches the figure's Q{i} chips",
                  [c[0] for c in hit_chips(named, k)], listed.get(i))

    # ---------------------------------------------------------------- app agreement
    app_parquet = os.path.join(a.app_dir, "results.parquet")
    if os.path.exists(app_parquet):
        app = pd.read_parquet(app_parquet)
        key = lambda d: set(zip(d["run"], d["analyte"], d["term"]))
        check("app snapshot robust count", len(rob), int(app["robust_FE"].sum()))
        check("app snapshot robust set identical to meta", True,
              key(app[app["robust_FE"]]) == key(rob))
        check("app rows = meta rows - thin_within", len(df) - thin, len(app))
        man = os.path.join(a.app_dir, "manifest.json")
        if os.path.exists(man):
            m = json.load(open(man))
            check("app manifest n_robust", len(rob), m["n_robust"])
            check("app manifest n_runs", df["run"].nunique(), m["n_runs"])
            check("app manifest n_rows", len(df) - thin, m["n_rows"])
    else:
        CHECKS.append(("app snapshot", "-", f"not found at {app_parquet}, skipped"))

    width = max(len(n) for n, _, _ in CHECKS)
    failed = 0
    for name, exp, obs in CHECKS:
        bad = isinstance(obs, str) and "MISMATCH" in obs
        failed += bad
        print(f"{'FAIL' if bad else 'ok  '}  {name:<{width}}  expected {exp!s:<28} got {obs}")
    print(f"\n{len(CHECKS) - failed}/{len(CHECKS)} checks pass")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
