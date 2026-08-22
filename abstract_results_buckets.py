#!/usr/bin/env python3
"""abstract_results_buckets.py — regroup the robust meta-analysis results onto the
four questions the stakeholders asked, which is how the AD/PD abstract now reports them.

The batch is organised by MODEL (logit / ols / cox / lmm); the stakeholders asked four
QUESTIONS that cut across those models:

  1. Baseline proteins differing by NSD status within each PPMI cohort, by genotype
     (LRRK2, GBA1), or by NSD-ISS stage.
  2. Baseline proteins separating genetic PD (LRRK2-PD, GBA-PD) from sporadic PD --
     the `_SAAadj` variants add NSD status to that comparison.
  3. Baseline proteins predicting progression (motor, cognition, milestones) or DaT.
  4. Proteins changing longitudinally across the same cohorts as (1).

Every run in batch.yaml is assigned to exactly one bucket, including the ones that fall
OUTSIDE the four questions -- the assignment asserts exhaustiveness, so a run added to
batch.yaml without a bucket fails here rather than silently vanishing from the abstract.

Robust is the meta-analysis definition and is not recomputed here: Bonferroni-significant
on the family threshold, replicated in both ancestry strata, direction-concordant and not
stratum-dominated (`robust_FE` in the META_*.csv).

Two counting rules that the abstract depends on:

  * Trajectory runs report two terms per analyte. `within` is the within-subject
    time x group interaction -- an actual difference in rate of change, which is what
    question 4 asks for. `between` is the subject-mean time term, which carries cohort
    composition rather than change. They are counted separately.
  * Rows are analyte-level, so one protein measured in CSF and plasma on two platforms
    can contribute up to four rows to one run. Protein counts are therefore always
    reported alongside row counts.

Reads the newest META_*.csv per run name -- the same discovery rule the app's
snapshot_data.py uses, so the abstract, the browser and this table cannot disagree.

Usage:  python3 abstract_results_buckets.py [--meta-dir meta] [-o abstract_results_buckets.md]
"""
from __future__ import annotations

import argparse
import glob
import os
import re

import pandas as pd

_META_RE = re.compile(r"^META_(?P<name>.+)-(?P<ts>\d{8}_\d{6})\.csv$")

Q1 = "Q1 baseline: NSD status, genotype, NSD-ISS stage"
Q2 = "Q2 baseline: genetic PD vs sporadic PD"
Q3 = "Q3 baseline predictors of progression"
Q4 = "Q4 longitudinal change"
OUT_CI = "Outside: prevalent cognitive impairment at baseline"
OUT_PRS = "Outside: PRS157 interaction runs"

BUCKETS: dict[str, list[str]] = {
    Q1: [
        "NSD_vs_HC", "NSD_vs_notNSD_sPD", "NSD_vs_notNSD_LRRK2", "NSD_vs_notNSD_GBA",
        "NSD_vs_notNSD_prodromal", "NSD_stage_2A_vs_2B", "NSD_stage_2A_vs_2B_PDonly",
        "NSD_stage_2A_vs_3", "NSD_stage_2A_vs_3_PDonly", "NSD_stage_2B_vs_3",
        "NSD_stage_2B_vs_3_PDonly", "NSD_stage_early_vs_late",
    ],
    Q2: ["sPD_vs_LRRK2", "sPD_vs_LRRK2_SAAadj", "sPD_vs_GBA", "sPD_vs_GBA_SAAadj"],
    Q3: [
        "slope_moca", "slope_updrs3_off", "slope_updrs3_off_LEDD_adj", "slope_lowput_ratio",
        "cox_moca_lt26", "cox_cogstate_worsen", "cox_pm_any", "cox_pm_cog_any",
        "cox_pm_mc_any", "cox_stage_d", "cox_nsd_2a_to_later", "cox_nsd_2b_to_later",
    ],
    Q4: [
        "trajectory_NSD_vs_HC", "trajectory_HC_vs_Prodromal", "trajectory_HC_vs_PD",
        "trajectory_Prodromal_vs_PD", "trajectory_HCNSDneg_vs_ProdNSDpos",
        "trajectory_HCNSDneg_vs_PDNSDpos", "trajectory_ProdNSDpos_vs_PDNSDpos",
        "trajectory_stage_2A_vs_2B", "trajectory_CI_MOCA", "trajectory_CI_PI",
    ],
    # Cross-sectional CI at baseline is neither a group contrast nor a prediction, so it
    # answers none of the four as asked.
    OUT_CI: ["CI_PI_baseline", "CI_MOCA_baseline"],
    # PRS157 carries LRRK2 and GBA variants, so a protein x PRS interaction is partly a
    # protein x carrier interaction. Kept out of the abstract's headline counts.
    OUT_PRS: [
        "slope_moca_x_PRS", "slope_updrs3_off_x_PRS", "slope_updrs3_off_x_PRS_LEDD_adj",
        "slope_lowput_ratio_x_PRS", "cox_moca_lt26_x_PRS", "cox_cogstate_worsen_x_PRS",
        "cox_pm_any_x_PRS", "cox_pm_cog_any_x_PRS", "cox_pm_mc_any_x_PRS",
        "cox_stage_d_x_PRS", "cox_nsd_2a_to_later_x_PRS", "cox_nsd_2b_to_later_x_PRS",
    ],
}
FOUR = [Q1, Q2, Q3, Q4]

# PRS157 contains LRRK2 and GBA variants, so it partly proxies carrier status: a PRS
# result read off a contrast whose groups are DEFINED by carrier status is circular.
# These are those runs; PRS results in them are excluded from anything reported.
CARRIER_RUNS = ["sPD_vs_LRRK2", "sPD_vs_LRRK2_SAAadj", "sPD_vs_GBA", "sPD_vs_GBA_SAAadj",
                "NSD_vs_notNSD_LRRK2", "NSD_vs_notNSD_GBA"]
PRS_COL = "p9005_Genetic_PRS_PRS157"


def prs_summary(rob: pd.DataFrame) -> dict:
    """Robust PRS157 main effects, split into carrier-confounded and reportable."""
    prs = rob[rob["analyte"] == PRS_COL]
    clean = prs[~prs["run"].isin(CARRIER_RUNS)]
    nsd = clean[clean["run"] == "NSD_vs_HC"]
    return {
        "total": len(prs),
        "carrier": len(prs) - len(clean),
        "clean": len(clean),
        "runs": sorted(clean["run"]),
        "beta_nsd_vs_hc": float(nsd["beta_FE"].iloc[0]) if len(nsd) else float("nan"),
        "p_nsd_vs_hc": float(nsd["P_FE"].iloc[0]) if len(nsd) else float("nan"),
    }


# ---------------------------------------------------------------- the numbers
# Plain-language rendering of each stakeholder question, how many hits Panel B has room
# for, and which term answers it. Order matters: it is the order the abstract's RESULTS
# section reads in.
#
# Q4's `term` is why it is here: a trajectory run reports a within-subject rate term and
# a between-subject level term, and only the within term is a protein CHANGING. The row's
# count covers every robust result those analyses returned; the named proteins are the
# ones that answer the question as it was asked.
QUESTION_TEXT = {
    Q1: ("Who differs at baseline",
         "NSD status within each cohort, LRRK2 / GBA1 genotype, and NSD-ISS stage", 6, None),
    Q2: ("Genetic PD vs sporadic PD",
         "LRRK2-PD and GBA-PD against sporadic PD, with and without NSD status", 3, None),
    Q3: ("What predicts later progression",
         "baseline protein vs motor, cognitive, milestone and DaT outcomes", 6, None),
    Q4: ("What changes over time",
         "trajectories across HC, prodromal and PD, and across NSD-ISS stages", 6, "within"),
}


def hit_chips(prot: "pd.DataFrame", k: int) -> list[tuple[str, int, str]]:
    """Top k proteins for one question: (symbol, robust results, fluid).

    Ranked by how many robust results the protein carries -- recurrence across blocks is
    what `robust` is trying to buy, so it leads -- with ties broken by best P. The single
    strongest result in the question is then forced in if that ranking dropped it, which
    is what keeps SPAG7 (P=1e-56, seen once) on the slide next to DDC (seen sixteen times).

    `fluid` is 'both' when the protein is robust in CSF and in plasma; the chip is
    coloured by it.
    """
    if prot.empty:
        return []
    g = prot.groupby("protein").agg(n=("protein", "size"), p=("P_FE", "min"),
                                    fluids=("fluid", lambda s: set(s)))
    g = g.sort_values(["n", "p"], ascending=[False, True])
    top = g.head(k)
    best = g["p"].idxmin()
    if best not in top.index:
        top = pd.concat([top.head(k - 1), g.loc[[best]]])
    out = []
    for name, r in top.iterrows():
        fl = "both" if {"CSF", "plasma"} <= r["fluids"] else next(iter(r["fluids"]))
        out.append((name, int(r["n"]), fl))
    return out


def protein_symbol(analyte: object) -> str:
    """`harmonized_olink_csf_DDC_NPX` -> `DDC`; non-proteomic names pass through."""
    s = re.sub(r"^harmonized_(olink|nulisa_cns|nulisa_inf)_(csf|plasma)_", "",
               str(analyte), flags=re.IGNORECASE)
    s = re.sub(r"^p\d+_(CNS_plasma|CNS_CSF|Inflammation_CSF|Inflammation_Plasma|"
               r"Inflammation_plasma|Neuro_CSF|Neuro_Plasma|olink_plasma|CSF|Plasma)_", "", s)
    return re.sub(r"_(NPX|NPQ)$", "", s)


def fluid_of(analyte: object) -> str:
    s = str(analyte).lower()
    return "CSF" if "_csf_" in s else ("plasma" if "plasma" in s else "-")


def platform_of(analyte: object) -> str:
    s = str(analyte).lower()
    return "Olink" if "olink" in s else ("NULISA" if "nulisa" in s else "-")


def newest_per_run(meta_dir: str) -> dict[str, str]:
    """Latest timestamp per run name -- snapshot_data.py's discovery rule."""
    latest: dict[str, tuple[str, str]] = {}
    for path in sorted(glob.glob(os.path.join(meta_dir, "META_*.csv"))):
        m = _META_RE.match(os.path.basename(path))
        if not m:
            continue
        name, ts = m.group("name"), m.group("ts")
        if name not in latest or ts > latest[name][0]:
            latest[name] = (ts, path)
    return {n: p for n, (_, p) in latest.items()}


def load(meta_dir: str) -> pd.DataFrame:
    paths = newest_per_run(meta_dir)
    if not paths:
        raise SystemExit(f"ERROR: no META_*.csv found in {meta_dir}")
    frames = []
    for run, path in sorted(paths.items()):
        d = pd.read_csv(path, low_memory=False)
        d["run"] = run
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)

    assigned = {r: b for b, runs in BUCKETS.items() for r in runs}
    unmapped = sorted(set(df["run"]) - set(assigned))
    if unmapped:
        raise SystemExit("ERROR: runs with no stakeholder bucket: " + ", ".join(unmapped))
    df["bucket"] = df["run"].map(assigned)
    df["protein"] = df["analyte"].map(protein_symbol)
    df["fluid"] = df["analyte"].map(fluid_of)
    df["platform"] = df["analyte"].map(platform_of)
    return df


def summarise(rob: pd.DataFrame, df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for b, runs in BUCKETS.items():
        sub = rob[rob["bucket"] == b]
        prot = sub[sub["is_proteomic"]]
        rows.append({
            "Question": b,
            "Analyses": len(runs),
            "Analyses with a robust result": sub["run"].nunique(),
            "Robust results": len(sub),
            "Robust proteomic": len(prot),
            "Distinct proteins": prot["protein"].nunique(),
        })
    return pd.DataFrame(rows)


def md_table(d: pd.DataFrame) -> str:
    head = "| " + " | ".join(d.columns) + " |"
    rule = "|" + "|".join(["---" if i == 0 else "---:" for i in range(len(d.columns))]) + "|"
    body = ["| " + " | ".join(str(v) for v in r) + " |" for r in d.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def top_hits(sub: pd.DataFrame, n: int = 8) -> pd.DataFrame:
    d = sub.sort_values("P_FE").head(n).copy()
    d["Analyte"] = [
        p if f == "-" else f"{p} ({f}, {pl})"
        for p, f, pl in zip(d["protein"], d["fluid"], d["platform"])
    ]
    d["beta (FE)"] = d["beta_FE"].map(lambda v: f"{v:+.2f}")
    d["P (FE)"] = d["P_FE"].map(lambda v: f"{v:.0e}".replace("e-0", "e-"))
    d["N"] = d["N_total"].map(lambda v: f"{int(v):,}")
    cols = ["run", "Analyte", "term", "beta (FE)", "P (FE)", "N"]
    return d[cols].rename(columns={"run": "Analysis", "term": "Term"})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta-dir", default="meta")
    ap.add_argument("-o", "--output", default="abstract_results_buckets.md")
    ap.add_argument("--csv", default="abstract_results_buckets.csv")
    ap.add_argument("--top", type=int, default=8, help="hits listed per question")
    a = ap.parse_args()

    df = load(a.meta_dir)
    rob = df[df["robust_FE"]].copy()
    four = rob[rob["bucket"].isin(FOUR)]
    fourp = four[four["is_proteomic"]]

    out = ["# Robust results by stakeholder question",
           "",
           f"Generated by `abstract_results_buckets.py` from `{a.meta_dir}/META_*.csv` "
           f"(newest timestamp per run).",
           "",
           f"{len(df):,} associations tested across {df['run'].nunique()} analyses; "
           f"**{len(rob)}** are robust. **{len(four)}** of those "
           f"({len(fourp)} proteomic, {fourp['protein'].nunique()} distinct proteins) answer one "
           f"of the four questions, across {four['run'].nunique()} analyses.",
           "",
           md_table(summarise(rob, df)),
           ""]

    within = four[(four["bucket"] == Q4) & (four["term"] == "within")]
    between = four[(four["bucket"] == Q4) & (four["term"] == "between")]
    out += [f"Q4 splits into **{len(within)}** within-subject rate-of-change terms "
            f"({within['protein'].nunique()} proteins) and {len(between)} between-subject "
            f"terms ({between['protein'].nunique()} proteins). Only the within-subject term "
            "is a difference in how fast a protein moves.",
            ""]

    for b in FOUR:
        sub = rob[rob["bucket"] == b]
        prot = sub[sub["is_proteomic"]]
        out += [f"## {b}", ""]
        if sub.empty:
            out += ["No robust result.", ""]
            continue
        per_run = ", ".join(f"`{r}` {n}" for r, n in sub["run"].value_counts().items())
        recur = ", ".join(f"{p} ({n})" for p, n in prot["protein"].value_counts().head(6).items())
        empty = [r for r in BUCKETS[b] if r not in set(sub["run"])]
        out += [f"Robust per analysis: {per_run}.", ""]
        out += [f"Most recurrent proteins: {recur}.", ""]
        if empty:
            out += ["Returned nothing robust: " + ", ".join(f"`{r}`" for r in empty) + ".", ""]
        out += [md_table(top_hits(sub, a.top)), ""]

    with open(a.output, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))

    keep = ["bucket", "run", "analyte", "protein", "fluid", "platform", "term",
            "beta_FE", "SE_FE", "P_FE", "N_total", "I2", "beta_EUR", "P_EUR",
            "beta_AJ", "P_AJ", "is_proteomic"]
    rob.sort_values(["bucket", "run", "P_FE"])[keep].to_csv(a.csv, index=False)

    print("\n".join(out[:12]))
    print(f"\nWrote {a.output} and {a.csv} ({len(rob)} robust rows)")


if __name__ == "__main__":
    main()
