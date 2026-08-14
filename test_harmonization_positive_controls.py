"""Does harmonization work? Longitudinal positive-control test.

REFERENCE CONVENTION: the reference is the project with MORE ROWS, chosen per block
rather than by age. This maximises the share of the harmonized column that is
uncorrected native values and minimises exposure to correction error. On the current
data it makes p314 the reference for both Olink blocks (p314 has 2,098 CSF and 2,942
plasma rows against p277's 1,167 and p293's 784), while the four NULISA blocks keep the
older project, which is also the larger one.

THE TEST. Harmonization adds participants *and* timepoints. In a repeated-measures fit
that extra N should buy precision, so a working harmonization shows smaller p-values
than the reference project alone. If p-values do not improve, either the correction is
injecting noise or the added samples are not carrying the effect.

Model, per analyte: OLS  analyte_z ~ group + age_at_visit + SEX  with PATNO-clustered
standard errors, over ALL visits carrying a value (not one row per participant). The
outcome is z-scored so betas are per-SD and comparable across fits. A random-intercept
MixedLM would be the textbook choice but is not identifiable on these panels — they
average ~1.2 observations per participant, so the random-effect variance cannot be
estimated and the fit is singular.

Fitted three ways:
  reference     the larger project alone           (the baseline to beat)
  mapped        the smaller project alone          (independent estimate)
  harmonized    reference, with corrected mapped filling only where reference is absent

Positive controls: DDC (strongest established signal here), NEFL, GFAP, SNCA, CHI3L1,
CST3, MAPT.  Negative controls: CTRL / mCherry spike-ins, which must not gain
significance.
"""

from __future__ import annotations

import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "Project_9004_Unified_Emerging_Biomarkers.tab")
REPORT = os.path.join(HERE, "harmonization_positive_controls.md")

# (block, project_1, project_2) — reference is chosen at runtime as whichever has more rows
PAIRS = [
    ("olink_plasma",      "p293_olink_plasma",        "p314_Plasma"),
    ("olink_csf",         "p277_CSF",                 "p314_CSF"),
    ("nulisa_cns_plasma", "p288_CNS_plasma",          "p312_Neuro_Plasma"),
    ("nulisa_cns_csf",    "p282_CNS_CSF",             "p312_Neuro_CSF"),
    ("nulisa_inf_plasma", "p288_Inflammation_plasma", "p312_Inflammation_Plasma"),
    ("nulisa_inf_csf",    "p282_Inflammation_CSF",    "p312_Inflammation_CSF"),
]
POS = ["DDC", "NEFL", "GFAP", "SNCA", "CHI3L1", "CST3", "MAPT"]
NEG = ["CTRL", "mCherry"]
GROUP = "grp_NSD_vs_HC"
MIN_CELL = 10
L: list[str] = []


def w(s: str = "") -> None:
    L.append(s)
    print(s, flush=True)


def amap(header, prefix):
    return {re.sub(r"_(NPX|NPQ)$", "", c[len(prefix) + 1:]): c
            for c in header
            if c.startswith(prefix + "_") and "earliest_visit_PC" not in c
            and "PLATE_ID" not in c}


def lmm(df, col):
    """Fit over all visits carrying a value, with PATNO-clustered standard errors.

    A random-intercept MixedLM is the natural choice for repeated measures, but it is
    not identifiable here: these panels average only ~1.2 observations per participant,
    so the random-effect variance is unestimable and statsmodels fails with a singular
    matrix (p314_CSF DDC: 1,199 observations across 1,014 participants). Cluster-robust
    OLS handles the repeated measures correctly without estimating a variance component,
    and degrades gracefully when most clusters are singletons.
    """
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
    if not np.isfinite(m.bse.get(GROUP, np.nan)):
        return None
    return dict(n_obs=len(s), n_ppl=s["PATNO"].nunique(),
                obs_per_ppl=len(s) / s["PATNO"].nunique(),
                beta=float(m.params[GROUP]), se=float(m.bse[GROUP]),
                p=float(m.pvalues[GROUP]))


def main() -> None:
    header = pd.read_csv(DATA, sep="\t", nrows=0).columns.tolist()
    maps, need = {}, set()
    for _, p1, p2 in PAIRS:
        maps[p1], maps[p2] = amap(header, p1), amap(header, p2)
        for t in POS + NEG:
            for m in (maps[p1], maps[p2]):
                if t in m:
                    need.add(m[t])
    # A sample of each project's analytes, purely to establish which rows that project
    # assayed — the target analytes alone would understate coverage.
    presence = {p: sorted(m.values())[:40] for p, m in maps.items()}
    need |= {c for v in presence.values() for c in v}

    base = ["PATNO", "EVENT_ID", "COHORT", "age_at_visit", "SEX", GROUP]
    d = pd.read_csv(DATA, sep="\t", usecols=base + sorted(need), low_memory=False,
                    dtype={"PATNO": str, "EVENT_ID": str, "COHORT": str})
    d["cell"] = d["EVENT_ID"] + "|" + d["COHORT"].fillna("NA")

    # --- choose the reference per block: more rows wins ---
    refs = {}
    for name, p1, p2 in PAIRS:
        n1 = int(d[presence[p1]].notna().any(axis=1).sum())
        n2 = int(d[presence[p2]].notna().any(axis=1).sum())
        refs[name] = ((p1, p2, n1, n2) if n1 >= n2 else (p2, p1, n2, n1))

    w("# Does harmonization work? Longitudinal positive controls")
    w()
    w(f"Contrast `{GROUP}`. OLS `analyte ~ group + age_at_visit + SEX` with "
      "PATNO-clustered SEs, fitted over **all visits** carrying a value. Outcome "
      "z-scored so betas are per-SD. A random-intercept LMM is not identifiable here — "
      "these panels average ~1.2 observations per participant, so the variance "
      "component is unestimable and statsmodels returns a singular matrix.")
    w()
    w("**The reference is the project with more rows**, so the harmonized column is "
      "mostly uncorrected native values and correction error is minimised.")
    w()
    w("| block | **reference** | rows | mapped | rows |")
    w("|---|---|---|---|---|")
    for name, _, _ in PAIRS:
        R, M, nR, nM = refs[name]
        w(f"| `{name}` | **{R}** | {nR:,} | {M} | {nM:,} |")
    w()
    w("Harmonization should add both participants and timepoints, so a working "
      "correction gives **smaller p-values than the reference alone**. If p-values do "
      "not improve, the added samples are not carrying the effect.")
    w()

    rows = []
    for name, _, _ in PAIRS:
        R, M, _, _ = refs[name]
        tR, tM = maps[R], maps[M]
        for t in POS + NEG:
            if t not in tR or t not in tM:
                continue
            cR, cM = tR[t], tM[t]
            hR, hM = d[cR].notna(), d[cM].notna()
            nR = d[hR].groupby("cell").size()
            nM = d[hM].groupby("cell").size()
            cells = sorted(set(nR[nR >= MIN_CELL].index) & set(nM[nM >= MIN_CELL].index))
            if len(cells) < 2:
                continue
            gR, gM = d[hR].groupby("cell")[cR], d[hM].groupby("cell")[cM]
            wt = np.minimum(nR[cells], nM[cells]).astype(float)
            MR = np.average(gR.median().reindex(cells), weights=wt)
            MM = np.average(gM.median().reindex(cells), weights=wt)
            IR = np.average((gR.quantile(.75) - gR.quantile(.25)).reindex(cells), weights=wt)
            IM = np.average((gM.quantile(.75) - gM.quantile(.25)).reindex(cells), weights=wt)
            slope = IR / IM if IM > 1e-9 else 1.0
            inter = MR - slope * MM

            dd = d.copy()
            dd["_ref"] = dd[cR]
            dd["_map"] = dd[cM]
            dd["_harm"] = dd[cR].where(dd[cR].notna(), slope * dd[cM] + inter)

            out = {k: lmm(dd, c) for k, c in
                   (("reference", "_ref"), ("mapped", "_map"), ("harmonized", "_harm"))}
            if out["reference"] is None or out["harmonized"] is None:
                continue
            rec = dict(block=name, analyte=t, control="neg" if t in NEG else "pos",
                       reference_project=R, mapped_project=M,
                       slope=slope, intercept=inter)
            for k, v in out.items():
                for s in ("n_obs", "n_ppl", "obs_per_ppl", "beta", "se", "p"):
                    rec[f"{k}_{s}"] = v[s] if v else np.nan
            rows.append(rec)
            print(f"  {name:19s}{t:8s} ref P={rec['reference_p']:.2e} -> "
                  f"harm P={rec['harmonized_p']:.2e}")

    R_ = pd.DataFrame(rows)
    R_["z_het"] = ((R_.reference_beta - R_.mapped_beta)
                   / np.sqrt(R_.reference_se ** 2 + R_.mapped_se ** 2))
    R_["p_het"] = 2 * stats.norm.sf(R_.z_het.abs())
    R_["z_shift"] = ((R_.harmonized_beta - R_.reference_beta)
                     / np.sqrt(R_.harmonized_se ** 2 + R_.reference_se ** 2))
    R_.to_csv(os.path.join(HERE, "harmonization_positive_controls.tab"),
              sep="\t", index=False)

    for kind, label in (("pos", "Positive controls"), ("neg", "Negative controls")):
        S = R_[R_.control == kind]
        if S.empty:
            continue
        w(f"## {label}")
        w()
        w("| block | analyte | slope | **ref** obs / ppl | ref beta | ref P | "
          "**harm** obs / ppl | harm beta | harm P | P better? | mapped beta | mapped P |")
        w("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for _, r in S.sort_values(["block", "analyte"]).iterrows():
            better = "**yes**" if r.harmonized_p < r.reference_p else "no"
            w(f"| `{r.block}` | {r.analyte} | {r.slope:.2f} "
              f"| {r.reference_n_obs:.0f} / {r.reference_n_ppl:.0f} "
              f"| {r.reference_beta:+.3f} | {r.reference_p:.1e} "
              f"| {r.harmonized_n_obs:.0f} / {r.harmonized_n_ppl:.0f} "
              f"| {r.harmonized_beta:+.3f} | {r.harmonized_p:.1e} | {better} "
              f"| {r.mapped_beta:+.3f} | {r.mapped_p:.1e} |")
        w()

    P = R_[R_.control == "pos"]
    if not P.empty:
        sig = P[P.reference_p < 0.05]
        w("## Summary")
        w()
        w(f"- positive-control fits: **{len(P)}**")
        w(f"- observations: reference median {P.reference_n_obs.median():.0f} -> "
          f"harmonized {P.harmonized_n_obs.median():.0f} "
          f"(**+{P.harmonized_n_obs.median() - P.reference_n_obs.median():.0f}**)")
        w(f"- participants: reference median {P.reference_n_ppl.median():.0f} -> "
          f"harmonized {P.harmonized_n_ppl.median():.0f}")
        w(f"- median SE: {P.reference_se.median():.4f} -> **{P.harmonized_se.median():.4f}** "
          f"({(1 - P.harmonized_se.median() / P.reference_se.median()) * 100:.0f}% smaller)")
        w(f"- **p-value improved in {(P.harmonized_p < P.reference_p).sum()} of {len(P)} "
          "fits overall**")
        if len(sig):
            w(f"- among the {len(sig)} where the reference was already significant, "
              f"improved in **{(sig.harmonized_p < sig.reference_p).sum()}**")
        w(f"- beta shift beyond sampling noise (|z| > 1.96): "
          f"**{(P.z_shift.abs() > 1.96).sum()} of {len(P)}**")
        w(f"- constituent projects differ (P_het < 0.05): "
          f"**{(P.p_het < 0.05).sum()} of {len(P)}** "
          f"({0.05 * len(P):.1f} expected by chance)")
        w()

    with open(REPORT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\nwrote {REPORT}")


if __name__ == "__main__":
    sys.exit(main())
