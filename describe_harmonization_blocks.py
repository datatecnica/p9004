"""Descriptive statistics for the six harmonized proteomic blocks.

One table per question, covering every pair that will be combined:

  1. composition        — which projects, how many analytes, rows, participants
  2. descriptives       — mean / SD / median / Q1 / Q3 / IQR / min / max per project
  3. pair comparison    — location and dispersion differences, and the correction implied
  4. overlapping pairs  — samples appearing at the same PATNO_EVENT_ID in both projects,
                          with N, r and r-squared. Reported both over all analytes and
                          restricted to analytes with SD >= 1.0 in BOTH projects, because
                          low-variance analytes attenuate correlation and hide whether
                          the overlaps are genuine replicates.

Writes harmonization_block_descriptives.md.
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np
import pandas as pd

from build_common import DATASET_STEM, require_build

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = require_build(DATASET_STEM)
REPORT = os.path.join(HERE, "harmonization_block_descriptives.md")

BLOCKS = [
    ("olink_plasma",      "Olink Explore HT",    "Plasma", "NPX", "p293_olink_plasma",        "p314_Plasma"),
    ("olink_csf",         "Olink Explore HT",    "CSF",    "NPX", "p277_CSF",                 "p314_CSF"),
    ("nulisa_cns_plasma", "NULISA CNS/Neuro",    "Plasma", "NPQ", "p288_CNS_plasma",          "p312_Neuro_Plasma"),
    ("nulisa_cns_csf",    "NULISA CNS/Neuro",    "CSF",    "NPQ", "p282_CNS_CSF",             "p312_Neuro_CSF"),
    ("nulisa_inf_plasma", "NULISA Inflammation", "Plasma", "NPQ", "p288_Inflammation_plasma", "p312_Inflammation_Plasma"),
    ("nulisa_inf_csf",    "NULISA Inflammation", "CSF",    "NPQ", "p282_Inflammation_CSF",    "p312_Inflammation_CSF"),
]

SD_FLOOR = 1.0   # "sufficient variability" threshold for the restricted r-squared
L: list[str] = []


def w(s: str = "") -> None:
    L.append(s)
    print(s, flush=True)


def amap(header: list[str], prefix: str) -> dict[str, str]:
    return {re.sub(r"_(NPX|NPQ)$", "", c[len(prefix) + 1:]): c
            for c in header
            if c.startswith(prefix + "_") and "earliest_visit_PC" not in c
            and "PLATE_ID" not in c}


def deming(x: np.ndarray, y: np.ndarray, lam: float = 1.0) -> float:
    xm, ym = x.mean(), y.mean()
    Sxx = ((x - xm) ** 2).sum()
    Syy = ((y - ym) ** 2).sum()
    Sxy = ((x - xm) * (y - ym)).sum()
    if Sxy == 0:
        return np.nan
    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2 + 4 * lam * Sxy ** 2))
                 / (2 * Sxy))


def main() -> None:
    header = pd.read_csv(DATA, sep="\t", nrows=0).columns.tolist()
    maps = {}
    cols: set[str] = set()
    for *_, a, b in BLOCKS:
        maps[a] = amap(header, a)
        maps[b] = amap(header, b)
        cols |= set(maps[a].values()) | set(maps[b].values())

    print(f"reading {len(cols):,} analyte columns ...")
    d = pd.read_csv(DATA, sep="\t", usecols=["PATNO", "EVENT_ID"] + sorted(cols),
                    low_memory=False, dtype={"PATNO": str, "EVENT_ID": str})
    print(f"rows {len(d):,}\n")

    comp, desc, cmpr, ovl = [], [], [], []

    for name, plat, fluid, metric, A, B in BLOCKS:
        ta, tb = maps[A], maps[B]
        core = sorted(set(ta) & set(tb))
        ca = [ta[t] for t in core]
        cb = [tb[t] for t in core]
        hA = d[list(ta.values())].notna().any(axis=1)
        hB = d[list(tb.values())].notna().any(axis=1)
        ov = hA & hB

        comp.append(dict(block=name, platform=plat, fluid=fluid, metric=metric,
                         project_A=A, project_B=B,
                         analytes_A=len(ta), analytes_B=len(tb), core=len(core),
                         rows_A=int(hA.sum()), rows_B=int(hB.sum()),
                         ppl_A=int(d.loc[hA, "PATNO"].nunique()),
                         ppl_B=int(d.loc[hB, "PATNO"].nunique())))

        stats = {}
        for proj, cc, mask in ((A, ca, hA), (B, cb, hB)):
            v = d.loc[mask, cc].stack()
            q1, med, q3 = v.quantile([.25, .5, .75])
            stats[proj] = dict(n=len(v), mean=v.mean(), sd=v.std(), median=med,
                               q1=q1, q3=q3, iqr=q3 - q1, mn=v.min(), mx=v.max())
            desc.append(dict(block=name, project=proj, **stats[proj]))

        sA, sB = stats[A], stats[B]
        pa_med = d.loc[hA, ca].median().to_numpy() - d.loc[hB, cb].median().to_numpy()
        iA = (d.loc[hA, ca].quantile(.75) - d.loc[hA, ca].quantile(.25)).to_numpy()
        iB = (d.loc[hB, cb].quantile(.75) - d.loc[hB, cb].quantile(.25)).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            iqr_ratio = iA / iB
        loc_shift = float(np.nanmedian(pa_med))
        disp = float(np.nanmedian(iqr_ratio))
        needs_slope = abs(disp - 1) >= 0.10
        cmpr.append(dict(block=name,
                         mean_diff=sA["mean"] - sB["mean"],
                         median_diff=sA["median"] - sB["median"],
                         sd_ratio=sA["sd"] / sB["sd"],
                         iqr_ratio_pooled=sA["iqr"] / sB["iqr"],
                         analyte_median_diff=loc_shift,
                         analyte_iqr_ratio=disp,
                         iqr_ratio_q1=float(np.nanpercentile(iqr_ratio, 25)),
                         iqr_ratio_q3=float(np.nanpercentile(iqr_ratio, 75)),
                         correction="intercept + slope" if needs_slope else "intercept only"))

        # --- overlapping sample/visit pairs ---
        n_ov = int(ov.sum())
        row = dict(block=name, overlap_visits=n_ov,
                   overlap_ppl=int(d.loc[ov, "PATNO"].nunique()))
        if n_ov >= 10:
            sub = d[ov]
            rec = []
            for t in core:
                x, y = sub[ta[t]], sub[tb[t]]
                m = x.notna() & y.notna()
                if m.sum() < 10:
                    continue
                xv, yv = x[m].to_numpy(float), y[m].to_numpy(float)
                r = np.corrcoef(xv, yv)[0, 1]
                rec.append((min(xv.std(ddof=1), yv.std(ddof=1)), r, deming(xv, yv)))
            R = pd.DataFrame(rec, columns=["minsd", "r", "slope"]).dropna(subset=["r"])
            hi = R[R.minsd >= SD_FLOOR]
            row.update(n_analytes=len(R), median_r=R.r.median(), median_r2=(R.r ** 2).median(),
                       n_hi=len(hi),
                       median_r_hi=hi.r.median() if len(hi) else np.nan,
                       median_r2_hi=(hi.r ** 2).median() if len(hi) else np.nan,
                       pct_r2_gt_half_hi=(hi.r ** 2 > 0.5).mean() * 100 if len(hi) else np.nan,
                       deming_hi=hi.slope.median() if len(hi) else np.nan)
        else:
            row.update(n_analytes=0, median_r=np.nan, median_r2=np.nan, n_hi=0,
                       median_r_hi=np.nan, median_r2_hi=np.nan,
                       pct_r2_gt_half_hi=np.nan, deming_hi=np.nan)
        ovl.append(row)

    C, D_, M, O = (pd.DataFrame(x) for x in (comp, desc, cmpr, ovl))

    w("# Harmonized block descriptives — the six pairs")
    w()
    w(f"Generated by `describe_harmonization_blocks.py` against the rebuilt dataset "
      f"({len(d):,} rows). Each block combines two projects of the same "
      "platform x panel x biofluid. Values are on native NPQ/NPX scales throughout.")
    w()

    w("## 1. Block composition")
    w()
    w("| block | platform | fluid | metric | project A | project B | analytes A | analytes B | **core** | rows A | rows B | participants A | participants B |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in C.iterrows():
        w(f"| `{r.block}` | {r.platform} | {r.fluid} | {r.metric} | {r.project_A} | {r.project_B} "
          f"| {r.analytes_A:,} | {r.analytes_B:,} | **{r.core:,}** | {r.rows_A:,} | {r.rows_B:,} "
          f"| {r.ppl_A:,} | {r.ppl_B:,} |")
    w()

    w("## 2. Descriptive statistics per project")
    w()
    w("Pooled over all core analytes x all samples with data.")
    w()
    w("| block | project | n values | mean | SD | median | Q1 | Q3 | IQR | min | max |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in D_.iterrows():
        w(f"| `{r.block}` | {r.project} | {r.n:,} | {r['mean']:.2f} | {r.sd:.2f} | "
          f"{r['median']:.2f} | {r.q1:.2f} | {r.q3:.2f} | {r.iqr:.2f} | {r.mn:.2f} | {r.mx:.2f} |")
    w()

    w("## 3. Pair comparison and implied correction")
    w()
    w("`analyte median diff` and `analyte IQR ratio` are the medians of the per-analyte "
      "statistics — the robust analogues of offset and scale factor, and what a "
      "correction would actually be fitted on. A slope is flagged when the median "
      "per-analyte IQR ratio departs from 1 by >=10%.")
    w()
    w("| block | mean diff | median diff | SD ratio | IQR ratio (pooled) | **analyte median diff** | **analyte IQR ratio** | IQR ratio IQR | correction |")
    w("|---|---|---|---|---|---|---|---|---|")
    for _, r in M.iterrows():
        w(f"| `{r.block}` | {r.mean_diff:+.3f} | {r.median_diff:+.3f} | {r.sd_ratio:.2f} "
          f"| {r.iqr_ratio_pooled:.2f} | **{r.analyte_median_diff:+.3f}** "
          f"| **{r.analyte_iqr_ratio:.2f}** | [{r.iqr_ratio_q1:.2f}, {r.iqr_ratio_q3:.2f}] "
          f"| {r.correction} |")
    w()

    w("## 4. Samples present at the same PATNO_EVENT_ID in both projects")
    w()
    w(f"`r2 (all)` uses every core analyte; `r2 (SD>={SD_FLOOR})` restricts to analytes "
      f"with SD >= {SD_FLOOR} in **both** projects. Low-variance analytes cannot "
      "correlate, so the restricted figure is the one that says whether the overlapping "
      "samples are genuine replicates: under attenuation the restricted r2 rises, and "
      "if it does not, the pairs are not the same material.")
    w()
    w("| block | **overlapping visits** | participants | analytes | median r | **median r2** | analytes SD>=1 | median r (SD>=1) | **median r2 (SD>=1)** | r2>0.5 | Deming slope (SD>=1) |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in O.iterrows():
        if r.overlap_visits < 10:
            w(f"| `{r.block}` | **{r.overlap_visits}** | {r.overlap_ppl} | — | — | — "
              "| — | — | — | — | — |")
            continue
        w(f"| `{r.block}` | **{r.overlap_visits}** | {r.overlap_ppl} | {r.n_analytes:,} "
          f"| {r.median_r:.3f} | **{r.median_r2:.3f}** | {r.n_hi:,} | {r.median_r_hi:.3f} "
          f"| **{r.median_r2_hi:.3f}** | {r.pct_r2_gt_half_hi:.0f}% | {r.deming_hi:.3f} |")
    w()
    w("Blocks showing 0 overlapping visits have none by design: p282/p288 are the "
      "baseline arm and p312 the follow-up arm, so they share participants but never a "
      "visit. Their corrections are fitted from visit- and cohort-matched groups instead.")
    w()

    with open(REPORT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\nwrote {REPORT}")


if __name__ == "__main__":
    sys.exit(main())
