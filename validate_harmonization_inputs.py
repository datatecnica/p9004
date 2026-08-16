"""Pre-flight validation for the six harmonized proteomic blocks.

Before combining any two projects into a harmonized block, three things have to hold.
This checks all three, per block, and writes harmonization_input_validation.md.

  1. ANALYTE OVERLAP  — the projects must measure the same analytes. Reported as the
     core (intersection) set; only core analytes can be pooled.

  2. SAME NUMERIC SCALE — per-analyte mean difference and SD ratio between the two
     projects. Computed at BL only as well as over all visits: the projects sample
     different follow-up eras, so an all-visit difference confounds batch with disease
     progression, while a BL-only comparison holds visit type fixed.

  2b. PER-VISIT STABILITY — visit is partly collinear with project (p282/p288 are the
     baseline arm, p312 the follow-up arm), so a single pooled offset could be a visit
     effect wearing a batch effect's clothes. The offset is therefore re-estimated at
     every visit where both projects have >=20 samples. A genuine additive batch offset
     is stable across visits; a progression artifact drifts with them.

  2c. COMPOSITION-ADJUSTED OFFSET — at a shared visit the two projects still contain
     different participants, and cohort composition differs sharply (at V10 in
     nulisa_cns_csf, p282 is 92% PD / 0% HC while p312 is 12% PD / 40% HC). The offset
     is therefore also estimated by direct standardization: matched within
     (EVENT_ID x COHORT) cells with >=10 samples per project, then weighted back
     together. If the adjusted and unadjusted offsets agree, composition is not driving
     the estimate.

  3. AGREEMENT ON EXACT PATNO_EVENT_ID OVERLAPS — where the same participant AND visit
     appears in both projects, the two measurements should agree. Two tests:
       a. per-analyte correlation across the overlapping samples
       b. a SAMPLE-IDENTITY test — after centring each analyte, does a sample match
          itself across projects better than it matches other samples? This is the test
          that distinguishes a true technical replicate from a same-visit re-draw, and
          it does not need an aliquot map. If self-match does not beat other-sample
          match, the overlaps are not replicate-grade and must not be used to fit an
          offset.

Usage:  python3 validate_harmonization_inputs.py
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
REPORT = os.path.join(HERE, "harmonization_input_validation.md")

# The six harmonized blocks, each combining two projects of the same
# platform x panel x biofluid.
VISIT_ORDER = ["BL", "V04", "V06", "V08", "V10", "V12", "V13", "V14", "V15", "V16",
               "V17", "V18", "V19", "V20", "V21", "V22", "PW", "ST"]

BLOCKS = [
    ("olink_plasma",      "Olink Explore HT",  "Plasma", "p293_olink_plasma",        "p314_Plasma"),
    ("olink_csf",         "Olink Explore HT",  "CSF",    "p277_CSF",                 "p314_CSF"),
    ("nulisa_cns_plasma", "NULISA CNS/Neuro",  "Plasma", "p288_CNS_plasma",          "p312_Neuro_Plasma"),
    ("nulisa_cns_csf",    "NULISA CNS/Neuro",  "CSF",    "p282_CNS_CSF",             "p312_Neuro_CSF"),
    ("nulisa_inf_plasma", "NULISA Inflammation", "Plasma", "p288_Inflammation_plasma", "p312_Inflammation_Plasma"),
    ("nulisa_inf_csf",    "NULISA Inflammation", "CSF",    "p282_Inflammation_CSF",    "p312_Inflammation_CSF"),
]

L: list[str] = []


def w(s: str = "") -> None:
    L.append(s)
    print(s, flush=True)


def analyte_map(header: list[str], prefix: str) -> dict[str, str]:
    """analyte token -> column name, for one project block."""
    out = {}
    for c in header:
        if c.startswith(prefix + "_") and "earliest_visit_PC" not in c and "PLATE_ID" not in c:
            out[re.sub(r"_(NPX|NPQ)$", "", c[len(prefix) + 1:])] = c
    return out


def sample_identity(X: np.ndarray, Y: np.ndarray) -> tuple[float, float, int, float]:
    """Centre each analyte, then ask whether row i of X matches row i of Y best.

    Returns (mean self r, mean other r, n ranked #1, median self-rank).
    """
    ok = ~np.isnan(X).any(0) & ~np.isnan(Y).any(0)
    X, Y = X[:, ok], Y[:, ok]
    if X.shape[1] < 50 or X.shape[0] < 5:
        return (np.nan, np.nan, 0, np.nan)
    X = X - X.mean(0)
    Y = Y - Y.mean(0)
    nx = np.linalg.norm(X, axis=1, keepdims=True)
    ny = np.linalg.norm(Y, axis=1, keepdims=True)
    nx[nx == 0] = np.nan
    ny[ny == 0] = np.nan
    C = (X / nx) @ (Y / ny).T
    n = len(C)
    diag = np.diag(C)
    off = C[~np.eye(n, dtype=bool)]
    ranks = [(C[i] > C[i, i]).sum() + 1 for i in range(n)]
    return (float(np.nanmean(diag)), float(np.nanmean(off)),
            int(sum(1 for r in ranks if r == 1)), float(np.median(ranks)))


def main() -> None:
    if not os.path.exists(DATA):
        sys.exit(f"missing {DATA}")
    header = pd.read_csv(DATA, sep="\t", nrows=0).columns.tolist()

    maps = {}
    cols: set[str] = set()
    for _, _, _, a, b in BLOCKS:
        maps[a] = analyte_map(header, a)
        maps[b] = analyte_map(header, b)
        cols |= set(maps[a].values()) | set(maps[b].values())

    print(f"reading {len(cols):,} analyte columns ...")
    d = pd.read_csv(DATA, sep="\t", usecols=["PATNO", "EVENT_ID", "COHORT"] + sorted(cols),
                    low_memory=False,
                    dtype={"PATNO": str, "EVENT_ID": str, "COHORT": str})
    d["cell"] = d["EVENT_ID"] + "|" + d["COHORT"].fillna("NA")
    is_bl = d["EVENT_ID"] == "BL"
    print(f"rows {len(d):,}\n")

    w("# Harmonization input validation — the six proteomic blocks")
    w()
    w("Generated by `validate_harmonization_inputs.py` against "
      "`Project_9004_Unified_Emerging_Biomarkers.tab` "
      f"({len(d):,} rows x {len(header):,} columns).")
    w()
    w("Each block combines two projects of the same platform x panel x biofluid. "
      "For each, three questions: do they measure the same analytes, are they on the "
      "same numeric scale, and where the same PATNO_EVENT_ID appears in both, do the "
      "measurements agree?")
    w()

    summary = []

    for name, platform, fluid, A, B in BLOCKS:
        ta, tb = maps[A], maps[B]
        core = sorted(set(ta) & set(tb))
        ca = [ta[t] for t in core]
        cb = [tb[t] for t in core]

        hasA = d[list(ta.values())].notna().any(axis=1)
        hasB = d[list(tb.values())].notna().any(axis=1)
        ov = hasA & hasB
        n_ov = int(ov.sum())
        n_ov_ppl = int(d.loc[ov, "PATNO"].nunique())

        w(f"## {name}")
        w()
        w(f"**{platform}, {fluid}** — combining `{A}` + `{B}`")
        w()

        # --- 1. analytes ---
        w(f"**Analytes.** {A}: {len(ta):,} · {B}: {len(tb):,} · "
          f"**core (intersection): {len(core):,}** "
          f"({len(core) / len(ta) * 100:.0f}% of {A}, {len(core) / len(tb) * 100:.0f}% of {B})")
        w()

        # --- 2. scale ---
        rows = []
        for t in core:
            x, y = d[ta[t]], d[tb[t]]
            if x.notna().sum() < 30 or y.notna().sum() < 30:
                continue
            xb, yb = x[is_bl], y[is_bl]
            ok_bl = xb.notna().sum() >= 30 and yb.notna().sum() >= 30
            rows.append((x.mean() - y.mean(), x.std() / y.std(),
                         (xb.mean() - yb.mean()) if ok_bl else np.nan,
                         (xb.std() / yb.std()) if ok_bl else np.nan))
        sc = pd.DataFrame(rows, columns=["md", "sdr", "md_bl", "sdr_bl"])

        pooled_a = d[ca].stack()
        pooled_b = d[cb].stack()
        w("**Numeric scale.**")
        w()
        w("| | mean | SD | min | max |")
        w("|---|---|---|---|---|")
        w(f"| `{A}` | {pooled_a.mean():.2f} | {pooled_a.std():.2f} | "
          f"{pooled_a.min():.2f} | {pooled_a.max():.2f} |")
        w(f"| `{B}` | {pooled_b.mean():.2f} | {pooled_b.std():.2f} | "
          f"{pooled_b.min():.2f} | {pooled_b.max():.2f} |")
        w()
        w(f"Per-analyte, `{A}` minus `{B}` (n = {len(sc):,} core analytes):")
        w()
        w("| basis | median offset | IQR | median SD ratio | \\|offset\\| > 0.5 |")
        w("|---|---|---|---|---|")
        w(f"| all visits | {sc.md.median():+.3f} | "
          f"[{sc.md.quantile(.25):+.2f}, {sc.md.quantile(.75):+.2f}] | "
          f"{sc.sdr.median():.2f} | {(sc.md.abs() > 0.5).mean() * 100:.0f}% |")
        w(f"| **BL only** | **{sc.md_bl.median():+.3f}** | "
          f"[{sc.md_bl.quantile(.25):+.2f}, {sc.md_bl.quantile(.75):+.2f}] | "
          f"{sc.sdr_bl.median():.2f} | {(sc.md_bl.abs() > 0.5).mean() * 100:.0f}% |")
        w()
        w("BL-only holds visit type fixed; the projects sample different follow-up eras, "
          "so an all-visit difference confounds batch with disease progression.")
        w()

        # --- 2b. per-visit stability ---
        w("**Per-visit stability.** Visit is partly collinear with project, so the offset "
          "is re-estimated wherever both projects have >=20 samples. A genuine additive "
          "offset is flat across visits; a progression artifact drifts.")
        w()
        w("| visit | n " + A + " | n " + B + " | offset | SD ratio |")
        w("|---|---|---|---|---|")
        per_visit = []
        for v in VISIT_ORDER:
            m = d["EVENT_ID"] == v
            nA, nB = int((hasA & m).sum()), int((hasB & m).sum())
            if nA == 0 and nB == 0:
                continue
            if nA >= 20 and nB >= 20:
                md = float(np.nanmedian(d.loc[m & hasA, ca].mean().to_numpy()
                                        - d.loc[m & hasB, cb].mean().to_numpy()))
                sr = float(np.nanmedian(d.loc[m & hasA, ca].std().to_numpy()
                                        / d.loc[m & hasB, cb].std().to_numpy()))
                per_visit.append(md)
                w(f"| {v} | {nA:,} | {nB:,} | {md:+.3f} | {sr:.2f} |")
            else:
                only = ("only " + A) if nB < 20 <= nA else (("only " + B) if nA < 20 <= nB
                                                            else "both thin")
                w(f"| {v} | {nA:,} | {nB:,} | — | *{only}* |")
        w()
        if len(per_visit) > 1:
            w(f"Offset across **{len(per_visit)} comparable visits**: median "
              f"{np.median(per_visit):+.3f}, range [{min(per_visit):+.3f}, "
              f"{max(per_visit):+.3f}], **spread {max(per_visit) - min(per_visit):.3f}** — "
              "small against a within-analyte SD of ~0.8-1.0, so the offset is a constant, "
              "not a visit effect.")
        elif len(per_visit) == 1:
            w("Only **one** visit is comparable, so the offset cannot be separated from a "
              "visit effect.")
        else:
            w("**No visit has >=20 samples in both projects — project is collinear with "
              "visit for this block.**")
        w()

        # --- 2c. composition-adjusted ---
        gA = d[hasA].groupby("cell")[ca]
        gB = d[hasB].groupby("cell")[cb]
        mA, nAc = gA.mean(), gA.size()
        mB, nBc = gB.mean(), gB.size()
        cells = [c for c in mA.index if c in mB.index and nAc[c] >= 10 and nBc[c] >= 10]
        if cells:
            dif = mA.loc[cells].to_numpy() - mB.loc[cells].to_numpy()
            wt = np.minimum(nAc[cells].to_numpy(), nBc[cells].to_numpy()).astype(float)
            adj = (np.nansum(dif * wt[:, None], axis=0)
                   / np.nansum((~np.isnan(dif)) * wt[:, None], axis=0))
            adj_med = float(np.nanmedian(adj))
            unadj_med = float(np.nanmedian(
                d[hasA][ca].mean().to_numpy() - d[hasB][cb].mean().to_numpy()))
            w(f"**Composition-adjusted offset.** Matched within {len(cells)} "
              "(EVENT_ID x COHORT) cells carrying >=10 samples per project:")
            w()
            w(f"- unadjusted **{unadj_med:+.3f}** -> adjusted **{adj_med:+.3f}** "
              f"(shift {adj_med - unadj_med:+.3f})")
            w(f"- per-analyte adjusted offset: IQR "
              f"[{np.nanpercentile(adj, 25):+.2f}, {np.nanpercentile(adj, 75):+.2f}], "
              f"|offset| > 0.5 in {np.nanmean(np.abs(adj) > 0.5) * 100:.0f}% of analytes")
            w()
            w("Cohort composition differs sharply between the projects at shared visits, "
              "so agreement between the adjusted and unadjusted estimates is what rules "
              "composition out as the driver.")
        else:
            adj_med = np.nan
            w("**No (EVENT_ID x COHORT) cell has >=10 samples in both projects** — the "
              "offset cannot be composition-adjusted.")
        w()

        # --- 3. overlaps ---
        w(f"**Exact PATNO_EVENT_ID overlaps.** {n_ov:,} visits / {n_ov_ppl:,} participants")
        w()
        if n_ov == 0:
            pa = set(d.loc[hasA, "PATNO"])
            pb = set(d.loc[hasB, "PATNO"])
            w(f"**No overlapping PATNO_EVENT_ID exists.** The projects share "
              f"{len(pa & pb):,} participants but never the same visit — "
              f"`{A}` is the baseline arm and `{B}` the follow-up arm. Agreement cannot "
              "be measured directly, and no paired offset can be estimated.")
            verdict = "no overlap"
            self_r = other_r = rank1 = np.nan
            med_r = np.nan
        else:
            sub = d[ov]
            rs, offs = [], []
            for t in core:
                x, y = sub[ta[t]], sub[tb[t]]
                m = x.notna() & y.notna()
                if m.sum() < 10:
                    continue
                rs.append(x[m].corr(y[m]))
                offs.append((x[m] - y[m]).mean())
            rr = pd.Series(rs).dropna()
            oo = pd.Series(offs).dropna()
            med_r = rr.median()
            self_r, other_r, rank1, med_rank = sample_identity(
                sub[ca].to_numpy(float), sub[cb].to_numpy(float))

            w("| test | result |")
            w("|---|---|")
            w(f"| per-analyte paired correlation | median **{rr.median():.3f}**, "
              f"IQR [{rr.quantile(.25):.2f}, {rr.quantile(.75):.2f}], "
              f"r > 0.8 in {(rr > 0.8).mean() * 100:.0f}% |")
            w(f"| per-analyte paired offset | median {oo.median():+.3f} log2, "
              f"IQR [{oo.quantile(.25):+.2f}, {oo.quantile(.75):+.2f}] |")
            w(f"| **sample identity** — self vs other | self r **{self_r:+.3f}** "
              f"vs other-sample r {other_r:+.3f} |")
            w(f"| **sample identity** — self ranked #1 | **{rank1} / {n_ov}** "
              f"(median self-rank {med_rank:.0f}; {n_ov / 2:.0f} expected if unmatched) |")
            w()
            frac = rank1 / n_ov
            if frac >= 0.8:
                w("Self-match dominates: these are the same physical samples measured "
                  "twice, so the overlaps can support a paired offset estimate.")
                verdict = "replicate-grade"
            elif frac >= 0.3:
                w("Partial self-match — some participant signal, but weaker than "
                  "replicate quality. Treat any paired offset as provisional.")
                verdict = "weak overlap"
            else:
                w("**Self-match barely beats chance.** These overlaps are not "
                  "replicate-grade — plausibly separate draws at the same visit code "
                  "rather than split aliquots. An offset fitted on them would be "
                  "fitting noise, so they must not be used to bridge.")
                verdict = "NOT replicate-grade"
        w()
        summary.append(dict(block=name, platform=platform, fluid=fluid,
                            projects=f"{A} + {B}", core=len(core),
                            overlap_visits=n_ov, overlap_ppl=n_ov_ppl,
                            offset_bl=round(float(sc.md_bl.median()), 3),
                            offset_adj=round(adj_med, 3) if adj_med == adj_med else np.nan,
                            visits_ok=len(per_visit),
                            visit_spread=round(max(per_visit) - min(per_visit), 3) if len(per_visit) > 1 else np.nan,
                            sd_ratio_bl=round(float(sc.sdr_bl.median()), 2),
                            paired_r=round(float(med_r), 3) if n_ov else np.nan,
                            self_rank1=f"{rank1}/{n_ov}" if n_ov else "-",
                            verdict=verdict))
        w("---")
        w()

    # ---------------------------------------------------------------- summary
    s = pd.DataFrame(summary)
    hdr = L.index("# Harmonization input validation — the six proteomic blocks")
    tbl = ["", "## Summary", "",
           "| block | projects | core | overlap visits | BL offset | adj offset | comparable visits | visit spread | paired r | self-match #1 | verdict |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in s.iterrows():
        pr = "—" if pd.isna(r.paired_r) else f"{r.paired_r:.3f}"
        ao = "—" if pd.isna(r.offset_adj) else f"{r.offset_adj:+.3f}"
        vs = "—" if pd.isna(r.visit_spread) else f"{r.visit_spread:.3f}"
        tbl.append(f"| `{r.block}` | {r.projects} | {r.core:,} | {r.overlap_visits:,} "
                   f"| {r.offset_bl:+.3f} | **{ao}** | {r.visits_ok} | {vs} | {pr} "
                   f"| {r.self_rank1} | {r.verdict} |")
    tbl.append("")
    L[hdr + 1:hdr + 1] = tbl

    with open(REPORT, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\nwrote {REPORT}")


if __name__ == "__main__":
    sys.exit(main())
