"""Deep per-analyte comparison of each paired project, with proposed corrections.

For every core analyte in each of the six blocks this emits: gross and per-analyte
location/dispersion in both projects, the agreement statistics where the same
PATNO_EVENT_ID appears in both, a proposed linear correction, and a HELD-OUT check of
whether that correction actually works.

Correction form. In this diagnostic the older project is the reference and the newer one
is mapped onto it. NOTE this is the opposite orientation to build_step_9_harmonize.py for
the two Olink blocks: the build references whichever project has MORE ROWS, which is p314
for both, so slopes and intercepts here are reciprocal to the ones actually applied. The
four NULISA blocks coincide under either rule. The applied coefficients live in the data
dictionary; this table is a comparison diagnostic, not the source of the correction:

    B_corrected = slope * B + intercept

fitted robustly on (EVENT_ID x COHORT) matched cells rather than on the raw pooled
values, because cohort composition differs sharply between the projects at shared visits
(at V10 in nulisa_cns_csf, p282 is 92% PD / 0% HC against p312's 12% PD / 40% HC).
Within each matched cell the median and IQR of each project are taken, then combined
across cells weighted by the smaller of the two cell counts:

    slope     = IQR_A / IQR_B          (robust scale ratio)
    intercept = median_A - slope * median_B

Held-out validation. Cells are split into alternating folds by index; the correction is
fitted on fold 0 and the residual difference is evaluated on fold 1, and vice versa. A
correction that only reproduces its own training cells is not a correction.

Writes harmonization_corrections.tab (one row per analyte) and
harmonization_corrections.md (gross summary).
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
OUT_TAB = os.path.join(HERE, "harmonization_corrections.tab")
OUT_MD = os.path.join(HERE, "harmonization_corrections.md")

BLOCKS = [
    ("olink_plasma",      "Olink",  "Plasma", "p293_olink_plasma",        "p314_Plasma"),
    ("olink_csf",         "Olink",  "CSF",    "p277_CSF",                 "p314_CSF"),
    ("nulisa_cns_plasma", "NULISA", "Plasma", "p288_CNS_plasma",          "p312_Neuro_Plasma"),
    ("nulisa_cns_csf",    "NULISA", "CSF",    "p282_CNS_CSF",             "p312_Neuro_CSF"),
    ("nulisa_inf_plasma", "NULISA", "Plasma", "p288_Inflammation_plasma", "p312_Inflammation_Plasma"),
    ("nulisa_inf_csf",    "NULISA", "CSF",    "p282_Inflammation_CSF",    "p312_Inflammation_CSF"),
]

MIN_CELL = 10        # min samples per project within a matched cell
SD_FLOOR = 1.0       # "sufficient variability" for trusting a paired Deming slope
SLOPE_TOL = 0.10     # |slope-1| below this is treated as no scale difference
L: list[str] = []


def w(s: str = "") -> None:
    L.append(s)
    print(s, flush=True)


def amap(header, prefix):
    return {re.sub(r"_(NPX|NPQ)$", "", c[len(prefix) + 1:]): c
            for c in header
            if c.startswith(prefix + "_") and "earliest_visit_PC" not in c
            and "PLATE_ID" not in c}


def deming(x, y, lam=1.0):
    xm, ym = x.mean(), y.mean()
    Sxx = ((x - xm) ** 2).sum(); Syy = ((y - ym) ** 2).sum()
    Sxy = ((x - xm) * (y - ym)).sum()
    if Sxy == 0:
        return np.nan
    return float((Syy - lam * Sxx + np.sqrt((Syy - lam * Sxx) ** 2 + 4 * lam * Sxy ** 2))
                 / (2 * Sxy))


def cell_stats(d, cols, mask, cells):
    """median, IQR and n per (cell x analyte). Returns 3 arrays of shape (cells, analytes)."""
    sub = d.loc[mask, ["cell"] + cols]
    g = sub.groupby("cell")
    med = g[cols].median().reindex(cells)
    q1 = g[cols].quantile(.25).reindex(cells)
    q3 = g[cols].quantile(.75).reindex(cells)
    n = g.size().reindex(cells).fillna(0)
    return med.to_numpy(), (q3 - q1).to_numpy(), n.to_numpy()


def fit(medA, iqrA, medB, iqrB, wt):
    """Weighted robust slope/intercept mapping B onto A. Arrays are (cells, analytes)."""
    ok = np.isfinite(medA) & np.isfinite(medB) & np.isfinite(iqrA) & np.isfinite(iqrB)
    W = np.where(ok, wt[:, None], 0.0)
    tot = W.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mA = np.nansum(np.where(ok, medA, 0) * W, axis=0) / tot
        mB = np.nansum(np.where(ok, medB, 0) * W, axis=0) / tot
        iA = np.nansum(np.where(ok, iqrA, 0) * W, axis=0) / tot
        iB = np.nansum(np.where(ok, iqrB, 0) * W, axis=0) / tot
        slope = iA / iB
    slope = np.where(np.isfinite(slope) & (iB > 1e-9), slope, 1.0)
    intercept = mA - slope * mB
    return slope, intercept, mA, mB, iA, iB, tot


def main() -> None:
    header = pd.read_csv(DATA, sep="\t", nrows=0).columns.tolist()
    maps, cols = {}, set()
    for _, _, _, a, b in BLOCKS:
        maps[a] = amap(header, a); maps[b] = amap(header, b)
        cols |= set(maps[a].values()) | set(maps[b].values())
    print(f"reading {len(cols):,} analyte columns ...")
    d = pd.read_csv(DATA, sep="\t", usecols=["PATNO", "EVENT_ID", "COHORT"] + sorted(cols),
                    low_memory=False, dtype={"PATNO": str, "EVENT_ID": str, "COHORT": str})
    d["cell"] = d["EVENT_ID"] + "|" + d["COHORT"].fillna("NA")
    print(f"rows {len(d):,}\n")

    rows, gross = [], []

    for name, plat, fluid, A, B in BLOCKS:
        ta, tb = maps[A], maps[B]
        core = sorted(set(ta) & set(tb))
        ca, cb = [ta[t] for t in core], [tb[t] for t in core]
        hA = d[list(ta.values())].notna().any(axis=1)
        hB = d[list(tb.values())].notna().any(axis=1)
        ov = hA & hB

        # cells usable for fitting
        nA_cell = d.loc[hA].groupby("cell").size()
        nB_cell = d.loc[hB].groupby("cell").size()
        cells = sorted(set(nA_cell[nA_cell >= MIN_CELL].index)
                       & set(nB_cell[nB_cell >= MIN_CELL].index))
        print(f"{name}: {len(core):,} core analytes, {len(cells)} matched cells, "
              f"{int(ov.sum())} overlapping visits")

        medA, iqrA, nA = cell_stats(d, ca, hA, cells)
        medB, iqrB, nB = cell_stats(d, cb, hB, cells)
        wt = np.minimum(nA, nB).astype(float)

        slope, inter, mA, mB, iA, iB, tot = fit(medA, iqrA, medB, iqrB, wt)

        # held-out: fit on alternating cells, evaluate residual on the others
        resid = np.full(len(core), np.nan)
        if len(cells) >= 4:
            f0 = np.arange(len(cells)) % 2 == 0
            outs = []
            for tr in (f0, ~f0):
                te = ~tr
                s, i, *_ = fit(medA[tr], iqrA[tr], medB[tr], iqrB[tr], wt[tr])
                corrected = s[None, :] * medB[te] + i[None, :]
                outs.append(np.nanmedian(medA[te] - corrected, axis=0))
            resid = np.nanmean(np.vstack(outs), axis=0)
        # naive residual with no correction, same held-out cells, for comparison
        raw_gap = np.nanmedian(medA - medB, axis=0)

        # paired agreement where overlaps exist
        r = np.full(len(core), np.nan); r2 = np.full(len(core), np.nan)
        dslope = np.full(len(core), np.nan); n_ov = np.zeros(len(core), int)
        sdA_ov = np.full(len(core), np.nan); sdB_ov = np.full(len(core), np.nan)
        if ov.sum() >= 10:
            sub = d[ov]
            for j, t in enumerate(core):
                x, y = sub[ta[t]], sub[tb[t]]
                m = x.notna() & y.notna()
                n_ov[j] = int(m.sum())
                if m.sum() < 10:
                    continue
                xv, yv = x[m].to_numpy(float), y[m].to_numpy(float)
                sdA_ov[j], sdB_ov[j] = xv.std(ddof=1), yv.std(ddof=1)
                rr = np.corrcoef(xv, yv)[0, 1]
                r[j], r2[j] = rr, rr * rr
                dslope[j] = deming(xv, yv)

        applied_slope = np.where(np.abs(slope - 1) >= SLOPE_TOL, slope, 1.0)
        applied_inter = np.where(np.abs(slope - 1) >= SLOPE_TOL, inter, mA - mB)

        for j, t in enumerate(core):
            rows.append(dict(
                block=name, platform=plat, fluid=fluid, analyte=t,
                col_reference=ca[j], col_mapped=cb[j],
                n_ref=int(d.loc[hA, ca[j]].notna().sum()),
                n_map=int(d.loc[hB, cb[j]].notna().sum()),
                median_ref=np.round(mA[j], 4), median_map=np.round(mB[j], 4),
                iqr_ref=np.round(iA[j], 4), iqr_map=np.round(iB[j], 4),
                raw_median_gap=np.round(raw_gap[j], 4),
                fitted_slope=np.round(slope[j], 4),
                fitted_intercept=np.round(inter[j], 4),
                applied_slope=np.round(applied_slope[j], 4),
                applied_intercept=np.round(applied_inter[j], 4),
                heldout_residual=np.round(resid[j], 4),
                n_matched_cells=len(cells),
                n_overlap=n_ov[j],
                overlap_r=np.round(r[j], 4), overlap_r2=np.round(r2[j], 4),
                overlap_deming=np.round(dslope[j], 4),
                overlap_sd_ref=np.round(sdA_ov[j], 3), overlap_sd_map=np.round(sdB_ov[j], 3),
            ))

        gross.append(dict(
            block=name, projects=f"{A} -> {B}", core=len(core), cells=len(cells),
            overlap=int(ov.sum()),
            raw_gap=float(np.nanmedian(raw_gap)),
            slope_med=float(np.nanmedian(slope)),
            slope_q1=float(np.nanpercentile(slope, 25)),
            slope_q3=float(np.nanpercentile(slope, 75)),
            inter_med=float(np.nanmedian(inter)),
            pct_slope_flagged=float(np.mean(np.abs(slope - 1) >= SLOPE_TOL) * 100),
            resid_med=float(np.nanmedian(np.abs(resid))),
            raw_med_abs=float(np.nanmedian(np.abs(raw_gap))),
        ))

    T = pd.DataFrame(rows)
    T.to_csv(OUT_TAB, sep="\t", index=False)
    G = pd.DataFrame(gross)

    w("# Paired-project comparison and proposed corrections")
    w()
    w(f"Per-analyte table: `harmonization_corrections.tab` ({len(T):,} rows). "
      "In THIS diagnostic the older project is the reference and the newer is mapped "
      "onto it via `B_corrected = slope * B + intercept`, fitted on (EVENT_ID x COHORT) "
      "matched cells so cohort composition cannot drive the estimate.")
    w()
    w("**Orientation differs from the shipped build.** `build_step_9_harmonize.py` takes "
      "the project with MORE ROWS as the reference, so that most of each harmonized "
      "column is uncorrected native value. For the four NULISA blocks the older project "
      "is also the larger one and the two rules coincide. For `olink_plasma` and "
      "`olink_csf` they invert: this table treats p293 and p277 as the references, "
      "whereas the build maps them ONTO p314. Every number below, including the columns "
      "named `applied_slope` and `applied_intercept`, is therefore in this diagnostic's "
      "orientation and is NOT what the build applied to those two blocks — for "
      "`harmonized_olink_plasma_A1BG_NPX` this table gives 1.847 where the build applied "
      "0.541, its reciprocal. The coefficients actually applied are recorded per analyte "
      "in the data dictionary's derivation notes and are verified against the exported "
      "columns in `harmonized_build_validation.md`. The four NULISA blocks are unaffected.")
    w()
    w("## Gross comparison")
    w()
    w("| block | projects | core | matched cells | overlapping visits | raw median gap | fitted slope (median) | slope IQR | % analytes slope-corrected | intercept (median) |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for _, g in G.iterrows():
        w(f"| `{g.block}` | {g.projects} | {g.core:,} | {g.cells} | {g.overlap} "
          f"| {g.raw_gap:+.3f} | **{g.slope_med:.3f}** | [{g.slope_q1:.2f}, {g.slope_q3:.2f}] "
          f"| {g.pct_slope_flagged:.0f}% | {g.inter_med:+.3f} |")
    w()
    w("## Does the correction hold out?")
    w()
    w("Fitted on half the matched cells, residual measured on the other half. "
      "`|gap| before` is the same quantity with no correction applied.")
    w()
    w("| block | median \\|gap\\| before | median \\|residual\\| after | reduction |")
    w("|---|---|---|---|")
    for _, g in G.iterrows():
        red = (1 - g.resid_med / g.raw_med_abs) * 100 if g.raw_med_abs else np.nan
        w(f"| `{g.block}` | {g.raw_med_abs:.3f} | **{g.resid_med:.3f}** | {red:.0f}% |")
    w()
    w("## Per-analyte distributions")
    w()
    w("| block | slope <0.9 | slope 0.9-1.1 | slope >1.1 | \\|intercept\\| >0.5 | analytes with overlap r2 > 0.5 |")
    w("|---|---|---|---|---|---|")
    for name in T.block.unique():
        s = T[T.block == name]
        hi = s[(s.overlap_sd_ref >= SD_FLOOR) & (s.overlap_sd_map >= SD_FLOOR)]
        r2txt = (f"{(hi.overlap_r2 > 0.5).sum():,} / {len(hi):,}"
                 if len(hi) else "— no usable overlaps")
        w(f"| `{name}` | {(s.fitted_slope < 0.9).sum():,} | "
          f"{((s.fitted_slope >= 0.9) & (s.fitted_slope <= 1.1)).sum():,} | "
          f"{(s.fitted_slope > 1.1).sum():,} | {(s.applied_intercept.abs() > 0.5).sum():,} "
          f"| {r2txt} |")
    w()
    w(f"Columns in `harmonization_corrections.tab`: {', '.join(T.columns)}")
    w()

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"\nwrote {OUT_TAB}\nwrote {OUT_MD}")


if __name__ == "__main__":
    sys.exit(main())
