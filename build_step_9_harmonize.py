"""Step 9 — build the six harmonized proteomic blocks.

Each block combines the two projects that ran the same platform x panel x biofluid. The
reference is whichever project has MORE ROWS, so most of each harmonized column is
uncorrected native value. The mapped project is placed on the reference's scale:

    mapped_corrected = slope * mapped + intercept

fitted per analyte on (EVENT_ID x COHORT) matched cells with >=10 samples per project,
using medians and IQRs so cohort composition cannot drive the estimate:

    slope     = IQR_reference / IQR_mapped
    intercept = median_reference - slope * median_mapped

A slope is only applied where a bootstrap over cells shows the estimate is real —
32-62% of nominally-flagged slopes have CIs spanning the tolerance band, and the
upstream package's own validation shows spurious scale corrections at low cell counts
actively hurt. Slope defaults to 1 otherwise, leaving an intercept-only correction.

The reference value always wins; the corrected mapped value fills only rows the
reference did not cover. Values are never averaged, so each cell traces to one assay run.
Values stay on native NPQ/NPX scales — see "Why not z-score" in the README.

Also emits `collection_era`, the categorical collection period. p277's post-2020 samples
carry 2 healthy controls against 372 NSD+, which inflates its whole panel (lambda 2.49 vs
p314's 1.11). No samples are dropped: applying era as a covariate is an analysis-time
decision. See the README's p277 section.

Emits:
  build_intermediates/harmonized-<ts>.tab
  build_intermediates/harmonization_fit-<ts>.tab    per-analyte slope/intercept/support
"""

from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
import pandas as pd

from build_common import TIMESTAMP, log

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")

# (block, platform, fluid, true metric, project_1, project_2)
# The reference is chosen at runtime as whichever project has more rows.
BLOCKS = [
    ("olink_plasma",      "Olink Explore HT",    "Plasma", "NPX",
     "p293_olink_plasma",        "p314_Plasma"),
    ("olink_csf",         "Olink Explore HT",    "CSF",    "NPX",
     "p277_CSF",                 "p314_CSF"),
    ("nulisa_cns_plasma", "NULISA CNS/Neuro",    "Plasma", "NPQ",
     "p288_CNS_plasma",          "p312_Neuro_Plasma"),
    ("nulisa_cns_csf",    "NULISA CNS/Neuro",    "CSF",    "NPQ",
     "p282_CNS_CSF",             "p312_Neuro_CSF"),
    ("nulisa_inf_plasma", "NULISA Inflammation", "Plasma", "NPQ",
     "p288_Inflammation_plasma", "p312_Inflammation_Plasma"),
    ("nulisa_inf_csf",    "NULISA Inflammation", "CSF",    "NPQ",
     "p282_Inflammation_CSF",    "p312_Inflammation_CSF"),
]

MIN_CELL = 10        # min samples per project within a matched cell
SLOPE_TOL = 0.10     # |slope-1| below this is treated as no scale difference
N_BOOT = 400
ERA_CUT = 2020       # collection_era boundary; see the README's p277 section


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} — run the prior step first")
    return hits[-1]


def amap(header, prefix):
    return {re.sub(r"_(NPX|NPQ)$", "", c[len(prefix) + 1:]): c
            for c in header
            if c.startswith(prefix + "_") and "earliest_visit_PC" not in c
            and "PLATE_ID" not in c}


def main() -> None:
    log(f"=== Step 9: harmonization ===  timestamp {TIMESTAMP}")
    src = latest("merged-*.tab")
    header = pd.read_csv(src, sep="\t", nrows=0).columns.tolist()

    maps, cols = {}, set()
    for _, _, _, _, p1, p2 in BLOCKS:
        maps[p1], maps[p2] = amap(header, p1), amap(header, p2)
        cols |= set(maps[p1].values()) | set(maps[p2].values())

    log(f"reading {len(cols):,} analyte columns ...")
    d = pd.read_csv(src, sep="\t",
                    usecols=["key", "PATNO", "EVENT_ID", "COHORT", "visit_date"] + sorted(cols),
                    low_memory=False,
                    dtype={"key": str, "PATNO": str, "EVENT_ID": str, "COHORT": str})
    n_rows = len(d)
    d["cell"] = d["EVENT_ID"] + "|" + d["COHORT"].fillna("NA")
    log(f"rows {n_rows:,}")

    # --- collection_era -----------------------------------------------------
    yr = pd.to_datetime(d["visit_date"], errors="coerce").dt.year
    era = pd.Series(pd.NA, index=d.index, dtype="string")
    era[yr.notna() & (yr < ERA_CUT)] = f"pre_{ERA_CUT}"
    era[yr.notna() & (yr >= ERA_CUT)] = f"{ERA_CUT}_plus"
    log(f"\ncollection_era: {dict(era.value_counts(dropna=False))}")

    out = pd.DataFrame({"key": d["key"], "collection_era": era})
    fits: list[dict] = []
    rng = np.random.default_rng(0)

    for name, plat, fluid, metric, p1, p2 in BLOCKS:
        n1 = int(d[list(maps[p1].values())].notna().any(axis=1).sum())
        n2 = int(d[list(maps[p2].values())].notna().any(axis=1).sum())
        R, M = (p1, p2) if n1 >= n2 else (p2, p1)
        tR, tM = maps[R], maps[M]
        core = sorted(set(tR) & set(tM))
        cR = [tR[t] for t in core]
        cM = [tM[t] for t in core]

        hR = d[list(tR.values())].notna().any(axis=1)
        hM = d[list(tM.values())].notna().any(axis=1)
        nR = d[hR].groupby("cell").size()
        nM = d[hM].groupby("cell").size()
        cells = sorted(set(nR[nR >= MIN_CELL].index) & set(nM[nM >= MIN_CELL].index))

        log(f"\n{name}: reference {R} ({max(n1, n2):,} rows), mapped {M} "
            f"({min(n1, n2):,} rows), {len(core):,} core analytes, {len(cells)} matched cells")
        if len(cells) < 2:
            log("  SKIPPED — fewer than 2 matched cells")
            continue

        gR, gM = d[hR].groupby("cell"), d[hM].groupby("cell")
        medR = gR[cR].median().reindex(cells).to_numpy()
        medM = gM[cM].median().reindex(cells).to_numpy()
        iqrR = (gR[cR].quantile(.75) - gR[cR].quantile(.25)).reindex(cells).to_numpy()
        iqrM = (gM[cM].quantile(.75) - gM[cM].quantile(.25)).reindex(cells).to_numpy()
        wt = np.minimum(nR[cells].to_numpy(), nM[cells].to_numpy()).astype(float)

        # Scope of the per-analyte weight mask below: an analyte absent from a cell has
        # no median there, and only those analytes are affected by masking the weights.
        gap = np.isnan(medR) | np.isnan(medM)
        log(f"  {int(gap.sum()):,} of {gap.size:,} analyte-cells carry no median; "
            f"{int(gap.any(axis=0).sum()):,} of {len(core):,} analytes are missing in "
            f"at least one cell")

        def combine(idx):
            # Cell-size-weighted mean down the cell axis, per analyte.
            #
            # The denominator must be summed over the SAME cells the numerator used.
            # Until 2026-08-15 this divided every statistic by the total weight of all
            # cells while np.nansum skipped the cells where that analyte had no median
            # — so an analyte measured in only some cells had its median and IQR pulled
            # toward zero in proportion to its missingness. Because the slope is IR/IM,
            # the bias partly cancels there, but it does not cancel in the intercept
            # (MR - slope*MM), and it does not cancel when the two projects have
            # different missingness patterns. Raised by the stats core; the weights are
            # now masked per analyte before both the numerator and the denominator.
            W = wt[idx][:, None]

            def wmean(mat):
                m = mat[idx]
                Wm = np.where(np.isnan(m), 0.0, W)      # broadcasts to (cells, analytes)
                s = Wm.sum(axis=0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    return np.nansum(m * Wm, axis=0) / np.where(s > 0, s, np.nan)

            with np.errstate(invalid="ignore", divide="ignore"):
                MR, MM = wmean(medR), wmean(medM)
                IR, IM = wmean(iqrR), wmean(iqrM)
                sl = IR / IM
            sl = np.where(np.isfinite(sl) & (IM > 1e-9), sl, 1.0)
            return sl, MR, MM

        slope, MR, MM = combine(np.arange(len(cells)))

        # bootstrap over cells: only trust a slope whose CI clears the tolerance band
        boots = np.vstack([combine(rng.integers(0, len(cells), len(cells)))[0]
                           for _ in range(N_BOOT)])
        lo = np.nanpercentile(boots, 2.5, axis=0)
        hi = np.nanpercentile(boots, 97.5, axis=0)
        flagged = np.abs(slope - 1) >= SLOPE_TOL
        supported = flagged & ((lo > 1 + SLOPE_TOL) | (hi < 1 - SLOPE_TOL))

        applied_slope = np.where(supported, slope, 1.0)
        applied_inter = MR - applied_slope * MM
        log(f"  slope flagged {int(flagged.sum()):,}, bootstrap-supported "
            f"{int(supported.sum()):,} ({supported.mean() * 100:.0f}% of analytes get a slope)")
        log(f"  median applied slope {np.nanmedian(applied_slope):.3f}, "
            f"median intercept {np.nanmedian(applied_inter):+.3f}")

        # --- build the harmonized columns ---
        vR = d[cR].to_numpy(dtype=float)
        vM = d[cM].to_numpy(dtype=float)
        corrected = applied_slope[None, :] * vM + applied_inter[None, :]
        harm = np.where(np.isnan(vR), corrected, vR)

        names = [f"harmonized_{name}_{t}_{metric}" for t in core]
        out = pd.concat([out, pd.DataFrame(harm, columns=names, index=d.index)], axis=1)

        # provenance: which project supplied this row's values
        any_r = ~np.isnan(vR).all(axis=1)
        any_m = ~np.isnan(corrected).all(axis=1)
        src_col = np.where(any_r, R, np.where(any_m, f"{M}_corrected", None))
        out[f"harmonized_{name}_src"] = pd.Series(src_col, index=d.index, dtype="string")

        n_from_ref = int(any_r.sum())
        n_from_map = int((~any_r & any_m).sum())
        log(f"  rows with data {n_from_ref + n_from_map:,} "
            f"({n_from_ref:,} from reference, {n_from_map:,} from corrected mapped "
            f"= {n_from_map / max(n_from_ref + n_from_map, 1) * 100:.0f}% corrected)")

        for j, t in enumerate(core):
            fits.append(dict(block=name, platform=plat, fluid=fluid, metric=metric,
                             analyte=t, reference_project=R, mapped_project=M,
                             col_reference=cR[j], col_mapped=cM[j],
                             harmonized_column=names[j],
                             fitted_slope=round(float(slope[j]), 4),
                             slope_ci_lo=round(float(lo[j]), 4),
                             slope_ci_hi=round(float(hi[j]), 4),
                             slope_supported=bool(supported[j]),
                             # full precision, deliberately unrounded: these are the
                             # machine-readable record of what was applied, so the export
                             # must be exactly reproducible from them. The dictionary
                             # formats them for display separately.
                             applied_slope=float(applied_slope[j]),
                             applied_intercept=float(applied_inter[j]),
                             n_matched_cells=len(cells),
                             pct_from_corrected=round(
                                 n_from_map / max(n_from_ref + n_from_map, 1) * 100, 1)))

    assert len(out) == n_rows, "row count changed during harmonization"
    assert out["key"].tolist() == d["key"].tolist(), "key order changed"
    assert not out.columns.duplicated().any(), "duplicate harmonized column names"

    f_out = os.path.join(OUT, f"harmonized-{TIMESTAMP}.tab")
    f_fit = os.path.join(OUT, f"harmonization_fit-{TIMESTAMP}.tab")
    log(f"\nwriting {os.path.basename(f_out)} ...")
    out.to_csv(f_out, sep="\t", index=False)
    pd.DataFrame(fits).to_csv(f_fit, sep="\t", index=False)
    log(f"wrote {os.path.basename(f_out)}  ({out.shape[0]:,} x {out.shape[1]:,})")
    log(f"wrote {os.path.basename(f_fit)}  ({len(fits):,} analyte fits)")


if __name__ == "__main__":
    sys.exit(main())
