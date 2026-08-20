#!/usr/bin/env python3
"""meta_analysis.py — IVW fixed-effects + DerSimonian-Laird random-effects
meta-analysis over per-stratum regression results from regressions.py.

Auto-discovers per-stratum CSVs in a results directory by grouping on the
`<run_name>-<stratum>-<timestamp>.csv` filename convention. For each run,
joins on (predictor, outcome), pools β/SE across studies, and emits a
timestamped META CSV with both FE and RE statistics plus heterogeneity
(Q, P_het, I²) and a per-study direction string.

Usage:
    python3 meta_analysis.py --results-dir results/ [--output-dir meta/]
    python3 meta_analysis.py --inputs X-EUR-*.csv X-AJ-*.csv --output META_X.csv
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import sys
from datetime import datetime
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm


log = logging.getLogger(__name__)
ALPHA = 0.05

# Filename parser: "<name>-<stratum>-<YYYYMMDD_HHMMSS>.csv"
# <name> may itself contain hyphens; stratum and timestamp are at the tail.
_FILENAME_RE = re.compile(r"^(?P<name>.+)-(?P<stratum>[^-]+)-(?P<ts>\d{8}_\d{6})\.csv$")


def setup_logging(log_file: str) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file)],
        force=True,
    )


# ============================================================================
# File discovery
# ============================================================================

# Must match regressions.NON_HT_PROTEOMICS_FAMILY: the Bonferroni family for every
# predictor without a high-throughput proteomics project prefix.
NON_HT_PROTEOMICS_FAMILY = "non_highthroughput_proteomics"


def discover_runs(results_dir: str) -> dict[str, list[tuple[str, str, str]]]:
    """Return {run_name: [(stratum, timestamp, path), ...]} from a directory.

    When multiple timestamps exist for the same (run_name, stratum) pair, only
    the latest one is kept — supports re-running a batch without stale files
    polluting the meta. Files that don't match the naming convention are
    skipped with a warning.
    """
    groups: dict[tuple[str, str], tuple[str, str]] = {}  # (name, stratum) -> (ts, path)
    for path in sorted(glob.glob(os.path.join(results_dir, "*.csv"))):
        fname = os.path.basename(path)
        if fname.startswith("META_"):
            continue  # skip prior meta outputs
        m = _FILENAME_RE.match(fname)
        if not m:
            log.warning(f"  skipping non-conforming filename: {fname}")
            continue
        key = (m.group("name"), m.group("stratum"))
        ts = m.group("ts")
        if key not in groups or groups[key][0] < ts:
            groups[key] = (ts, path)
    out: dict[str, list[tuple[str, str, str]]] = {}
    for (name, stratum), (ts, path) in groups.items():
        out.setdefault(name, []).append((stratum, ts, path))
    return out


# ============================================================================
# Meta math
# ============================================================================

def _finite(*arrs: np.ndarray) -> np.ndarray:
    m = np.ones_like(arrs[0], dtype=bool)
    for a in arrs:
        m &= np.isfinite(a)
    return m


def _direction_string(betas: np.ndarray) -> str:
    chars: list[str] = []
    for b in betas:
        if not np.isfinite(b):
            chars.append("?")
        elif b > 0:
            chars.append("+")
        elif b < 0:
            chars.append("-")
        else:
            chars.append("0")
    return "".join(chars)


def meta_row(betas: np.ndarray, ses: np.ndarray, ns: np.ndarray) -> dict:
    """Compute IVW-FE + DL-RE meta stats for one (predictor, outcome) tuple.

    `betas`, `ses`, `ns` are aligned arrays across studies; missing studies
    should be absent rather than NaN. Behavior:
      - 1 study → pass-through (FE = RE = that study; heterogeneity NaN)
      - 2+ studies with finite values → full FE + DL-RE calculation
      - 0 studies → returns all-NaN record
    """
    mask = _finite(betas, ses) & (ses > 0)
    betas = betas[mask]; ses = ses[mask]; ns = ns[mask]
    k = len(betas)
    out = {
        "n_studies": k,
        "N_total": int(np.nansum(ns)) if k else 0,
        "direction": _direction_string(betas),
        "beta_FE": np.nan, "SE_FE": np.nan, "Z_FE": np.nan, "P_FE": np.nan,
        "beta_RE": np.nan, "SE_RE": np.nan, "Z_RE": np.nan, "P_RE": np.nan,
        "Q": np.nan, "P_het": np.nan, "I2": np.nan, "tau2": np.nan,
    }
    if k == 0:
        return out

    w = 1.0 / (ses ** 2)
    beta_fe = float(np.sum(w * betas) / np.sum(w))
    se_fe = float(np.sqrt(1.0 / np.sum(w)))
    z_fe = beta_fe / se_fe
    out.update(beta_FE=beta_fe, SE_FE=se_fe, Z_FE=float(z_fe),
               P_FE=float(2 * norm.sf(abs(z_fe))))

    if k == 1:
        # Pass-through: RE = FE = the single study
        out.update(beta_RE=beta_fe, SE_RE=se_fe, Z_RE=float(z_fe), P_RE=out["P_FE"])
        return out

    # Heterogeneity (Cochran's Q)
    Q = float(np.sum(w * (betas - beta_fe) ** 2))
    dfree = k - 1
    P_het = float(chi2.sf(Q, dfree))
    I2 = float(max(0.0, (Q - dfree) / Q * 100)) if Q > 0 else 0.0

    # DerSimonian-Laird τ²
    sum_w = float(np.sum(w))
    sum_w2 = float(np.sum(w ** 2))
    denom = sum_w - (sum_w2 / sum_w) if sum_w > 0 else 0.0
    tau2 = max(0.0, (Q - dfree) / denom) if denom > 0 else 0.0

    # Random-effects weights
    w_star = 1.0 / (ses ** 2 + tau2)
    beta_re = float(np.sum(w_star * betas) / np.sum(w_star))
    se_re = float(np.sqrt(1.0 / np.sum(w_star)))
    z_re = beta_re / se_re

    out.update(Q=Q, P_het=P_het, I2=I2, tau2=tau2,
               beta_RE=beta_re, SE_RE=se_re, Z_RE=float(z_re),
               P_RE=float(2 * norm.sf(abs(z_re))))
    return out


# ============================================================================
# Combine one run's stratum files
# ============================================================================

def combine_run(run_name: str, files: list[tuple[str, str, str]]) -> pd.DataFrame:
    """Load each stratum CSV, outer-join on (predictor, outcome), meta-analyze."""
    strata_frames: dict[str, pd.DataFrame] = {}
    for stratum, _ts, path in files:
        df = pd.read_csv(path)
        required = {"predictor", "outcome", "beta", "SE", "P", "N"}
        missing = required - set(df.columns)
        if missing:
            log.warning(f"  {os.path.basename(path)} missing cols {missing}; skipping")
            continue
        strata_frames[stratum] = df

    if not strata_frames:
        return pd.DataFrame()

    # Build a unified index over all (predictor, outcome) pairs seen in any stratum
    index_df = pd.concat(
        [f[["predictor", "outcome"]] for f in strata_frames.values()], ignore_index=True
    ).drop_duplicates(subset=["predictor", "outcome"]).reset_index(drop=True)

    # Carry over is_proteomic / assay_family + per-stratum β/SE/P/N.
    # `term` and `analyte` also ride along: a decomposed LMM emits a within and a
    # between row per analyte, distinguished by the [term] suffix on `predictor`,
    # so the join key already separates them — but Bonferroni needs the analyte
    # count, and reporting needs the bare term label.
    is_proteomic_map: dict[tuple, bool] = {}
    assay_family_map: dict[tuple, str] = {}
    passthrough_maps: dict[str, dict[tuple, object]] = {
        c: {} for c in ("term", "analyte", "is_control", "degenerate", "n_multi_obs",
                        "model_type")
    }
    per_stratum: dict[str, pd.DataFrame] = {}
    for stratum, f in strata_frames.items():
        sub = f[["predictor", "outcome", "beta", "SE", "P", "N"]].rename(
            columns={"beta": f"beta_{stratum}", "SE": f"SE_{stratum}",
                     "P": f"P_{stratum}", "N": f"N_{stratum}"}
        )
        per_stratum[stratum] = sub
        if "is_proteomic" in f.columns:
            for _, row in f[["predictor", "outcome", "is_proteomic"]].drop_duplicates().iterrows():
                is_proteomic_map.setdefault((row["predictor"], row["outcome"]),
                                            bool(row["is_proteomic"]))
        if "assay_family" in f.columns:
            for _, row in f[["predictor", "outcome", "assay_family"]].drop_duplicates().iterrows():
                assay_family_map.setdefault((row["predictor"], row["outcome"]),
                                            str(row["assay_family"]))
        for col, store in passthrough_maps.items():
            if col not in f.columns:
                continue
            for _, row in f[["predictor", "outcome", col]].drop_duplicates().iterrows():
                store.setdefault((row["predictor"], row["outcome"]), row[col])

    merged = index_df
    for stratum, sub in per_stratum.items():
        merged = merged.merge(sub, on=["predictor", "outcome"], how="left")

    strata_order = sorted(strata_frames.keys())
    log.info(f"  {run_name}: {len(merged):,} (predictor, outcome) pairs across "
             f"strata {strata_order}")

    # Compute meta per row
    results: list[dict] = []
    for _, row in merged.iterrows():
        betas = np.array([row.get(f"beta_{s}", np.nan) for s in strata_order], dtype=float)
        ses = np.array([row.get(f"SE_{s}", np.nan) for s in strata_order], dtype=float)
        ns = np.array([row.get(f"N_{s}", np.nan) for s in strata_order], dtype=float)
        stats = meta_row(betas, ses, ns)
        stats["predictor"] = row["predictor"]
        stats["outcome"] = row["outcome"]
        key = (row["predictor"], row["outcome"])
        stats["is_proteomic"] = is_proteomic_map.get(key, False)
        stats["assay_family"] = assay_family_map.get(key, NON_HT_PROTEOMICS_FAMILY)
        for col, store in passthrough_maps.items():
            if store:
                stats[col] = store.get(key, np.nan)
        for s in strata_order:
            stats[f"beta_{s}"] = row.get(f"beta_{s}", np.nan)
            stats[f"SE_{s}"] = row.get(f"SE_{s}", np.nan)
            stats[f"P_{s}"] = row.get(f"P_{s}", np.nan)
            stats[f"N_{s}"] = row.get(f"N_{s}", np.nan)
        # Inverse-variance weight share of the single heaviest stratum, and whether
        # the weight ordering INVERTS the sample-size ordering.
        #
        # A large weight share is not by itself anomalous — with two strata one of
        # them always exceeds 50%, and the larger study is expected to dominate. The
        # meaningful pathology is the SMALLER stratum outweighing the larger, which
        # requires its SE to beat a bigger sample's. That happens when a fit's
        # variance collapses: FE weight is 1/SE², so one degenerate SE hijacks the
        # pooled estimate. On the 2026-07-30 batch AJ carried >50% of the weight in
        # 31% of trajectory hits despite ~1/8 of EUR's subjects.
        se_arr = np.array([row.get(f"SE_{s}", np.nan) for s in strata_order], dtype=float)
        n_arr = np.array([row.get(f"N_{s}", np.nan) for s in strata_order], dtype=float)
        ok = np.isfinite(se_arr) & (se_arr > 0) & np.isfinite(n_arr)
        if ok.sum() > 1:
            w = 1.0 / se_arr[ok] ** 2
            stats["max_stratum_weight_share"] = float(w.max() / w.sum())
            # True when the heaviest-weighted stratum is not the largest one.
            stats["weight_size_inverted"] = bool(np.argmax(w) != np.argmax(n_arr[ok]))
        else:
            stats["max_stratum_weight_share"] = 1.0 if ok.sum() else np.nan
            stats["weight_size_inverted"] = False
        results.append(stats)

    meta_df = pd.DataFrame(results)

    # Bonferroni on the meta set. Denominator is the number of ANALYTES tested
    # within each reported term, not the number of rows: a decomposed LMM emits a
    # within and a between row per analyte, and correcting on rows would halve the
    # threshold purely as an artefact of the output schema. The two terms are
    # separate questions, each corrected over its own analyte count.
    n = len(meta_df)
    if "term" in meta_df.columns and "analyte" in meta_df.columns and n:
        n_per_term = meta_df.groupby("term")["analyte"].transform("nunique")
        fam_n = meta_df.groupby(["term", "assay_family"])["analyte"].transform("nunique")
    else:
        n_per_term = pd.Series(n, index=meta_df.index)
        fam_n = (meta_df.groupby("assay_family")["P_FE"].transform("size")
                 if "assay_family" in meta_df.columns and n else n_per_term)
    meta_df["n_analytes_tested"] = n_per_term
    meta_df["bonferroni_threshold"] = ALPHA / n_per_term
    meta_df["significant_FE"] = meta_df["P_FE"] < meta_df["bonferroni_threshold"]
    meta_df["significant_RE"] = meta_df["P_RE"] < meta_df["bonferroni_threshold"]
    meta_df["bonferroni_threshold_family"] = ALPHA / fam_n
    meta_df["significant_FE_family"] = meta_df["P_FE"] < meta_df["bonferroni_threshold_family"]
    meta_df["significant_RE_family"] = meta_df["P_RE"] < meta_df["bonferroni_threshold_family"]

    # Robustness screen. FE is the primary estimate: with k=2 strata a DerSimonian-
    # Laird random-effects model is not a real random-effects model — τ² was exactly
    # 0 in 63% of rows on the 2026-07-30 batch (RE ≡ FE, no shrinkage) and estimated
    # from a single degree of freedom in the rest. RE is kept as a sensitivity
    # column, and heterogeneity is screened directly instead.
    _b = [f"beta_{s}" for s in strata_order if f"beta_{s}" in meta_df.columns]
    if len(_b) > 1:
        signs = np.sign(meta_df[_b].to_numpy())
        with np.errstate(invalid="ignore"):
            obs = np.where(np.isfinite(meta_df[_b].to_numpy()), signs, np.nan)
        concordant = pd.Series(
            [len({v for v in r if np.isfinite(v)}) <= 1 for r in obs], index=meta_df.index
        )
    else:
        concordant = pd.Series(True, index=meta_df.index)
    meta_df["strata_concordant"] = concordant
    # NOT "one stratum holds >50%" — with k=2 that is true of almost every row and
    # simply reflects the larger study carrying more weight, which is correct.
    # The flag fires only when the weight ordering inverts the sample-size ordering.
    meta_df["stratum_dominated"] = (
        (meta_df["n_studies"] > 1) & meta_df["weight_size_inverted"]
    )
    # Collapsed-variance guard, LMM runs only. The LMM outcome is z-scored, so
    # SE * sqrt(N) is comparable across analytes: median 1.45 on the 2026-08-01
    # batch against a 0.1st percentile of 0.28. A fit far below that has had its
    # residual variance collapse, which yields an implausibly small SE and an
    # extreme p-value — the four most significant trajectory fits sat at 0.028-0.085.
    # Deliberately NOT applied to non-LMM runs: SE is scale-dependent and the
    # slope_* models regress a raw clinical slope where small SEs are correct.
    LOW_VAR_SE_SCALED = 0.2
    is_lmm = (meta_df.get("model_type", pd.Series("", index=meta_df.index))
              .astype(str).str.startswith("LMM"))
    low_var = pd.Series(False, index=meta_df.index)
    for s in strata_order:
        bcol, ncol = f"SE_{s}", f"N_{s}"
        if bcol in meta_df.columns and ncol in meta_df.columns:
            scaled = meta_df[bcol] * np.sqrt(meta_df[ncol])
            low_var |= is_lmm & scaled.notna() & (scaled < LOW_VAR_SE_SCALED)
    meta_df["low_variance"] = low_var

    # Thin within-subject support. Only subjects with >=2 observations for an
    # analyte carry any within-subject information, so `n_multi_obs` — not N and not
    # n_groups — is the effective sample size for the WITHIN term. It does not apply
    # to the between term, where every subject contributes regardless of visit count.
    # Without this, the 2026-08-01 batch reported p312_Neuro_Plasma_NEFL_NPQ at
    # P = 2e-33 off 8 subjects' worth of repeat measurement.
    MIN_MULTI_OBS_WITHIN = 20
    nmo = pd.to_numeric(meta_df.get("n_multi_obs", pd.Series(np.nan, index=meta_df.index)),
                        errors="coerce")
    is_within = meta_df.get("term", pd.Series("", index=meta_df.index)).astype(str) == "within"
    meta_df["thin_within"] = is_within & nmo.notna() & (nmo < MIN_MULTI_OBS_WITHIN)

    bad = meta_df.get("degenerate", pd.Series(False, index=meta_df.index)).fillna(False).astype(bool)
    ctrl = meta_df.get("is_control", pd.Series(False, index=meta_df.index)).fillna(False).astype(bool)
    # Replication is a precondition, not one guard among several. Every other term
    # here is vacuously satisfied by a single-stratum row: `strata_concordant` is
    # trivially true with one direction to compare, `stratum_dominated` is explicitly
    # gated on n_studies > 1, `low_variance` is LMM-only and `thin_within` fires only
    # on within-terms. So without this an unreplicated result inherits the "robust"
    # label -- on the 2026-08-20 batch, 14 rows across the seven EUR-only run groups
    # (AJ too small to fit) came out robust_FE=True off a single stratum, including
    # four significant hits. n_studies is the count of strata that actually
    # contributed an estimate for THIS analyte, so it also catches per-analyte
    # single-stratum rows inside an otherwise two-stratum run group.
    meta_df["robust_FE"] = (
        (meta_df["n_studies"] >= 2)
        & meta_df["significant_FE_family"]
        & meta_df["strata_concordant"]
        & ~meta_df["stratum_dominated"]
        & ~bad & ~ctrl & ~low_var & ~meta_df["thin_within"]
    )

    # Column order: key → meta → heterogeneity → per-stratum → Bonferroni
    col_order = (
        ["predictor", "outcome", "analyte", "term", "n_studies", "N_total", "direction",
         "beta_FE", "SE_FE", "Z_FE", "P_FE",
         "beta_RE", "SE_RE", "Z_RE", "P_RE",
         "Q", "P_het", "I2", "tau2", "max_stratum_weight_share", "weight_size_inverted"]
        + [f"{c}_{s}" for s in strata_order for c in ("beta", "SE", "P", "N")]
        + ["is_proteomic", "assay_family", "is_control", "degenerate", "low_variance", "thin_within", "n_multi_obs", "model_type",
           "n_analytes_tested",
           "bonferroni_threshold", "significant_FE", "significant_RE",
           "bonferroni_threshold_family", "significant_FE_family", "significant_RE_family",
           "strata_concordant", "stratum_dominated", "robust_FE"]
    )
    col_order = [c for c in col_order if c in meta_df.columns]
    return meta_df[col_order]


# ============================================================================
# λ-GC (genomic inflation)
# ============================================================================

def lambda_gc(pvals: Iterable[float]) -> tuple[float, int]:
    p = pd.Series(pvals).dropna()
    p = p[(p > 0) & (p <= 1)]
    if len(p) == 0:
        return float("nan"), 0
    chi2_obs = chi2.ppf(1 - p, df=1)
    lam = float(np.median(chi2_obs) / chi2.ppf(0.5, df=1))
    return lam, len(p)


# ============================================================================
# Main driver
# ============================================================================

def run_auto(results_dir: str, output_dir: str, run_filter: str | None) -> None:
    groups = discover_runs(results_dir)
    if not groups:
        sys.exit(f"ERROR: no conforming CSVs found in {results_dir}")
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log.info(f"Discovered {len(groups)} run group(s) in {results_dir}")

    for name, files in sorted(groups.items()):
        if run_filter and name != run_filter:
            continue
        log.info("")
        log.info("=" * 70)
        log.info(f"META: {name}  ({len(files)} stratum files)")
        log.info("=" * 70)
        meta_df = combine_run(name, files)
        if meta_df.empty:
            log.warning(f"  {name}: no rows after combining; skipping")
            continue
        out = os.path.join(output_dir, f"META_{name}-{timestamp}.csv")
        meta_df.to_csv(out, index=False)
        n_sig_fe = int(meta_df["significant_FE"].sum())
        n_sig_re = int(meta_df["significant_RE"].sum())
        thr = float(meta_df["bonferroni_threshold"].iloc[0])
        log.info(f"  wrote {len(meta_df):,} rows "
                 f"({n_sig_fe} sig FE, {n_sig_re} sig RE @ {thr:.2e}) → {os.path.basename(out)}")
        # λ at three granularities. A pooled λ over a run averages panels with very
        # different depth and platform behaviour together — on the 2026-07-30 batch a
        # single trajectory run spanned λ 0.86 (p293) to 3.93 (p288_CNS_plasma) — and
        # for decomposed LMMs the within and between terms are separately calibrated.
        term_groups = (meta_df.groupby("term", sort=False)
                       if "term" in meta_df.columns else [("effect", meta_df)])
        for term, tdf in term_groups:
            lam_fe, n_fe = lambda_gc(tdf["P_FE"])
            lam_re, n_re = lambda_gc(tdf["P_RE"])
            log.info(f"  [{term}] λ_FE = {lam_fe:.4f}  (N={n_fe})     "
                     f"λ_RE = {lam_re:.4f}  (N={n_re})")
            if "assay_family" not in tdf.columns:
                continue
            for fam, fdf in tdf.groupby("assay_family", sort=True):
                lam_f, n_f = lambda_gc(fdf["P_FE"])
                if n_f < 20:
                    continue
                flag = "  <-- inflated" if lam_f > 1.5 else ""
                log.info(f"      {fam}: λ_FE = {lam_f:.4f}  N = {n_f}{flag}")


def run_explicit(inputs: list[str], output_path: str) -> None:
    # Explicit mode: user supplies the per-stratum files and output filename.
    # Stratum label is parsed from filename or set to the basename if non-conforming.
    files: list[tuple[str, str, str]] = []
    name = None
    for path in inputs:
        fname = os.path.basename(path)
        m = _FILENAME_RE.match(fname)
        if m:
            files.append((m.group("stratum"), m.group("ts"), path))
            if name is None:
                name = m.group("name")
        else:
            files.append((os.path.splitext(fname)[0], "00000000_000000", path))
    meta_df = combine_run(name or "manual", files)
    if meta_df.empty:
        sys.exit("ERROR: no rows after combining inputs")
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    meta_df.to_csv(output_path, index=False)
    log.info(f"Wrote {len(meta_df):,} rows → {output_path}")


# ============================================================================
# CLI
# ============================================================================

def main() -> None:
    p = argparse.ArgumentParser(description="Meta-analyze per-stratum regression results")
    p.add_argument("--results-dir", default="results",
                   help="Directory of per-stratum <name>-<stratum>-<ts>.csv files (auto-discover mode)")
    p.add_argument("--output-dir", default="meta",
                   help="Where to write META_<name>-<ts>.csv files (auto mode)")
    p.add_argument("--run", default=None,
                   help="Auto mode: only meta-analyze this named run group")
    p.add_argument("--inputs", nargs="+", default=None,
                   help="Explicit mode: per-stratum input files to combine")
    p.add_argument("--output", default=None,
                   help="Explicit mode: output CSV path (required when --inputs is used)")
    args = p.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = os.path.dirname(os.path.abspath(args.output or args.output_dir))
    os.makedirs(out_dir, exist_ok=True)
    log_file = os.path.join(out_dir, f"meta_analysis_{timestamp}.log")
    setup_logging(log_file)

    if args.inputs:
        if not args.output:
            sys.exit("ERROR: --output is required when using --inputs")
        run_explicit(args.inputs, args.output)
    else:
        run_auto(args.results_dir, args.output_dir, args.run)
    log.info("\nDone.")


if __name__ == "__main__":
    main()
