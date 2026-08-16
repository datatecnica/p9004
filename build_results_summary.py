#!/usr/bin/env python3
"""build_results_summary.py — produce a results summary Markdown.

Reads every per-stratum CSV in results/ and every META CSV in meta/, then
emits two Markdown tables:

  1. Coverage table  — per run: model type, EUR/AJ N ranges, whether a META
                        was produced, λ-GC for EUR / AJ / META FE / META RE.
                        λ is recomputed from each file's P column (validates
                        against rather than trusting the original run log).
  2. Hits table      — every Bonferroni-significant row across the entire
                        analysis, grouped by run, sorted by P within each
                        source (EUR / AJ / META FE / META RE).

The number of tests for Bonferroni purposes is read from the
`bonferroni_threshold` column already stored by regressions.py and
meta_analysis.py — no recomputation, since those values are what the runs
were actually evaluated against.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import chi2

_FILENAME_RE = re.compile(r"^(?P<name>.+)-(?P<stratum>[^-]+)-(?P<ts>\d{8}_\d{6})\.csv$")
_META_RE = re.compile(r"^META_(?P<name>.+)-(?P<ts>\d{8}_\d{6})\.csv$")

# The PRS as it appears in the `predictor` column of every results CSV. Must track
# regressions.py NEW_PREDICTORS: this is matched by equality, so a stale name yields
# zero PRS-as-predictor hits silently rather than raising. It read GP2_PRS_zscore
# until 2026-08-15, which the p9005 genetics block had already retired — so the
# "PRS as the looped predictor" count was reported as 0 regardless of the results.
PRS_PREDICTOR = "p9005_Genetic_PRS_PRS157"


def lambda_gc(pvals) -> float:
    p = pd.Series(pvals).dropna()
    p = p[(p > 0) & (p <= 1)]
    if len(p) == 0:
        return float("nan")
    chi2_obs = chi2.ppf(1 - p, df=1)
    return float(np.median(chi2_obs) / chi2.ppf(0.5, df=1))


def discover_per_stratum(results_dir: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    by_ts: dict[tuple[str, str], str] = {}
    for path in sorted(glob.glob(os.path.join(results_dir, "*.csv"))):
        fname = os.path.basename(path)
        if fname.startswith("META_"):
            continue
        m = _FILENAME_RE.match(fname)
        if not m:
            continue
        key = (m.group("name"), m.group("stratum"))
        ts = m.group("ts")
        if key not in by_ts or by_ts[key] < ts:
            by_ts[key] = ts
            out.setdefault(m.group("name"), {})[m.group("stratum")] = path
    return out


def discover_meta(meta_dir: str) -> dict[str, str]:
    out: dict[str, str] = {}
    by_ts: dict[str, str] = {}
    for path in sorted(glob.glob(os.path.join(meta_dir, "META_*.csv"))):
        m = _META_RE.match(os.path.basename(path))
        if not m:
            continue
        if m.group("name") not in by_ts or by_ts[m.group("name")] < m.group("ts"):
            by_ts[m.group("name")] = m.group("ts")
            out[m.group("name")] = path
    return out


def n_range(df: pd.DataFrame) -> str:
    if "N" not in df.columns or df["N"].dropna().empty:
        return "—"
    n = df["N"].dropna().astype(int)
    if n.min() == n.max():
        return f"{n.min():,}"
    return f"[{n.min():,}, {n.max():,}]"


def model_type(df: pd.DataFrame) -> str:
    if "model_type" not in df.columns or df["model_type"].dropna().empty:
        return "—"
    vc = df["model_type"].value_counts()
    primary = vc.index[0]
    if len(vc) > 1:
        return f"{primary} (+{len(vc)-1} variants)"
    return str(primary)


def _pct_over(cell: object, limit: float) -> bool:
    """True if any percentage in a '12% / 7%'-style cell exceeds `limit`."""
    if not isinstance(cell, str) or "%" not in cell:
        return False
    vals = re.findall(r"([\d.]+)%", cell)
    return any(float(v) > limit for v in vals)


def lambda_by_term(df: pd.DataFrame, pcol: str = "P") -> str:
    """λ per reported term. LMM runs emit a within and a between term whose
    calibration differs materially (1.14 vs 0.70 in the same panel on the
    2026-07-30 batch), so a single pooled λ hides both."""
    if df is None or pcol not in df.columns:
        return "—"
    if "term" not in df.columns or df["term"].nunique(dropna=True) <= 1:
        return f"{lambda_gc(df[pcol]):.3f}"
    parts = [f"{t}={lambda_gc(g[pcol]):.3f}"
             for t, g in df.groupby("term", sort=False)]
    return " / ".join(parts)


def control_type_i(df: pd.DataFrame, pcol: str = "P") -> str:
    """Empirical type-I rate on the spike-in / assay control analytes.

    These analytes cannot carry biological signal, so their rejection rate at
    α=0.05 is the calibration of the run. Degenerate fits are excluded because a
    numerically-zero coefficient is a separate failure mode, already flagged.
    """
    if df is None or "is_control" not in df.columns or pcol not in df.columns:
        return "—"
    c = df[df["is_control"] == True]
    if "degenerate" in c.columns:
        c = c[c["degenerate"] != True]
    if len(c) < 4:   # see calibration_check.py: only 4 non-degenerate controls exist
        return "—"
    if "term" in c.columns and c["term"].nunique(dropna=True) > 1:
        return " / ".join(
            f"{t}={100 * (g[pcol] < 0.05).mean():.0f}% (n={len(g)})" for t, g in c.groupby("term", sort=False)
        )
    return f"{100 * (c[pcol] < 0.05).mean():.0f}% (n={len(c)})"


def collect_sig(df: pd.DataFrame, run: str, source: str) -> pd.DataFrame:
    """Return per-row significant hits for a single source."""
    if df is None or df.empty:
        return pd.DataFrame()
    if source.startswith("META"):
        flag = "significant_FE" if source == "META_FE" else "significant_RE"
        beta_col = "beta_FE" if source == "META_FE" else "beta_RE"
        se_col = "SE_FE" if source == "META_FE" else "SE_RE"
        p_col = "P_FE" if source == "META_FE" else "P_RE"
    else:
        flag = "significant"
        beta_col = "beta"
        se_col = "SE"
        p_col = "P"
    if flag not in df.columns:
        return pd.DataFrame()

    # Two Bonferroni lenses, emitted by regressions.py / meta_analysis.py:
    #   run-wide -- 0.05 / every predictor fitted in the run
    #   family   -- 0.05 / predictors sharing the same assay panel
    # Adding Projects 312/314 roughly doubled the run-wide denominator, so a hit can be
    # family-significant while falling short run-wide. Keep the union and label which
    # lens each row passes, rather than silently dropping the family-only hits.
    fam_flag = f"{flag}_family"
    has_family = fam_flag in df.columns
    keep = df[flag] == True
    if has_family:
        keep = keep | (df[fam_flag] == True)

    sig = df[keep].copy()
    if sig.empty:
        return pd.DataFrame()
    def _col(name, default=pd.NA):
        return sig[name] if name in sig.columns else default

    out = pd.DataFrame({
        "run": run,
        "source": source,
        "predictor": sig.get("predictor"),
        "outcome": sig.get("outcome"),
        "assay_family": _col("assay_family"),
        "beta": sig.get(beta_col),
        "SE": sig.get(se_col),
        "P": sig.get(p_col),
        "threshold": sig.get("bonferroni_threshold"),
        "threshold_family": _col("bonferroni_threshold_family"),
        "sig_runwide": sig[flag].astype(bool),
        "sig_family": sig[fam_flag].astype(bool) if has_family else sig[flag].astype(bool),
        "model_type": sig.get("model_type", "META"),
        # Provenance/quality flags carried through so the hits table can separate
        # real findings from controls and from structurally weak fits.
        "term": _col("term", "effect"),
        "is_control": _col("is_control", False),
        "degenerate": _col("degenerate", False),
        "n_multi_obs": _col("n_multi_obs"),
        "robust_FE": _col("robust_FE"),
        "strata_concordant": _col("strata_concordant"),
        "stratum_dominated": _col("stratum_dominated"),
    })
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", default="results")
    p.add_argument("--meta-dir", default="meta")
    p.add_argument("--output", default="results_summary.md")
    args = p.parse_args()

    per_stratum = discover_per_stratum(args.results_dir)
    metas = discover_meta(args.meta_dir)
    run_names = sorted(set(per_stratum.keys()) | set(metas.keys()))

    coverage_rows: list[dict] = []
    hits_frames: list[pd.DataFrame] = []
    panel_lambda_rows: list[dict] = []

    for name in run_names:
        eur_path = per_stratum.get(name, {}).get("EUR")
        aj_path = per_stratum.get(name, {}).get("AJ")
        meta_path = metas.get(name)

        eur_df = pd.read_csv(eur_path) if eur_path else None
        aj_df = pd.read_csv(aj_path) if aj_path else None
        meta_df = pd.read_csv(meta_path) if meta_path else None

        # Pick a model_type to display from whichever source has it
        mt = "—"
        for d in (eur_df, aj_df):
            if d is not None and "model_type" in d.columns and d["model_type"].dropna().size:
                mt = model_type(d)
                break

        coverage_rows.append({
            "run": name,
            "model": mt,
            "EUR_N": n_range(eur_df) if eur_df is not None else "—",
            "AJ_N": n_range(aj_df) if aj_df is not None else "—",
            "meta": "✓" if meta_df is not None else "—",
            "λ_EUR": lambda_by_term(eur_df),
            "λ_AJ":  lambda_by_term(aj_df),
            "λ_FE":  lambda_by_term(meta_df, "P_FE"),
            "λ_RE":  lambda_by_term(meta_df, "P_RE"),
            # Empirical calibration: rejection rate on analytes that cannot carry signal.
            "ctrl_EUR": control_type_i(eur_df),
            "ctrl_AJ":  control_type_i(aj_df),
            "n_eur_rows": int(len(eur_df)) if eur_df is not None else 0,
            "n_aj_rows":  int(len(aj_df))  if aj_df  is not None else 0,
            "n_meta_rows": int(len(meta_df)) if meta_df is not None else 0,
        })

        # λ per (term, assay panel). Prefer the meta FE column; fall back to EUR
        # for runs that produced only one stratum.
        src_df, pcol, src_label = (
            (meta_df, "P_FE", "META_FE") if meta_df is not None and "P_FE" in meta_df.columns
            else (eur_df, "P", "EUR")
        )
        if src_df is not None and "assay_family" in src_df.columns and pcol in src_df.columns:
            tgroups = (src_df.groupby("term", sort=False)
                       if "term" in src_df.columns else [("effect", src_df)])
            for term, tdf in tgroups:
                for fam, fdf in tdf.groupby("assay_family", sort=True):
                    if len(fdf) < 20:
                        continue
                    panel_lambda_rows.append({
                        "run": name, "term": term, "source": src_label,
                        "panel": fam, "lam": lambda_gc(fdf[pcol]), "n": len(fdf),
                    })

        if eur_df is not None: hits_frames.append(collect_sig(eur_df, name, "EUR"))
        if aj_df  is not None: hits_frames.append(collect_sig(aj_df,  name, "AJ"))
        if meta_df is not None:
            hits_frames.append(collect_sig(meta_df, name, "META_FE"))
            hits_frames.append(collect_sig(meta_df, name, "META_RE"))

    cov = pd.DataFrame(coverage_rows)
    hits_raw = pd.concat([h for h in hits_frames if not h.empty], ignore_index=True) if hits_frames else pd.DataFrame()

    # Dedupe: collapse multiple rows for the same (run, predictor, outcome)
    # tuple into one entry that keeps the row with the smallest P, plus a
    # comma-separated `sources` field listing every place it was flagged
    # significant. Then sort by P globally so the strongest signals are at
    # the top of the table regardless of which run they came from.
    hits = pd.DataFrame()
    if not hits_raw.empty:
        keep_idx = hits_raw.groupby(
            ["run", "predictor", "outcome"], sort=False
        )["P"].idxmin()
        hits = hits_raw.loc[keep_idx].copy()
        sources_per = (
            hits_raw.groupby(["run", "predictor", "outcome"])["source"]
            .apply(lambda s: ", ".join(sorted(set(s))))
            .rename("sources")
            .reset_index()
        )
        hits = hits.merge(sources_per, on=["run", "predictor", "outcome"], how="left")
        # A tuple counts as significant under a lens if ANY of its collapsed source
        # rows passed that lens.
        lens_any = (hits_raw.groupby(["run", "predictor", "outcome"])[["sig_runwide", "sig_family"]]
                    .any().reset_index())
        hits = hits.drop(columns=["sig_runwide", "sig_family"]).merge(
            lens_any, on=["run", "predictor", "outcome"], how="left")
        hits = hits.sort_values("P").reset_index(drop=True)
        hits["β"] = hits["beta"].map(lambda x: f"{x:+.3g}")
        hits["SE"] = hits["SE"].map(lambda x: f"{x:.3g}")
        hits["P"] = hits["P"].map(lambda x: f"{x:.2e}")
        hits["α (Bonf)"] = hits["threshold"].map(lambda x: f"{x:.2e}" if pd.notna(x) else "—")
        hits["α (family)"] = hits["threshold_family"].map(
            lambda x: f"{x:.2e}" if pd.notna(x) else "—")
        hits["lens"] = [
            "both" if rw and fam else ("run-wide" if rw else "family")
            for rw, fam in zip(hits["sig_runwide"], hits["sig_family"])
        ]
    n_unique_hits = len(hits)
    n_raw_rows = len(hits_raw)
    n_runwide = int(hits["sig_runwide"].sum()) if not hits.empty else 0
    n_family = int(hits["sig_family"].sum()) if not hits.empty else 0
    n_both = int((hits["sig_runwide"] & hits["sig_family"]).sum()) if not hits.empty else 0

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    n_hits = n_unique_hits if not hits.empty else 0
    n_runs_with_hits = hits["run"].nunique() if not hits.empty else 0

    lines: list[str] = []
    lines.append("# Proteomics Data Mine — Results Summary\n")
    lines.append(f"Generated: {timestamp}\n")
    lines.append(f"Inputs: `{args.results_dir}/` ({sum(d['n_eur_rows'] + d['n_aj_rows'] for d in coverage_rows)} per-stratum rows across "
                 f"{sum(1 for d in coverage_rows if d['n_eur_rows']) + sum(1 for d in coverage_rows if d['n_aj_rows'])} CSVs), "
                 f"`{args.meta_dir}/` ({sum(d['n_meta_rows'] for d in coverage_rows)} META rows across "
                 f"{sum(1 for d in coverage_rows if d['n_meta_rows'])} CSVs).\n")
    lines.append(f"Total runs surveyed: **{len(coverage_rows)}**.  "
                 f"Bonferroni-significant unique hits: **{n_hits}** "
                 f"`(run, predictor, outcome)` tuples across "
                 f"**{n_runs_with_hits}** runs"
                 + (f" (deduplicated from {n_raw_rows} per-source rows)." if n_raw_rows > n_hits else ".")
                 + "\n")
    lines.append(
        f"Two Bonferroni lenses are reported. **Run-wide** divides 0.05 by every predictor "
        f"fitted in the run ({n_runwide} hits); **family** divides by the predictors in the "
        f"same assay panel ({n_family} hits), which keeps each panel comparable with earlier "
        f"cuts now that Projects 312/314 have roughly doubled the run-wide denominator. "
        f"{n_both} hits pass both. The `lens` column marks which.\n")

    lines.append("\n## Coverage table\n")
    lines.append("`λ` (genomic inflation) is recomputed here from the P column of each CSV "
                 "(validates the regression-script log values). LMM runs decompose time into "
                 "a **within**-subject (trajectory) and a **between**-subject (cross-sectional) "
                 "term, reported as `within=… / between=…`; the two are separate questions and "
                 "are calibrated separately.\n")
    lines.append("| Run | Model | EUR N | AJ N | Meta? | λ EUR | λ AJ | λ FE | λ RE |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for r in coverage_rows:
        lines.append(f"| `{r['run']}` | {r['model']} | {r['EUR_N']} | {r['AJ_N']} | {r['meta']} | "
                     f"{r['λ_EUR']} | {r['λ_AJ']} | {r['λ_FE']} | {r['λ_RE']} |")

    # ---------- λ by assay panel ----------
    # A run-level λ averages panels that differ enormously in repeat-measure depth
    # and platform behaviour. On the 2026-07-30 batch one trajectory run spanned
    # λ 0.86 (p293_olink_plasma) to 3.93 (p288_CNS_plasma) while reporting 1.38
    # overall, so the panel driving the inflation was invisible at run level.
    if panel_lambda_rows:
        pl = pd.DataFrame(panel_lambda_rows)
        panels = sorted(pl["panel"].unique())
        lines.append("\n## λ by assay panel\n")
        lines.append(
            "Genomic inflation per **run × term × assay panel**, computed from the "
            "meta fixed-effect P where a meta exists and from EUR otherwise. A "
            "run-level λ averages panels with very different repeat-measure depth "
            "and platform behaviour, so the panel actually driving an inflated run "
            "is only visible here. Cells are blank where a panel has <20 tests in "
            "that run. **Bold** marks λ > 1.5.\n")
        lines.append(
            "Read small panels with the right error bar: λ estimated from *n* tests has "
            "a sampling SE of roughly `2.33/√n` under the null, so the ~130–250-test "
            "panels (p282, p288, p312) carry SE ≈ 0.15–0.21 against ≈ 0.03 for the "
            "~5,400-test panels (p277, p293, p314). That is wide, but not wide enough "
            "to explain a λ of 2+ — those sit 5–8 SE above 1. Note also that λ > 1 is "
            "*not* by itself evidence of miscalibration: widespread true signal inflates "
            "λ too, and it grows with power. Use the negative controls and the "
            "permutation null in `calibration_check.py` to tell the two apart.\n")
        lines.append("| Run | Term | Source | " + " | ".join(p.replace("_", " ") for p in panels) + " |")
        lines.append("|---|---|---|" + "---|" * len(panels))
        for (run, term, src), g in pl.groupby(["run", "term", "source"], sort=False):
            by_panel = dict(zip(g["panel"], g["lam"]))
            cells = []
            for p in panels:
                v = by_panel.get(p)
                if v is None or not np.isfinite(v):
                    cells.append("")
                elif v > 1.5:
                    cells.append(f"**{v:.2f}**")
                else:
                    cells.append(f"{v:.2f}")
            lines.append(f"| `{run}` | {term} | {src} | " + " | ".join(cells) + " |")
        n_infl = int((pl["lam"] > 1.5).sum())
        lines.append(f"\n**{n_infl} of {len(pl)}** run × term × panel cells exceed λ = 1.5.\n")

    # ---------- Calibration against the built-in negative controls ----------
    lines.append("\n## Calibration — negative controls\n")
    lines.append(
        "The panels ship spike-in / assay control analytes (`*_CTRL_*`, `*_mCherry_*`) that "
        "cannot carry biological signal. They are **deliberately retained** in all outputs and "
        "flagged via `is_control`, because their rejection rate at α=0.05 is the empirical "
        "calibration of each run. A well-calibrated run sits near **5%**; anything above ~10% "
        "means the reported p-values are not trustworthy for that run, regardless of how small "
        "they are. Numerically-degenerate fits (`degenerate`, |β| < 1e-10) are excluded from "
        "this rate since they are a separate, separately-flagged failure mode.\n")
    flagged = [r for r in coverage_rows
               if any(_pct_over(r.get(k), 10.0) for k in ("ctrl_EUR", "ctrl_AJ"))]
    lines.append("| Run | control type-I EUR | control type-I AJ | status |")
    lines.append("|---|---|---|---|")
    for r in coverage_rows:
        bad = any(_pct_over(r.get(k), 10.0) for k in ("ctrl_EUR", "ctrl_AJ"))
        lines.append(f"| `{r['run']}` | {r['ctrl_EUR']} | {r['ctrl_AJ']} | "
                     f"{'⚠️ >2× nominal' if bad else 'ok'} |")
    lines.append(f"\n**{len(flagged)} of {len(coverage_rows)} runs** exceed 2× the nominal "
                 f"rate on at least one stratum.\n")

    # ---------- PRS involvement subsection ----------
    if not hits_raw.empty:
        # Use the deduped tuples for "unique hits"
        prs_interaction = hits[hits["run"].str.contains("_x_PRS", regex=False)] \
                              if not hits.empty else pd.DataFrame()
        prs_predictor   = hits[hits["predictor"] == PRS_PREDICTOR] \
                              if not hits.empty else pd.DataFrame()
        prs_union = pd.concat([prs_interaction, prs_predictor]).drop_duplicates(
            subset=["run", "predictor", "outcome"]) if not hits.empty else pd.DataFrame()
        lines.append("\n## PRS involvement\n")
        lines.append(f"- **PRS-interaction runs** (name contains `_x_PRS`): "
                     f"**{len(prs_interaction)}** unique hits "
                     f"(reported β = `predictor:{PRS_PREDICTOR}` interaction term).")
        lines.append(f"- **PRS as the looped predictor** (`{PRS_PREDICTOR}` main effect): "
                     f"**{len(prs_predictor)}** unique hits.")
        lines.append(f"- **Union** (any PRS involvement): **{len(prs_union)}** unique hits.\n")
        if not prs_interaction.empty:
            lines.append("**Interaction-run breakdown:**\n")
            lines.append("| Run | # sig hits |")
            lines.append("|---|---|")
            for run, n in prs_interaction["run"].value_counts().items():
                lines.append(f"| `{run}` | {n} |")
            lines.append("")
        if not prs_predictor.empty:
            lines.append("**PRS-as-predictor breakdown:**\n")
            lines.append("| Run | # sig hits |")
            lines.append("|---|---|")
            for run, n in prs_predictor["run"].value_counts().items():
                lines.append(f"| `{run}` | {n} |")
            lines.append("")

    lines.append("\n## Bonferroni-significant hits\n")
    if hits.empty:
        lines.append("_No significant hits across any run/source._\n")
    else:
        is_ctrl = hits.get("is_control", pd.Series(False, index=hits.index)).fillna(False).astype(bool)
        is_deg = hits.get("degenerate", pd.Series(False, index=hits.index)).fillna(False).astype(bool)
        real_hits, ctrl_hits = hits[~is_ctrl], hits[is_ctrl]

        lines.append("Each row is one unique `(run, predictor, outcome)` tuple, deduped "
                     "across EUR / AJ / META FE / META RE — the `Sources` column lists "
                     "every place it was flagged significant. β/SE/P shown are from the "
                     "**most-significant source** for that tuple. `α (Bonf)` is the "
                     "run-wide threshold and `α (family)` the within-assay-panel one; "
                     "`Lens` says which the hit clears (`both`, `run-wide`, or `family`). "
                     "Sorted by P globally so the strongest signals are at the top.\n")
        lines.append("For LMM runs the `Term` column says whether the hit is a **within**-subject "
                     "trajectory difference or a **between**-subject cross-sectional one. "
                     "`⚠` marks a numerically-degenerate fit (|β| < 1e-10), retained for "
                     "transparency but never a finding.\n")
        if len(ctrl_hits):
            lines.append(f"**{len(ctrl_hits)} control-analyte rows are excluded from the table "
                         f"below and listed separately** — see *Control analytes reaching "
                         f"significance*.\n")

        def _emit(frame):
            lines.append("| Run | Term | Predictor | Family | Outcome | β | SE | P | α (Bonf) | "
                         "α (family) | Lens | Sources | Model |")
            lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for _, r in frame.iterrows():
                fam = r.get("assay_family")
                fam = f"`{fam}`" if pd.notna(fam) else "—"
                term = r.get("term")
                term = str(term) if pd.notna(term) else "—"
                if bool(r.get("degenerate", False)):
                    term += " ⚠"
                lines.append(f"| `{r['run']}` | {term} | `{r['predictor']}` | {fam} | `{r['outcome']}` | "
                             f"{r['β']} | {r['SE']} | {r['P']} | {r['α (Bonf)']} | "
                             f"{r['α (family)']} | {r['lens']} | "
                             f"{r['sources']} | {r.get('model_type','—')} |")

        _emit(real_hits)

        if len(ctrl_hits):
            lines.append("\n## Control analytes reaching significance\n")
            lines.append(
                "These are spike-in / assay controls and **are not findings**. They are listed "
                "because a control clearing Bonferroni is direct evidence that the run's "
                "p-values are miscalibrated — the most useful diagnostic the pipeline "
                "produces. An entry here should be read as a warning about its run, not as "
                "a result.\n")
            _emit(ctrl_hits)

    out_path = args.output
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {out_path}: {len(coverage_rows)} runs in coverage, {n_hits} hits in summary.")


if __name__ == "__main__":
    main()
