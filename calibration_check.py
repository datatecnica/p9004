#!/usr/bin/env python3
"""calibration_check.py — is a run's p-value scale trustworthy?

Every other script in this pipeline asks "what is significant?". This one asks the
prior question: "are these p-values calibrated at all?" It answers it two ways, both
cheap, and both of which would have caught the 2026-07-31 audit findings at the time
they were introduced rather than months later.

  1. NEGATIVE CONTROLS (free — reads the result CSVs)
     The panels ship spike-in / assay control analytes (`*_CTRL_*`, `*_mCherry_*`).
     They cannot carry biological signal, so their rejection rate at alpha is the
     empirical type-I error of the run. Nominal is 5%. The 2026-07-30 batch ran at
     20% in the trajectory LMMs -- a 4x error that no amount of multiple-testing
     correction can repair, because Bonferroni scales a threshold, not a broken
     reference distribution.

  2. PERMUTATION NULL (costs fits — samples analytes and refits)
     Shuffle the group label BETWEEN SUBJECTS, within blocks of (n_obs, mean-time),
     then refit. Blocking matters: a free shuffle also destroys the association
     between group and follow-up structure, which is itself a source of inflation.
     In the 2026-07-30 spec, a free shuffle looked calibrated (lambda 1.05) while the
     blocked shuffle did not (lambda 1.29, type-I 9.5%) -- the imbalance alone was
     manufacturing signal. Under the current spec the blocked null returns 5.3%.

Usage:
    # fast: negative controls only, over an existing results dir
    python3 calibration_check.py --results-dir results

    # add a permutation null for specific runs (slow: refits n_sample analytes x2)
    python3 calibration_check.py --results-dir results --config batch.yaml \
        --permute trajectory_HC_vs_PD --n-sample 150
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CONTROL_RE = re.compile(r"(?:_CTRL_NPX|_mCherry_NPQ)$")
_FILENAME_RE = re.compile(r"^(?P<name>.+)-(?P<stratum>[A-Za-z0-9]+)-(?P<ts>\d{8}_\d{6})\.csv$")

# A run whose controls reject at more than this is not reporting usable p-values.
CONTROL_TYPE_I_LIMIT = 10.0  # percent, i.e. 2x nominal


def _flag(g: pd.DataFrame, col: str) -> pd.Series:
    """A results-CSV significance column as a clean boolean mask.

    The column may be absent, or round-trip through the CSV as the strings
    "True"/"False", so a bare truth test would count "False" as significant.
    """
    if col not in g.columns:
        return pd.Series(False, index=g.index)
    s = g[col]
    if s.dtype == object:
        return s.astype(str).str.strip().str.lower().eq("true")
    return s.fillna(False).astype(bool)


def lambda_gc(p: pd.Series | np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if len(p) < 20:
        return float("nan")
    return float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1))


# ============================================================================
# 1. Negative controls
# ============================================================================

def control_report(results_dir: str) -> pd.DataFrame:
    rows: list[dict] = []
    for path in sorted(glob.glob(os.path.join(results_dir, "*.csv"))):
        fname = os.path.basename(path)
        if fname.startswith("META_"):
            continue
        m = _FILENAME_RE.match(fname)
        if not m:
            continue
        df = pd.read_csv(path)
        if "P" not in df.columns:
            continue
        if "is_control" in df.columns:
            ctrl = df["is_control"] == True
        else:  # older CSVs without the flag: fall back to name matching
            lab = (df["outcome"].astype(str)
                   .where(df["outcome"].astype(str).str.contains("NPX|NPQ"),
                          df["predictor"].astype(str)))
            ctrl = lab.str.contains(CONTROL_RE, na=False)
        c = df[ctrl]
        if "degenerate" in c.columns:
            c = c[c["degenerate"] != True]
        if c.empty:
            continue
        terms = c.groupby("term") if "term" in c.columns else [("effect", c)]
        for term, g in terms:
            # 4, not 5: half the 8 control analytes are mCherry, which are
            # numerically degenerate and excluded above, leaving 4 usable _CTRL_
            # analytes. At a higher minimum the LMM runs — the ones most in need of
            # a calibration check — produce no evaluable cell at all. n_controls is
            # reported alongside so the (very wide) precision is visible.
            if len(g) < 4:
                continue
            rows.append({
                "run": m.group("name"),
                "stratum": m.group("stratum"),
                "term": term,
                "n_controls": len(g),
                "pct_p05": 100 * float((g["P"] < 0.05).mean()),
                "lambda": lambda_gc(g["P"]),
                # Controls significant under EITHER Bonferroni lens. The union is taken
                # over the ROWS and then counted; until 2026-08-15 this bitwise-OR'd the
                # two COUNTS, which is arithmetic on integers rather than a set union
                # (4 controls run-wide | 2 family-wide reported 6, and 5 | 3 reported 7,
                # neither of which can exceed n_controls). Raised twice by the stats core.
                "n_bonferroni": int(
                    (_flag(g, "significant") | _flag(g, "significant_family")).sum()
                ),
            })
    return pd.DataFrame(rows)


# ============================================================================
# 2. Blocked permutation null
# ============================================================================

def permutation_null(config_path: str, run_name: str, stratum: str,
                     n_sample: int, seed: int, n_perm: int = 1) -> pd.DataFrame:
    """Refit `n_sample` analytes with the group label shuffled between subjects.

    The shuffle is blocked on (n_obs, mean-time tertile) so the permuted groups keep
    the real design's follow-up imbalance. Any inflation that survives is the model,
    not the biology.

    Repeated `n_perm` times, returning one row per fit with a `perm` column. A single
    permutation gives a point estimate with no sense of its own spread — the stats core
    asked how much of a verdict rests on which shuffle happened to be drawn, and with
    replicates that question is answerable from the output rather than assumed away.

    THE ANALYTE SAMPLE IS HELD FIXED across replicates (drawn once, from `seed`) while
    only the label shuffle varies (`seed + i`). Redrawing analytes each replicate would
    fold analyte-sampling variance into the spread and overstate permutation
    instability, which is the opposite of what the replicates are measuring.
    """
    import yaml
    import regressions as R

    warnings.filterwarnings("ignore")
    cfg = yaml.safe_load(open(config_path))
    defaults = cfg.get("defaults", {})
    spec = next((s for s in cfg["runs"] if s.get("name") == run_name), None)
    if spec is None:
        sys.exit(f"ERROR: run '{run_name}' not in {config_path}")

    cwd = os.path.dirname(os.path.abspath(config_path))
    df = R.load_data(R.resolve_input_path(defaults.get("input_glob", "unified_PPMI-*.tab"), cwd))
    proteomic = R.list_proteomic_columns(df)
    assay_of = R.proteomic_assay_map(proteomic)

    # The permuted label replaces the real grouping column everywhere it appears.
    report_term = spec.get("report_term", "")
    parts = [p.strip() for p in report_term.split(":")]
    time_col = spec.get("time_col", defaults.get("time_col", "YEAR"))
    group_col = next((p for p in parts if p != time_col), None)
    if group_col is None or group_col not in df.columns:
        sys.exit(f"ERROR: could not identify the group column from report_term '{report_term}'")

    strata_col = defaults.get("strata_col")
    work = df[df[strata_col] == stratum] if strata_col else df
    sample_filter = spec.get("sample_filter")
    if sample_filter:
        work = work.query(sample_filter)
    work = work.copy()

    # Fixed across replicates — see the docstring.
    cols = [c for c in proteomic if c in work.columns]
    sample = list(np.random.default_rng(seed).choice(
        cols, size=min(n_sample, len(cols)), replace=False))
    spec = dict(spec)
    spec["predictors"] = sample

    per = work.groupby("PATNO").agg(lab=(group_col, "first"),
                                    n_obs=(time_col, "size"),
                                    mt=(time_col, "mean"))
    blk = (per.n_obs.clip(upper=4).astype(str) + "|"
           + pd.qcut(per.mt, 3, labels=False, duplicates="drop").astype(str))
    blocks = list(per.groupby(blk).groups.values())

    frames: list[pd.DataFrame] = []
    for i in range(n_perm):
        rng = np.random.default_rng(seed + i)
        permuted = per.lab.copy()
        for idx in blocks:
            permuted.loc[idx] = rng.permutation(per.lab.loc[idx].to_numpy())
        rep = work.copy()
        rep[group_col] = rep.PATNO.map(permuted.to_dict())

        out = R.run_one(rep, spec, None, defaults, set(proteomic), assay_of)
        out = R.filter_numerical_failures(out, f"NULL:{run_name}|{stratum}|perm{i}")
        if not out.empty:
            out = out.copy()
            out["perm"] = i
            frames.append(out)
        print(f"    permutation {i + 1}/{n_perm}: {len(out):,} fits"
              f"{'' if not out.empty else '  (none survived filtering)'}", flush=True)

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description="Calibration check for a results directory")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--config", help="batch.yaml — required for --permute")
    ap.add_argument("--permute", nargs="*", default=[], help="run name(s) to permutation-test")
    ap.add_argument("--stratum", default="EUR")
    ap.add_argument("--n-sample", type=int, default=150)
    ap.add_argument("--n-perm", type=int, default=10,
                    help="permutation replicates per run (default 10). Cost is linear: "
                         "each replicate refits --n-sample analytes. 1 reproduces the "
                         "pre-2026-08-15 single-shuffle behaviour.")
    ap.add_argument("--seed", type=int, default=20260731)
    ap.add_argument("--output", default="calibration_report.md")
    args = ap.parse_args()

    lines: list[str] = ["# Calibration report\n"]

    ctrl = control_report(args.results_dir)
    lines.append("## Negative controls\n")
    if ctrl.empty:
        lines.append("_No control analytes found in the results directory._\n")
        print("No control analytes found.")
    else:
        ctrl = ctrl.sort_values("pct_p05", ascending=False)
        bad = ctrl[ctrl.pct_p05 > CONTROL_TYPE_I_LIMIT]
        lines.append(f"Control analytes cannot carry signal, so their rejection rate at "
                     f"α=0.05 is the run's empirical type-I error. Nominal is **5%**; "
                     f"anything above **{CONTROL_TYPE_I_LIMIT:.0f}%** is flagged.\n")
        lines.append(f"**{len(bad)} of {len(ctrl)}** run × stratum × term cells exceed the limit.\n")
        lines.append("| Run | Stratum | Term | n controls | % P<0.05 | λ | Bonf hits | status |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for _, r in ctrl.iterrows():
            status = "⚠️ miscalibrated" if r.pct_p05 > CONTROL_TYPE_I_LIMIT else "ok"
            lines.append(f"| `{r['run']}` | {r['stratum']} | {r['term']} | {r['n_controls']} | "
                         f"{r['pct_p05']:.1f}% | {r['lambda']:.3f} | {r['n_bonferroni']} | {status} |")
        print(f"Negative controls: {len(bad)}/{len(ctrl)} cells exceed "
              f"{CONTROL_TYPE_I_LIMIT:.0f}% type-I.")

    if args.permute:
        if not args.config:
            sys.exit("ERROR: --permute requires --config")
        lines.append("\n## Blocked permutation null\n")
        lines.append("Group labels shuffled between subjects within blocks of "
                     "(n_obs, mean-time tertile), preserving the real follow-up imbalance. "
                     "A calibrated run returns λ ≈ 1 and ~5% at P<0.05.\n")
        lines.append(f"Each run is permuted **{args.n_perm}×** over the same fixed sample of "
                     f"{args.n_sample} analytes, so the spread below is permutation "
                     f"variability alone, not analyte sampling. A single replicate reports "
                     f"a point with no error bar; the **range** column is what says whether "
                     f"that point could be trusted.\n")
        lines.append("| Run | Term | perms | n fits | λ mean ± sd | λ range | "
                     "% P<0.05 mean ± sd | % range | warn |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        warn_total = 0
        for run_name in args.permute:
            res = permutation_null(args.config, run_name, args.stratum,
                                   args.n_sample, args.seed, args.n_perm)
            if res.empty:
                lines.append(f"| `{run_name}` | — | 0 | 0 | — | — | — | — | — |")
                continue
            terms = res["term"].unique() if "term" in res.columns else ["effect"]
            for term in terms:
                g = res[res["term"] == term] if "term" in res.columns else res
                per_perm = g.groupby("perm")["P"].agg(
                    lam=lambda p: lambda_gc(p),
                    pct=lambda p: 100 * float((p < 0.05).mean()),
                    n="size")
                lam, pct = per_perm["lam"], per_perm["pct"]
                # A replicate over the limit is a warning; replicates straddling it mean
                # the verdict depends on which shuffle was drawn, which is exactly the
                # instability a single permutation cannot reveal.
                n_over = int((pct > CONTROL_TYPE_I_LIMIT).sum())
                straddles = 0 < n_over < len(pct)
                warn_total += n_over
                flag = ("⚠️ unstable" if straddles else
                        "⚠️ inflated" if n_over == len(pct) else "ok")
                lines.append(
                    f"| `{run_name}` | {term} | {len(per_perm)} | {int(per_perm['n'].sum()):,} | "
                    f"{lam.mean():.3f} ± {lam.std():.3f} | {lam.min():.3f}–{lam.max():.3f} | "
                    f"{pct.mean():.1f}% ± {pct.std():.1f} | {pct.min():.1f}–{pct.max():.1f}% | "
                    f"{flag} |")
                print(f"  {run_name} [{term}]: λ={lam.mean():.3f}±{lam.std():.3f}  "
                      f"{pct.mean():.1f}%±{pct.std():.1f} at P<0.05  "
                      f"({n_over}/{len(pct)} replicates over "
                      f"{CONTROL_TYPE_I_LIMIT:.0f}%){'  UNSTABLE' if straddles else ''}")
        lines.append(f"\n**{warn_total}** permutation replicate(s) exceeded the "
                     f"{CONTROL_TYPE_I_LIMIT:.0f}% type-I limit across all runs tested. "
                     f"A run flagged `unstable` had replicates on both sides of it.\n")

    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
