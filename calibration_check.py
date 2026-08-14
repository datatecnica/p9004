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
                "n_bonferroni": int(
                    (g.get("significant", pd.Series(False, index=g.index)) == True).sum()
                    | (g.get("significant_family", pd.Series(False, index=g.index)) == True).sum()
                ),
            })
    return pd.DataFrame(rows)


# ============================================================================
# 2. Blocked permutation null
# ============================================================================

def permutation_null(config_path: str, run_name: str, stratum: str,
                     n_sample: int, seed: int) -> pd.DataFrame:
    """Refit `n_sample` analytes with the group label shuffled between subjects.

    The shuffle is blocked on (n_obs, mean-time tertile) so the permuted groups keep
    the real design's follow-up imbalance. Any inflation that survives is the model,
    not the biology.
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

    rng = np.random.default_rng(seed)
    per = work.groupby("PATNO").agg(lab=(group_col, "first"),
                                    n_obs=(time_col, "size"),
                                    mt=(time_col, "mean"))
    blk = (per.n_obs.clip(upper=4).astype(str) + "|"
           + pd.qcut(per.mt, 3, labels=False, duplicates="drop").astype(str))
    permuted = per.lab.copy()
    for _, idx in per.groupby(blk).groups.items():
        permuted.loc[idx] = rng.permutation(per.lab.loc[idx].to_numpy())
    work[group_col] = work.PATNO.map(permuted.to_dict())

    cols = [c for c in proteomic if c in work.columns]
    sample = list(rng.choice(cols, size=min(n_sample, len(cols)), replace=False))
    spec = dict(spec)
    spec["predictors"] = sample

    out = R.run_one(work, spec, None, defaults, set(proteomic), assay_of)
    return R.filter_numerical_failures(out, f"NULL:{run_name}|{stratum}")


# ============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description="Calibration check for a results directory")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--config", help="batch.yaml — required for --permute")
    ap.add_argument("--permute", nargs="*", default=[], help="run name(s) to permutation-test")
    ap.add_argument("--stratum", default="EUR")
    ap.add_argument("--n-sample", type=int, default=150)
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
        lines.append("| Run | Term | n fits | λ (null) | % P<0.05 |")
        lines.append("|---|---|---|---|---|")
        for run_name in args.permute:
            res = permutation_null(args.config, run_name, args.stratum,
                                   args.n_sample, args.seed)
            if res.empty:
                lines.append(f"| `{run_name}` | — | 0 | — | — |")
                continue
            groups = res.groupby("term") if "term" in res.columns else [("effect", res)]
            for term, g in groups:
                lines.append(f"| `{run_name}` | {term} | {len(g)} | {lambda_gc(g['P']):.3f} | "
                             f"{100 * (g['P'] < 0.05).mean():.1f}% |")
                print(f"  {run_name} [{term}]: λ={lambda_gc(g['P']):.3f}  "
                      f"{100 * (g['P'] < 0.05).mean():.1f}% at P<0.05")

    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
