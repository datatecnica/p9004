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

    # add a permutation null for specific runs (slow: refits n_sample x n_perm analytes)
    python3 calibration_check.py --results-dir results --config batch.yaml \
        --permute trajectory_HC_vs_PD --n-sample 100 --n-perm 10 --workers 14

The permutation defaults are 100 analytes x 10 shuffles = 1,000 fits per run, drawn
from whatever `defaults.proteomic_scope` the config sets (harmonized, i.e. 11,567
analytes, for batch.yaml) so the null is built from the same universe the batch fits.
--workers applies to the trajectory LMMs only; every other model already runs serially
in about the time the fork would cost.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import warnings

# Pin BLAS to one thread per process BEFORE numpy imports, for the same reason and
# with the same values as regressions.py — see the note there. It has to be repeated
# here rather than inherited: `import regressions` happens inside permutation_null,
# long after this module has already imported numpy, so by the time regressions sets
# these OpenBLAS has read its thread count. Without this a --workers pool forks N
# processes each spawning 16 BLAS threads.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

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

_LHS_RE = re.compile(r"^\s*(.+?)\s*~")


def load_permutation_inputs(config_path: str) -> dict:
    """The release frame plus scope metadata, read ONCE for a whole sweep.

    Hoisted out of `permutation_null` on 2026-08-20, when the sweep went from the 10
    trajectory runs to all 52. Loading per run cost ~25 s and a 3.89 GB read peak each,
    so a full sweep would have spent 22 minutes in pure I/O and re-peaked memory 52
    times. The frame is read-only here — every run takes its own filtered copy.

    Read at the SAME analyte scope the batch ran at, resolved before the read the way
    regressions.main does. Until 2026-08-20 this loaded unscoped, which was wrong twice
    over. Statistically: it sampled from all 34,902 analytes while batch.yaml fits
    11,567, so the null characterised a universe the batch never touched. And in
    memory, which is what actually bit -- unscoped is 35,566 columns for a 5.15 GB
    frame at a ~11 GB read peak, against 12,231 columns / 1.77 GB / 3.89 GB under
    scope=harmonized. The trajectory LMMs are the tallest subsets in the batch and the
    replicate loop holds two more copies of one at a time, so the full-width read is
    what put the 31 GB box over.
    """
    import yaml
    import regressions as R

    cfg = yaml.safe_load(open(config_path))
    defaults = cfg.get("defaults", {})
    cwd = os.path.dirname(os.path.abspath(config_path))
    scope = defaults.get("proteomic_scope", "all")
    assay_prefixes = R.scoped_assay_prefixes(scope)
    df = R.load_data(
        R.resolve_input_path(defaults.get("input_glob", "unified_PPMI-*.tab"), cwd),
        assay_prefixes)
    proteomic = R.list_proteomic_columns(df, assay_prefixes)
    return {"cfg": cfg, "defaults": defaults, "scope": scope, "df": df,
            "proteomic": proteomic, "assay_of": R.proteomic_assay_map(proteomic)}


def permutation_target(spec: dict, defaults: dict, columns) -> tuple[list[str], str]:
    """Which column(s) a run's null shuffles, and how they travel.

    Returns `(cols, kind)`, or `([], reason)` for a run that cannot be permuted.

    Every run shuffles BETWEEN SUBJECTS; what differs is which column carries the
    assignment that has to be broken. Before 2026-08-20 only the first case existed
    and everything else exited, which is why the sweep covered 14 of 52 runs:

      interaction  `report_term` names a second term beside the time column -- the
                   trajectory LMMs (`YEAR:grp_*`) and the *_x_PRS OLS runs
                   (`PROTEIN:PRS157`). That term IS the reported contrast, so it is
                   what the null must destroy.
      outcome      `report_term` is a bare `PROTEIN`, meaning the protein itself is
                   the tested term. The association to break is then against the
                   formula's left-hand side -- the logit `grp_*` and OLS `slope_*`
                   runs.
      survival     Cox specs carry no formula at all; the outcome is the
                   (duration_col, event_col) pair. THE PAIR MUST MOVE TOGETHER -- a
                   shuffle that separates them manufactures follow-up times that were
                   never observed, and the resulting null is not the study's design.

    `PROTEIN` is a substitution token, not a column, so the `p in columns` test is what
    makes the interaction case fall through to `outcome` rather than matching on it.
    """
    model = str(spec.get("model", "")).lower()

    if model == "cox":
        dur = spec.get("duration_col", defaults.get("duration_col"))
        ev = spec.get("event_col", defaults.get("event_col"))
        pair = [c for c in (dur, ev) if c]
        missing = [c for c in pair if c not in columns]
        if len(pair) != 2 or missing:
            return [], f"cox spec needs duration_col+event_col present (missing {missing or pair})"
        return pair, "survival"

    time_col = spec.get("time_col", defaults.get("time_col", "YEAR"))
    parts = [p.strip() for p in str(spec.get("report_term", "")).split(":")]
    inter = [p for p in parts if p != time_col and p in columns]
    if inter:
        return inter[:1], "interaction"

    m = _LHS_RE.match(str(spec.get("formula", "")))
    if m and m.group(1) in columns:
        return [m.group(1)], "outcome"

    return [], (f"no permutable column from report_term "
                f"'{spec.get('report_term', '')}' or formula '{spec.get('formula', '')}'")


def permutation_null(ctx: dict, run_name: str, stratum: str,
                     n_sample: int, seed: int, n_perm: int = 1,
                     n_workers: int = 1) -> tuple[pd.DataFrame, str]:
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

    The cost parameters, stated once here because the reported false-positive range is
    only interpretable against them:

      n_sample   analytes drawn, from the `proteomic_scope` universe (default 100).
                 The rate is a binomial proportion over n_sample x n_perm fits, so
                 this and n_perm together set its precision -- 100 x 10 = 1,000 fits
                 gives roughly +/-1.4 points at a true 5%.
      n_perm     label shuffles (default 10). Sets the spread REPORTED in the range
                 column; n_perm=1 reports a point with no error bar.
      n_workers  processes over the analyte loop. LMM runs only, matching
                 R.PARALLEL_MODELS -- everything else is fast enough serially and
                 would only pay fork overhead. Results are unaffected: `imap`
                 preserves order and the fits are independent.
    """
    import regressions as R

    warnings.filterwarnings("ignore")
    defaults, df = ctx["defaults"], ctx["df"]
    spec = next((s for s in ctx["cfg"]["runs"] if s.get("name") == run_name), None)
    if spec is None:
        return pd.DataFrame(), f"no run named '{run_name}' in the config"

    target, kind = permutation_target(spec, defaults, df.columns)
    if not target:
        return pd.DataFrame(), kind  # `kind` carries the reason when target is empty

    strata_col = defaults.get("strata_col")
    work = df[df[strata_col] == stratum] if strata_col else df
    sample_filter = spec.get("sample_filter")
    if sample_filter:
        work = work.query(sample_filter)
    work = work.copy()

    min_n = spec.get("min_n", defaults.get("min_n", 20))
    if len(work) < min_n:
        return pd.DataFrame(), f"only {len(work)} rows after sample_filter (min_n={min_n})"

    # Fixed across replicates — see the docstring.
    cols = [c for c in ctx["proteomic"] if c in work.columns]
    if not cols:
        return pd.DataFrame(), "no in-scope analyte columns"
    sample = list(np.random.default_rng(seed).choice(
        cols, size=min(n_sample, len(cols)), replace=False))
    spec = dict(spec)
    spec["predictors"] = sample

    # One row per subject: the values to be shuffled, plus the block key. Blocking is on
    # (n_obs, mean-time tertile) for EVERY model, not just the longitudinal ones. The
    # single-visit runs still carry all their visits in `work` at this point -- the
    # one-row-per-subject reduction happens inside run_one, per predictor -- so n_obs
    # and mean-time remain the real follow-up structure there too, which is exactly the
    # imbalance a free shuffle would destroy.
    time_col = spec.get("time_col", defaults.get("time_col", "YEAR"))
    per = work.groupby("PATNO").agg(n_obs=(time_col, "size"), mt=(time_col, "mean"))
    lab = work.groupby("PATNO")[target].first()
    blk = (per.n_obs.clip(upper=4).astype(str) + "|"
           + pd.qcut(per.mt, 3, labels=False, duplicates="drop").astype(str))
    blocks = [idx for idx in per.groupby(blk).groups.values() if len(idx) > 1]

    # Same gate run_all applies: fork only the models whose analyte loop is worth it.
    run_workers = n_workers if str(spec["model"]).lower() in R.PARALLEL_MODELS else 1
    print(f"    {run_name}|{stratum} [{spec['model']}/{kind}]: {len(work):,} rows x "
          f"{len(sample)} analytes x {n_perm} perms on {run_workers} worker(s); "
          f"shuffling {'+'.join(target)} across {len(blocks)} blocks", flush=True)

    frames: list[pd.DataFrame] = []
    for i in range(n_perm):
        rng = np.random.default_rng(seed + i)
        permuted = lab.copy()
        for idx in blocks:
            # ONE subject order per block, applied to every target column. That is what
            # keeps a Cox (duration, event) pair together: shuffling the two columns
            # independently would pair a subject's follow-up time with another's event
            # status and invent observations the study never made.
            order = rng.permutation(len(idx))
            permuted.loc[idx, target] = lab.loc[idx, target].to_numpy()[order]
        rep = work.copy()
        for c in target:
            rep[c] = rep.PATNO.map(permuted[c])

        out = R.run_one(rep, spec, None, defaults, set(ctx["proteomic"]), ctx["assay_of"],
                        n_workers=run_workers)
        out = R.filter_numerical_failures(out, f"NULL:{run_name}|{stratum}|perm{i}")
        if not out.empty:
            out = out.copy()
            out["perm"] = i
            frames.append(out)
        print(f"    permutation {i + 1}/{n_perm}: {len(out):,} fits"
              f"{'' if not out.empty else '  (none survived filtering)'}", flush=True)

    if not frames:
        return pd.DataFrame(), "every fit was filtered out as a numerical failure"
    return pd.concat(frames, ignore_index=True), kind


# ============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description="Calibration check for a results directory")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--config", help="batch.yaml — required for --permute")
    ap.add_argument("--permute", nargs="*", default=[],
                    help="run name(s) to permutation-test, or 'all' for every run in "
                         "the config. A run with no permutable column is reported as "
                         "skipped rather than aborting the sweep.")
    ap.add_argument("--stratum", default="EUR")
    ap.add_argument("--n-sample", type=int, default=100,
                    help="analytes sampled per run (default 100), drawn from the "
                         "config's proteomic_scope universe. With --n-perm 10 that is "
                         "1,000 fits per run, which pins the reported false-positive "
                         "rate to about +/-1.4 points at a true 5%%.")
    ap.add_argument("--n-perm", type=int, default=10,
                    help="permutation replicates per run (default 10). Cost is linear: "
                         "each replicate refits --n-sample analytes. 1 reproduces the "
                         "pre-2026-08-15 single-shuffle behaviour.")
    ap.add_argument("--workers", type=int, default=1, metavar="N",
                    help="Fit analytes across N worker processes. LMM (trajectory) runs "
                         "only, matching regressions.py --workers. 1 (default) is the "
                         "serial path. 0 means auto: all cores but two, capped at 14. "
                         "Results are unaffected by the count.")
    ap.add_argument("--seed", type=int, default=20260731)
    ap.add_argument("--output", default="calibration_report.md")
    args = ap.parse_args()

    n_workers = args.workers
    if n_workers == 0:
        n_workers = max(1, min(14, (os.cpu_count() or 2) - 2))
    if n_workers < 0:
        ap.error("--workers must be >= 0")

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
        ctx = load_permutation_inputs(args.config)
        names = [s.get("name") for s in ctx["cfg"]["runs"] if s.get("name")]
        targets = names if args.permute == ["all"] else args.permute
        lines.append(f"Sweeping **{len(targets)} of {len(names)}** configured runs at "
                     f"stratum **{args.stratum}**, scope **{ctx['scope']}** "
                     f"({len(ctx['proteomic']):,} in-scope analytes), seed {args.seed}.\n")

        lines.append("| Run | Model | Shuffled | Term | perms | n fits | λ mean ± sd | λ range | "
                     "% P<0.05 mean ± sd | % range | warn |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        warn_total = 0
        skipped: list[tuple[str, str]] = []
        for run_name in targets:
            spec = next((s for s in ctx["cfg"]["runs"] if s.get("name") == run_name), {})
            model = spec.get("model", "?")
            res, note = permutation_null(ctx, run_name, args.stratum,
                                         args.n_sample, args.seed, args.n_perm, n_workers)
            if res.empty:
                lines.append(f"| `{run_name}` | {model} | — | — | 0 | 0 | — | — | — | — | skipped |")
                skipped.append((run_name, note))
                print(f"  {run_name}: SKIPPED — {note}")
                continue
            shuffled = "+".join(permutation_target(spec, ctx["defaults"], ctx["df"].columns)[0])
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
                    f"| `{run_name}` | {model} | `{shuffled}` | {term} | {len(per_perm)} | "
                    f"{int(per_perm['n'].sum()):,} | "
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

        # Listed explicitly rather than silently absent: a sweep that covered 40 of 52
        # runs and a sweep that covered 52 look identical in the table above, and the
        # difference is exactly what a reader would need to know before quoting a range.
        if skipped:
            lines.append(f"\n### Skipped ({len(skipped)} of {len(targets)})\n")
            lines.append("| Run | Why |")
            lines.append("|---|---|")
            for run_name, note in skipped:
                lines.append(f"| `{run_name}` | {note} |")
            lines.append("")

    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
