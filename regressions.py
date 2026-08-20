#!/usr/bin/env python3
"""regressions.py — generalized batch regressions over the unified PPMI dataset.

Runs one or more regression analyses (OLS / Logit / LMM / Cox) defined in a YAML
batch config. Each run loops a list of predictors (proteomic analytes, the 14
non-proteomic biomarkers, or both), fits one model per predictor per ancestry
stratum, and writes a timestamped CSV of summary statistics.

Covariates are the user-supplied `default_covariates` ± `add_covariates` /
`drop_covariates`. When the looped predictor is a proteomic column, the 5
assay-specific baseline PC columns are auto-appended to the RHS. Per-run
Bonferroni is 0.05 / N_fitted after numerical-failure filtering.

Usage:
    python3 regressions.py --config batch.yaml [--run run_name] [--resume]
                          [--workers N]

Design mirrors ../LRRK2-April2026/lrrk2_regressions.py for OLS/Logit/LMM math
(z-scoring, 3-tier LMM fallback, same numerical-failure filter, λ per run/assay);
Cox is added via lifelines.
"""

from __future__ import annotations

import argparse
import ctypes
import glob
import logging
import multiprocessing as mp
import os
import re
import sys
import warnings
from datetime import datetime
from typing import Any

# Pin BLAS to one thread per process, BEFORE numpy/scipy import — OpenBLAS reads
# these at load time, and a forked worker inherits whatever the parent set up.
#
# Every fit here is small (n ~ 2,600 rows x 13 parameters), far below the size at
# which threaded GEMM pays for its synchronisation: 16 threads measured 2.93 s/fit
# against 2.61 s/fit at one thread. The 2026-08-15 batch ran at ~2 of 16 cores for
# exactly this reason — OpenBLAS spilling onto a second core for no gain. Cores are
# worth far more spent on independent analytes (see --workers), and 14 workers each
# spawning 16 BLAS threads would oversubscribe the machine 14x.
#
# An explicit setting in the environment wins, so this can still be overridden.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

import statsmodels.formula.api as smf
from scipy.stats import chi2, zscore
from scipy.stats import t as student_t

try:
    from lifelines import CoxPHFitter
    from lifelines.exceptions import ConvergenceWarning as _LifelinesConvergenceWarning
    _HAS_LIFELINES = True
    # Silence expected within-stratum low-variance warnings from GP2 PCs
    # (PCs are computed within-ancestry upstream, so within a stratum their
    # scale collapses — does not affect β/SE/P).
    warnings.filterwarnings("ignore", category=_LifelinesConvergenceWarning)
except ImportError:
    CoxPHFitter = None  # type: ignore
    _HAS_LIFELINES = False

import yaml

# ============================================================================
# Constants
# ============================================================================

# Every proteomic analyte column must be covered by one of these prefixes:
# list_proteomic_columns() selects on the _NPX/_NPQ suffix, and proteomic_assay_map()
# then maps each column to its assay via these prefixes. A column matching no prefix
# raises KeyError when its assay PCs are injected.
ASSAY_PREFIXES = [
    "p277_CSF",
    "p282_CNS_CSF",
    "p282_Inflammation_CSF",
    "p288_CNS_plasma",
    "p288_Inflammation_plasma",
    "p293_olink_plasma",
    "p312_Inflammation_CSF",
    "p312_Inflammation_Plasma",
    "p312_Neuro_CSF",
    "p312_Neuro_Plasma",
    "p314_CSF",
    "p314_Plasma",
    # Harmonized blocks (phase 1 wrap). Each pools the two projects that ran the same
    # platform x panel x biofluid onto the larger project's scale. These are ADDITIONAL
    # to the project-specific prefixes above, not replacements — the project columns are
    # retained so any harmonized result stays checkable against the unpooled data.
    # A harmonized column contains its reference project's values, so the two are heavily
    # correlated and must NOT be counted as independent tests: see N_INDEPENDENT_PROTEINS.
    "harmonized_olink_plasma",
    "harmonized_olink_csf",
    "harmonized_nulisa_cns_plasma",
    "harmonized_nulisa_cns_csf",
    "harmonized_nulisa_inf_plasma",
    "harmonized_nulisa_inf_csf",
]

# The six harmonized blocks, in the same order as above. Each pools the two projects
# that ran the same platform x panel x biofluid onto the larger project's scale, so a
# harmonized column covers MORE PARTICIPANTS than either source project alone.
HARMONIZED_PREFIXES = [p for p in ASSAY_PREFIXES if p.startswith("harmonized_")]

# `defaults.proteomic_scope` selects which prefixes count as proteomic predictors:
#
#   all         every prefix above (34,902 analytes on the 2026-08-15 release)
#   harmonized  the six harmonized blocks only (11,567), plus the non-proteomic
#               predictors, which are unaffected either way
#
# `harmonized` is a superset of the source projects in PARTICIPANTS but a subset in
# ANALYTES: each block's core is the INTERSECTION of its two projects, so an analyte
# only one project measured has no harmonized counterpart. On the 2026-08-15 release
# that orphans 201 analytes -- 94 from p312_Neuro_Plasma and 84 from p312_Neuro_CSF
# (the Neuro panel is wider than the CNS panel it pairs with), 15 from
# p293_olink_plasma, and Abeta-38/40/42 from p288_CNS_plasma. Use `all` to keep them.
PROTEOMIC_SCOPES = ("all", "harmonized")

# Effective number of independent proteins across the whole panel set, used as the
# denominator for the proteome-wide Bonferroni lens. A raw analyte count is the wrong
# denominator: the same protein is measured by several projects and again in the
# harmonized blocks, so ~35k analyte columns represent far fewer independent tests.
# Overridable per-config via defaults.n_independent_proteins.
N_INDEPENDENT_PROTEINS = 9500

# Bonferroni family for every predictor without a high-throughput proteomics project
# prefix: the targeted biomarkers (abeta, asyn, tau, ptau, NFL, BMP species, urate),
# genetics (APOE_e4, p9005_Genetic_PRS_PRS157) and clinical measures (upsit, hemohi). Several of
# these are proteins too -- they are just measured by targeted/clinical assays rather
# than a high-throughput panel, which is why the family is named for the platform and
# not for whether the analyte is a protein.
NON_HT_PROTEOMICS_FAMILY = "non_highthroughput_proteomics"

# Non-proteomic biomarkers to analyze alongside the proteomic analytes.
# Extend this list to include new markers; can also be overridden per-config
# via defaults.new_predictors in the YAML.
NEW_PREDICTORS = [
    "APOE_e4",
    # Project 9005 (GP2 release 12) PRS, carried under the p9001_ prefix until
    # 2026-08-15. Replaces GP2_PRS_zscore, which it correlates with at r=0.992 while
    # covering 3,676 participants against 2,620 — so PRS results stay comparable to
    # earlier batches. PRS152/PRS149 in the same file are a different construct
    # (r~0.5) and are NOT interchangeable; see the Project 9005 documentation.
    "p9005_Genetic_PRS_PRS157",
    "upsit",
    "abeta",
    "asyn",
    "bd_tau_plasma",
    "hemohi",
    "NFL_CSF",       # dataset uses uppercase; recode cleanup normalized the dict to match
    "nfl_serum",
    "ptau",
    "tau",
    "total_di_18_1_BMP",
    "total_di_22_6_BMP",  # word doc listed 22:1 but dataset has 22:6 — the present di-BMP species
    "urate",
]

DEFAULT_COVARIATES = [
    "age_at_visit",
    "EDUCYRS",
    "SEX",
    "GP2_PC1",
    "GP2_PC2",
    "GP2_PC3",
    "GP2_PC4",
    "GP2_PC5",
]

PROTEIN_TOKEN = "PROTEIN"  # placeholder in formula / report_term

ALPHA = 0.05
MIN_N_DEFAULT = 20

# Assay control / spike-in analytes. These are NOT human proteins and cannot carry
# biological signal, so they double as a built-in negative control: their type-I rate
# is the empirical calibration of the whole pipeline. They are deliberately KEPT in
# the output (flagged, not dropped) so that calibration stays visible -- audit of
# 2026-07-31 found the trajectory LMMs rejecting these at 20% against a nominal 5%.
CONTROL_ANALYTE_RE = re.compile(r"(?:_CTRL_NPX|_mCherry_NPQ)$")

# A fitted coefficient this small is numerically zero, not a measurement. Every fit
# caught by this threshold in the 2026-07-30 batch was a control analyte (41/41, 0 real
# proteins), so it is a FLAG rather than a drop -- dropping would delete exactly the
# control rows that make the failure mode visible. Do NOT add a lower bound on SE
# instead: SE is scale-dependent, and an `SE < 0.1/sqrt(N)` rule deletes 24% of the
# legitimate slope_lowput_ratio fits (median |Z| 0.64, nominal 4.6% type-I).
DEGENERATE_BETA_THRESHOLD = 1e-10

# An LMM needs enough independent clusters to support its fixed effects. Requiring
# n_groups >= 2 x (number of fixed-effect parameters) eliminates saturated fits and
# drops the median Wald-z vs t p-value inflation from 20.5x to 1.9x while retaining
# 83% of fits; 3-4x buys little more and costs real EUR fits.
MIN_GROUPS_PER_PARAM_DEFAULT = 2

# The Cox analogue: a survival model is carried by its EVENTS, not its row count, so
# requiring events >= 2 x (number of design parameters) is the same rule as above
# expressed in the currency Cox actually spends. Set per-run/defaults via
# `min_events_per_param`; 0 disables the gate.
#
# 2 is deliberately permissive and matches MIN_GROUPS_PER_PARAM_DEFAULT. It removes
# what cannot be estimated, not everything underpowered: on the 2026-08-20 batch it
# drops the four cox_nsd_2a/2b_to_later|AJ stratum-runs (9-18 events against 15-17
# parameters, i.e. 0.5-1.2 events per parameter) and keeps 86% of proteomic rows.
# Raising it prunes hard -- 3 -> 68%, 5 -> 46%, 10 (the classic EPV rule) -> 12%, which
# would leave AJ represented by cox_pm_any alone. Rows carry `events_per_param` so a
# stricter cut can be applied downstream without refitting.
MIN_EVENTS_PER_PARAM_DEFAULT = 2

# Models whose analyte loop is worth forking across workers, and how many analytes to
# hand a worker at a time. LMM is the expensive one (a full REML optimisation per
# analyte, 13-17 h per stratum serially). Cox is far cheaper per fit (~46 ms) but there
# are 11,581 of them per stratum-run, which is still ~9 minutes each and ~5 h for a
# 32-run Cox batch; batching amortises the IPC that would otherwise dominate at that
# fit cost. OLS/Logit stay serial -- they finish in ~2 minutes, less than the fork and
# copy-on-write warm-up would cost.
PARALLEL_MODELS = ("lmm", "cox")
CHUNKSIZE_BY_MODEL = {"lmm": 1, "cox": 64}

# Collapsed-variance guard for LMM fits, on the SE * sqrt(N) scale. Because the LMM
# outcome is z-scored this quantity is comparable across analytes: the 2026-08-01
# batch had a median of 1.45 and a 0.1st percentile of 0.28, so 0.2 sits below
# essentially the entire observed distribution and flags only genuine collapses.
# LMM ONLY -- SE is scale-dependent, and applying a floor to the slope_* runs (which
# regress a raw clinical slope) would discard a quarter of their valid fits.
LOW_VARIANCE_SE_SCALED = 0.2

# Z-score defaults — False means the predictor is left on its native scale.
# `binary`:  {0, 1}-valued predictors (e.g. hemohi, any indicator flag).
# `dosage`:  integer {0, 1, 2}-valued predictors (e.g. APOE_e4 allele count)
#            which matches genotype-dosage convention in PD genetics work.
# Any predictor not matching those two patterns is z-scored (continuous +
# higher-cardinality ordinals like UPSIT 0–40).
ZSCORE_BINARY_DEFAULT = False
ZSCORE_DOSAGE_DEFAULT = False

# ============================================================================
# Logging
# ============================================================================

log = logging.getLogger(__name__)


def setup_logging(log_file: str) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file),
        ],
        force=True,
    )


# ============================================================================
# Config / data loading
# ============================================================================

def load_config(path: str) -> dict:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    if "runs" not in cfg or not isinstance(cfg["runs"], list):
        sys.exit(f"ERROR: config {path} must contain a 'runs' list")
    cfg.setdefault("defaults", {})
    return cfg


def resolve_input_path(input_glob: str, cwd: str) -> str:
    candidates = sorted(glob.glob(os.path.join(cwd, input_glob)))
    if not candidates:
        sys.exit(f"ERROR: no files match input_glob='{input_glob}' in {cwd}")
    return candidates[-1]


# Rows per chunk in the streaming load. 2,000 keeps the tokenizer's full-width buffer
# small while leaving few enough chunks that the concat is cheap; see load_data.
LOAD_CHUNK_ROWS = 2000


def _release_free_memory() -> None:
    """Ask glibc to return free arena memory to the OS.

    Python freeing a buffer does not shrink RSS -- glibc keeps it in the arena. After
    the load that is ~1.8 GB of held-but-unused memory, which every forked LMM worker
    would otherwise inherit. No-op on platforms without glibc's malloc_trim.
    """
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass


def read_header(path: str) -> list[str]:
    """Column names only, without reading a single data row."""
    with open(path) as fh:
        return fh.readline().rstrip("\n").split("\t")


def load_data(path: str, prefixes: list[str] | None = None) -> pd.DataFrame:
    """Read the release into memory, keeping only the in-scope analyte columns.

    Memory, not CPU, is the binding constraint: the release is ~19,451 rows x 35,566
    columns, of which 34,902 are analytes, and a full load is a ~5.6 GB frame with a
    peak well above that. That ceiling is what held the 2026-08-15 batch to a single
    process on a 31 GB machine while it used 2 of 16 cores.

    `prefixes` (from `defaults.proteomic_scope`) drops the analyte columns no run will
    look at. EVERY non-analyte column is kept unconditionally -- there are only 664 of
    them and they carry the clinical scaffold, covariates, outcomes, strata and the
    assay/genetic PCs, so narrowing by scope cannot cost a run a column it needs. Under
    scope=harmonized that is 12,231 columns instead of 35,566, a 2.91x cut.

    Read in row chunks and concatenated, which is what keeps the PEAK down: the C
    tokenizer buffers the full file width regardless of `usecols`, so a single-shot
    read spikes far above the resulting frame. Measured on the full release at
    scope=harmonized -- 19,450 rows x 12,231 columns, a 1.77 GB frame:

        single-shot   peak 14.30 GB   33 s
        chunked       peak  3.89 GB   23 s   <- this path

    The chunked frame is bit-identical to the single-shot one: same 12,231 columns,
    zero dtype differences, zero value or NaN-placement differences. (Chunked reads
    can otherwise infer dtypes per chunk, so this was verified, not assumed.)

    Finally `malloc_trim`, because glibc holds freed arena memory rather than
    returning it: without it the loaded process sits at 3.74 GB against a 1.77 GB
    frame, with it at 1.95 GB. That matters because the LMM worker pool forks from
    this process -- every GB resident here is a GB each worker inherits.

    NOT the pyarrow CSV engine, despite it being the obvious lever. Its default 1 MB
    `block_size` is smaller than ~9 rows of a file this wide (~120 KB per row), so it
    builds ~2,200 blocks x 35,566 columns of Arrow chunks, each with 64-byte-aligned
    buffers -- tens of GB of allocator padding on the full file. On a 2,000-row slice
    it measured 3.4x the peak memory of the C engine and 14x the wall time, and it is
    what took the WSL VM down on 2026-08-18.

    NOT float32: halving the analyte dtype would halve memory again, but it would also
    change fitted coefficients in the last significant figures, so runs from before and
    after the change could not be pooled. Dropping unread columns is free of that.
    """
    log.info(f"Loading dataset: {os.path.basename(path)}")
    usecols = None
    if prefixes is not None:
        cols = read_header(path)
        in_scope = set(list_proteomic_columns_from_names(cols, prefixes))
        usecols = [c for c in cols
                   if c in in_scope or not _is_analyte(c)]
        log.info(f"  reading {len(usecols):,} of {len(cols):,} columns "
                 f"({len(cols) - len(usecols):,} out-of-scope analytes skipped)")
    chunks = pd.read_csv(path, sep="\t", low_memory=False, usecols=usecols,
                         chunksize=LOAD_CHUNK_ROWS)
    df = pd.concat(list(chunks), ignore_index=True)
    del chunks
    _release_free_memory()
    log.info(f"  {df.shape[0]:,} rows x {df.shape[1]:,} cols "
             f"({df.memory_usage(deep=False).sum() / 1024**3:.2f} GB)")
    return df


# ============================================================================
# Predictor & covariate resolution
# ============================================================================

def scoped_assay_prefixes(scope: str) -> list[str]:
    """Assay prefixes a given `defaults.proteomic_scope` admits. See PROTEOMIC_SCOPES."""
    if scope == "all":
        return list(ASSAY_PREFIXES)
    if scope == "harmonized":
        return list(HARMONIZED_PREFIXES)
    sys.exit(f"ERROR: unknown proteomic_scope '{scope}'; expected one of {PROTEOMIC_SCOPES}")


def _is_analyte(col: str) -> bool:
    """The _NPX/_NPQ suffix is what makes a column a proteomic analyte."""
    return col.endswith("_NPX") or col.endswith("_NPQ")


def list_proteomic_columns_from_names(names: list[str],
                                      prefixes: list[str] | None = None) -> list[str]:
    """Analyte columns among `names`, optionally restricted to `prefixes`.

    Name-based so `load_data` can apply the same rule to a bare header line, before
    any frame exists -- the scope filter has to be decided before the read, not after.
    """
    cols = [c for c in names if _is_analyte(c)]
    if prefixes is None:
        return cols
    return [c for c in cols if any(c.startswith(p + "_") for p in prefixes)]


def list_proteomic_columns(df: pd.DataFrame,
                           prefixes: list[str] | None = None) -> list[str]:
    """Proteomic analyte columns of `df`, optionally restricted to `prefixes`.

    The prefix filter is what `proteomic_scope` uses to narrow which assays are in
    play. With prefixes=None every analyte column is returned (historical behaviour).
    """
    return list_proteomic_columns_from_names(list(df.columns), prefixes)


def proteomic_assay_map(cols: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for col in cols:
        for prefix in sorted(ASSAY_PREFIXES, key=len, reverse=True):
            if col.startswith(prefix + "_"):
                out[col] = prefix
                break
    return out


def resolve_predictors(df: pd.DataFrame, spec: Any, new_predictors: list[str],
                       prefixes: list[str] | None = None) -> list[str]:
    proteomic = list_proteomic_columns(df, prefixes)
    if spec is None:
        spec = "both"
    if isinstance(spec, str):
        key = spec.lower()
        if key == "proteomic":
            return proteomic
        if key == "new":
            return [c for c in new_predictors if c in df.columns]
        if key == "both":
            return proteomic + [c for c in new_predictors if c in df.columns]
        sys.exit(f"ERROR: unknown predictors spec '{spec}'")
    if isinstance(spec, list):
        missing = [c for c in spec if c not in df.columns]
        if missing:
            sys.exit(f"ERROR: predictor columns missing from dataset: {missing[:10]}")
        return list(spec)
    sys.exit(f"ERROR: predictors must be a string or list, got {type(spec).__name__}")


def resolve_covariates(
    defaults: list[str], add: list[str] | None, drop: list[str] | None
) -> list[str]:
    out = list(defaults)
    for c in add or []:
        if c not in out:
            out.append(c)
    for c in drop or []:
        if c in out:
            out.remove(c)
    return out


def assay_pcs(assay_prefix: str) -> list[str]:
    # Renamed from *_baseline_PC* on 2026-07-29: the PCA sample is each participant's
    # earliest visit carrying data for that assay, which is not necessarily BL.
    return [f"{assay_prefix}_earliest_visit_PC{i}" for i in range(1, 6)]


# ============================================================================
# Formula building
# ============================================================================

_VALID_IDENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def safe_col_name(col: str) -> str:
    return col if _VALID_IDENT_RE.match(col) else f"Q('{col}')"


def substitute_protein(formula: str, predictor_col: str) -> str:
    """Replace PROTEIN token with a safe-quoted column reference."""
    return re.sub(
        rf"\b{PROTEIN_TOKEN}\b",
        safe_col_name(predictor_col),
        formula,
    )


def substitute_protein_term(term: str, predictor_col: str) -> str:
    """Same substitution for a report_term (may be 'PROTEIN', 'PROTEIN:X', etc.)."""
    return re.sub(
        rf"\b{PROTEIN_TOKEN}\b",
        safe_col_name(predictor_col),
        term,
    )


def append_covariates(formula: str, covariates: list[str]) -> str:
    if not covariates:
        return formula
    cov_rhs = " + ".join(safe_col_name(c) for c in covariates)
    sep = " + " if "~" in formula else ""
    return formula + sep + cov_rhs


def build_full_formula(
    base_formula: str,
    predictor_col: str,
    covariates: list[str],
    proteomic_pcs: list[str],
) -> str:
    f = substitute_protein(base_formula, predictor_col)
    all_cov = covariates + [c for c in proteomic_pcs if c not in covariates]
    return append_covariates(f, all_cov)


# ============================================================================
# Model fitters — each returns (beta, SE, P, N, model_type, extras) or None
# ============================================================================

def _lookup_term(params_index: pd.Index, term: str) -> str | None:
    """Find a coefficient name in the model's index. Handles Q('x') wrapping."""
    if term in params_index:
        return term
    # If term is bare, try wrapped; if wrapped, try bare
    if term.startswith("Q('") and term.endswith("')") and term[3:-2] in params_index:
        return term[3:-2]
    candidate = f"Q('{term}')"
    if candidate in params_index:
        return candidate
    # Handle interaction terms: sides may each be wrapped/unwrapped independently
    if ":" in term:
        parts = term.split(":")
        variants = [":".join(p) for p in _interaction_variants(parts, params_index)]
        for v in variants:
            if v in params_index:
                return v
    return None


def _interaction_variants(parts: list[str], params_index: pd.Index) -> list[list[str]]:
    """Generate alternate orderings + Q-wrapping variants for interaction terms."""
    def variants_of(p: str) -> list[str]:
        return [p, f"Q('{p}')", p[3:-2] if p.startswith("Q('") and p.endswith("')") else p]
    from itertools import product, permutations
    out = []
    for perm in permutations(parts):
        for combo in product(*(variants_of(p) for p in perm)):
            out.append(list(combo))
    return out


def fit_ols(df: pd.DataFrame, formula: str, report_term: str) -> dict | None:
    try:
        fit = smf.ols(formula, data=df).fit()
    except Exception:
        return None
    key = _lookup_term(fit.params.index, report_term)
    if key is None:
        return None
    return {
        "beta": float(fit.params[key]),
        "SE": float(fit.bse[key]),
        "P": float(fit.pvalues[key]),
        "N": int(fit.nobs),
        "model_type": "OLS",
    }


def fit_logit(df: pd.DataFrame, formula: str, report_term: str) -> dict | None:
    try:
        fit = smf.logit(formula, data=df).fit(disp=0, maxiter=100)
    except Exception:
        return None
    key = _lookup_term(fit.params.index, report_term)
    if key is None:
        return None
    return {
        "beta": float(fit.params[key]),
        "SE": float(fit.bse[key]),
        "P": float(fit.pvalues[key]),
        "N": int(fit.nobs),
        "model_type": "LOGIT",
    }


def fit_lmm(
    df: pd.DataFrame,
    formula: str,
    report_terms: list[tuple[str, str]],
    random_group: str,
    *,
    random_slope: str | None = None,
    min_groups_per_param: int = MIN_GROUPS_PER_PARAM_DEFAULT,
) -> list[dict] | None:
    """LMM returning one row per reported term.

    `report_terms` is a list of (term_label, formula_term) pairs -- e.g.
    [("within", "YEAR_dev:grp_HC_vs_PD"), ("between", "YEAR_mean_c:grp_HC_vs_PD")].
    All terms come from ONE fit; this is two coefficients, not two models.

    Fallback ladder (the rung reached is recorded in `model_type` so it can be
    audited downstream):
        LMM_RS_{solver}  random slope on `random_slope`  -- correct spec for a
                         slope contrast; restores nominal type-I error (9.5% -> 5.3%
                         under a follow-up-structure-preserving permutation null)
        LMM_RI_{solver}  random intercept only -- where slope variance is not
                         estimable (shallow panels), this is all the data supports
        OLS_fallback     no mixed model converged

    Inference uses a t reference with df = n_groups - n_fixed_params rather than
    statsmodels' Wald-z. The normal approximation is anti-conservative at small
    cluster counts: on the 2026-07-30 batch the AJ trajectory hits had a median of
    29 subjects and a median p-value inflation of 24,312x against a t reference.
    """
    fit_result = None
    model_type = None
    re_singular = False

    ladder: list[tuple[str | None, str]] = []
    if random_slope:
        ladder.append((f"~{random_slope}", "RS"))
    ladder.append((None, "RI"))

    for re_formula, tag in ladder:
        for method in ("lbfgs", "powell"):
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    kw = {"re_formula": re_formula} if re_formula else {}
                    candidate = smf.mixedlm(
                        formula, data=df, groups=df[random_group], **kw
                    ).fit(reml=True, method=method)
                    singular = any(
                        "singular" in str(w.message).lower()
                        or "boundary" in str(w.message).lower()
                        for w in caught
                    )
                fit_result = candidate
                model_type = f"LMM_{tag}_{method}"
                re_singular = singular
                break
            except Exception:
                fit_result = None
        if fit_result is not None:
            break

    if fit_result is None:
        # OLS fallback on baseline rows only (matches LRRK2 script behavior)
        try:
            bl = df[df["EVENT_ID"] == "BL"].copy()
            if len(bl) < 5:
                return None
            fit_result = smf.ols(formula, data=bl).fit()
            model_type = "OLS_fallback"
        except Exception:
            return None

    is_lmm = model_type.startswith("LMM")
    params = fit_result.fe_params if is_lmm else fit_result.params
    bse = fit_result.bse_fe if is_lmm else fit_result.bse

    n_groups = int(fit_result.model.n_groups) if is_lmm else int(df[random_group].nunique())
    n_params = int(len(params))

    # Cluster-count gate: an LMM cannot support more fixed effects than it has
    # independent clusters to estimate them from.
    if is_lmm and n_groups < min_groups_per_param * n_params:
        return None

    resid_df = max(n_groups - n_params, 1)
    n_multi = int((df.groupby(random_group).size() >= 2).sum())

    out: list[dict] = []
    for term_label, formula_term in report_terms:
        key = _lookup_term(params.index, formula_term)
        if key is None:
            continue
        try:
            beta = float(params[key])
            se = float(bse[key])
        except Exception:
            continue
        if not np.isfinite(se) or se <= 0:
            continue
        p = float(2 * student_t.sf(abs(beta / se), resid_df))
        # SE on a 1/sqrt(N) scale. The LMM outcome is z-scored, so this is
        # comparable across analytes and a collapsed variance is visible as a value
        # far below the bulk. On the 2026-08-01 batch the median was 1.45 while the
        # four most extreme hits sat at 0.028-0.085 -- 17-50x below typical, and 86%
        # of fits under 0.1 had a singular random-effects covariance. Flagged, not
        # dropped, and deliberately NOT applied to non-LMM runs: the slope_* models
        # regress a raw clinical slope where small SEs are correct.
        se_scaled = float(se * np.sqrt(fit_result.nobs))
        out.append({
            "term": term_label,
            "beta": beta,
            "SE": se,
            "P": p,
            "N": int(fit_result.nobs),
            "n_groups": n_groups,
            "n_multi_obs": n_multi,
            "resid_df": resid_df,
            "re_singular": bool(re_singular),
            "se_scaled": se_scaled,
            "low_variance": bool(se_scaled < LOW_VARIANCE_SE_SCALED),
            "model_type": model_type,
        })
    return out or None


def fit_cox(
    df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    predictor_col: str,
    covariates: list[str],
    skip_zscore: bool = False,
    interaction_with: str | None = None,
    entry_col: str | None = None,
    min_events_per_param: float = 0.0,
) -> dict | None:
    """Cox PH via lifelines with 3-tier fallback.

    Tier 1: unpenalized CoxPHFitter (primary).                     -> COX
    Tier 2: L2-penalized fit, penalizer=0.01 (handles mild         -> COX_ridge_small
            collinearity / quasi-separation).
    Tier 3: L2-penalized fit, penalizer=0.1 (stronger regulariz-   -> COX_ridge_large
            ation; last resort before skipping).

    Penalized β is biased toward zero; the model_type column records
    which tier succeeded so these can be filtered or flagged
    downstream (e.g. in meta-analysis or visualization).

    If called directly, the predictor is z-scored here. When called from
    `_fit_one_predictor`, z-scoring (or the binary/dosage skip) has already
    happened upstream, so `skip_zscore=True` avoids a double transform.

    `entry_col` gives LEFT TRUNCATION (delayed entry). tte_* is measured from study
    baseline, but the predictor comes from each participant's earliest visit carrying
    it -- year 4 for Project 312. Without delayed entry those participants contribute
    immortal person-time, and the 2026-07-31 audit found 57-89% of p312 events
    occurring BEFORE the analyte was measured, i.e. post-event proteomics predicting
    the event. Participants not still at risk at entry are dropped: that loss is the
    honest answer about what a late-sampled panel can predict, not a defect.
    """
    if not _HAS_LIFELINES:
        sys.exit("ERROR: lifelines not installed. Run: pip install -r requirements.txt")
    if interaction_with == predictor_col:
        # Self-interaction (e.g. the PRS crossed with itself) is degenerate
        # — skip rather than emit a quadratic-term row.
        return None
    cols = [duration_col, event_col, predictor_col] + [c for c in covariates if c != predictor_col]
    if interaction_with and interaction_with not in cols:
        cols.append(interaction_with)
    if entry_col and entry_col in df.columns and entry_col not in cols:
        cols.append(entry_col)
    # Only keep cols that exist
    cols = [c for c in cols if c in df.columns]
    sub = df[cols].dropna().copy()

    n_pre_entry = 0
    use_entry = entry_col if (entry_col and entry_col in sub.columns) else None
    if use_entry:
        # Only participants still event-free at entry are in the risk set.
        at_risk = sub[duration_col] > sub[use_entry]
        n_pre_entry = int((~at_risk).sum())
        sub = sub.loc[at_risk].copy()
        if sub[use_entry].max() <= 0:
            use_entry = None  # every entry at time 0: plain Cox, no truncation needed

    if len(sub) < MIN_N_DEFAULT:
        return None
    if sub[event_col].sum() < 5:  # too few events
        return None

    if not skip_zscore:
        try:
            sub[predictor_col] = zscore(sub[predictor_col])
        except Exception:
            return None

    # Build interaction column inline so lifelines' design matrix sees it.
    # Reported term flips to predictor:interaction_with when requested.
    report_col = predictor_col
    if interaction_with:
        if interaction_with not in sub.columns:
            return None
        interaction_name = f"{predictor_col}:{interaction_with}"
        sub[interaction_name] = sub[predictor_col].astype(float) * sub[interaction_with].astype(float)
        report_col = interaction_name
        cols.append(interaction_name)

    # Drop all-constant columns that would break the fit. The entry column is a
    # timing variable, not a covariate, so it is preserved like duration/event.
    protected = {duration_col, event_col} | ({use_entry} if use_entry else set())
    keep = [c for c in cols if c in protected or sub[c].nunique(dropna=True) > 1]
    sub = sub[keep]
    if report_col not in sub.columns:
        return None

    # lifelines needs a fully numeric design matrix. A categorical covariate reaches
    # CoxPHFitter as-is -- `collection_era` is the string 'pre_2020'/'2020_plus' -- and
    # raises "could not convert string to float", which the tier loop below swallows,
    # so the analyte silently produces no row. Until 2026-08-20 the only Cox fits that
    # ever succeeded were those whose subsample happened to span a SINGLE era, where the
    # constant-column filter above removed the string before it reached the fit; those
    # results were era-unadjusted without recording it. Every proteomic analyte spans
    # both eras, so the whole proteomic side of the Cox arm was empty.
    #
    # drop_first keeps the dummies from being collinear with the baseline hazard, which
    # a Cox model has no intercept to absorb.
    cat_cols = [c for c in sub.columns
                if c not in protected and not pd.api.types.is_numeric_dtype(sub[c])]
    if cat_cols:
        sub = pd.get_dummies(sub, columns=cat_cols, drop_first=True, dtype=float)
        if report_col not in sub.columns:
            return None

    # Events-per-parameter gate, the Cox counterpart to the LMM cluster-count gate.
    # A Cox model's information comes from its events, so N is the wrong yardstick:
    # cox_nsd_2b_to_later|AJ reaches N=26 but only 18 events against 17 parameters,
    # and lifelines converges on it happily while the coefficients are meaningless.
    # Applied to the FINAL design, after constant columns are dropped and the dummies
    # are expanded, so the parameter count is the one actually estimated.
    n_params = int(sum(1 for c in sub.columns if c not in protected))
    n_events = int(sub[event_col].sum())
    if min_events_per_param > 0 and n_events < min_events_per_param * n_params:
        _FIT_STATUS["reason"] = "thin_events"
        return None

    fit_kwargs = {"duration_col": duration_col, "event_col": event_col, "show_progress": False}
    if use_entry:
        fit_kwargs["entry_col"] = use_entry

    fit_result = None
    model_type = None
    for penalizer, tag in [(0.0, "COX"), (0.01, "COX_ridge_small"), (0.1, "COX_ridge_large")]:
        try:
            cph = CoxPHFitter(penalizer=penalizer)
            cph.fit(sub, **fit_kwargs)
            if report_col in cph.summary.index:
                fit_result = cph
                model_type = tag
                break
        except Exception:
            continue

    if fit_result is None:
        return None

    if report_col not in fit_result.summary.index:
        return None
    row = fit_result.summary.loc[report_col]
    return {
        "beta": float(row["coef"]),
        "SE": float(row["se(coef)"]),
        "P": float(row["p"]),
        "N": int(sub.shape[0]),
        "n_events": n_events,
        # Design parameters actually estimated, and the events available per one of
        # them. Emitted so a stricter cut than `min_events_per_param` can be applied
        # downstream without refitting: the classic EPV rule is 10, the gate default
        # is 2, and everything in between is a filter on this column.
        "n_params": n_params,
        "events_per_param": round(n_events / n_params, 3) if n_params else np.nan,
        # Participants excluded because their event preceded the predictor
        # measurement. A large value means the panel was sampled too late to
        # predict this outcome from study baseline.
        "n_dropped_pre_entry": n_pre_entry,
        "delayed_entry": bool(use_entry),
        "model_type": model_type,
    }


# ============================================================================
# Per-predictor orchestration
# ============================================================================

def _predictor_kind(s: pd.Series) -> str:
    """Classify a numeric predictor: 'binary' ({0,1}), 'dosage' ({0,1,2}), 'other'."""
    vals = set(s.dropna().unique())
    if not vals:
        return "other"
    if vals <= {0.0, 1.0}:
        return "binary"
    if vals <= {0.0, 1.0, 2.0}:
        return "dosage"
    return "other"


def _should_zscore(kind: str, zscore_binary: bool, zscore_dosage: bool) -> bool:
    if kind == "binary":
        return zscore_binary
    if kind == "dosage":
        return zscore_dosage
    return True


def z_scored_copy(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c in out.columns and pd.api.types.is_numeric_dtype(out[c]):
            s = out[c].astype(float)
            if s.std(skipna=True) > 0:
                out[c] = zscore(s, nan_policy="omit")
    return out


def _earliest_visit_index(work: pd.DataFrame, pred: str, patno_values) -> pd.Index:
    """Row index of each participant's earliest visit where `pred` is non-null.

    `work` is pre-sorted by (PATNO, visit priority, EVENT_ID), so PATNO is
    non-decreasing and the first surviving row per participant is their earliest
    visit. Operates on two numpy arrays rather than the DataFrame, because this runs
    once per predictor (~23k times per run-stratum).
    """
    mask = work[pred].notna().to_numpy()
    if not mask.any():
        return work.index[:0]
    idx = work.index.to_numpy()[mask]
    _, first = np.unique(patno_values[mask], return_index=True)
    return pd.Index(idx[first])


# Read-only per-run state for the analyte loop. A forked worker inherits this by
# virtue of it being module-level and set before the pool is created, so the
# multi-GB frame is shared copy-on-write rather than pickled once per task.
# Never mutated by a worker.
_LOOP_STATE: dict[str, Any] = {}

# Why the most recent `_fit_one_predictor` call produced no rows. Written by that
# function, read once by `_fit_and_label`, which ships it back to the parent so the
# tally is identical on the serial and fork paths. Single-slot and reset per analyte.
_FIT_STATUS: dict[str, str] = {}

# Outcomes tallied per run/stratum and reported when the CSV is written. Without this,
# an analyte that never converged is indistinguishable from one that was never
# applicable: both simply produce no row. The 2026-08-15 NSD_vs_notNSD_LRRK2|EUR run
# lost 34,377 of 34,916 analytes this way and still logged "wrote 519 rows".
FIT_OUTCOMES = ("ok", "ok_no_assay_pcs", "below_min_n", "thin_events", "fit_failed")


def _fit_and_label(pred: str) -> dict:
    """Fit one analyte and attach its labelling/bookkeeping fields.

    The whole body of the analyte loop, factored out so the serial and parallel
    paths run *identical* code — the only difference being which process calls it.

    Takes only the analyte name: everything else is read from `_LOOP_STATE`, so a
    task costs a short string to ship rather than a copy of the design frame.

    Returns {"rows": [...], "status": one of FIT_OUTCOMES}. `rows` is empty for every
    status but `ok`/`ok_no_assay_pcs`.
    """
    st = _LOOP_STATE
    work = st["work"]
    spec = st["spec"]
    _FIT_STATUS.clear()

    is_proteomic = pred in st["proteomic_set"]
    # Auto-inject assay PCs for proteomic predictors
    pc_cols = assay_pcs(st["assay_of"][pred]) if is_proteomic else []

    # One row per participant at their earliest visit carrying this predictor.
    # Only the row INDEX is computed here -- subsetting the full 23,800-column
    # frame per predictor would copy ~1.5 GB each time. _fit_one_predictor
    # narrows to the ~20 design columns first, then applies this index.
    if st["visit_mode"] == "earliest_with_predictor":
        row_index = _earliest_visit_index(work, pred, st["patno_values"])
        if len(row_index) < st["min_n"]:
            return {"rows": [], "status": "below_min_n"}
    else:
        row_index = None

    rows = _fit_one_predictor(
        work,
        spec=spec,
        model=st["model"],
        predictor=pred,
        covariates=st["covariates"],
        proteomic_pcs=pc_cols,
        random_group=st["random_group"],
        min_n=st["min_n"],
        min_events_per_param=st["min_events_per_param"],
        zscore_binary=st["zscore_binary"],
        zscore_dosage=st["zscore_dosage"],
        row_index=row_index,
    )
    if not rows:
        return {"rows": [], "status": _FIT_STATUS.get("reason", "fit_failed")}

    pred_label, out_label = _predictor_outcome_labels(spec, pred)
    for row in rows:
        # For a decomposed LMM the reported term differs per row, so the
        # predictor label has to name which component it is.
        term = row.get("term")
        label = pred_label if term in (None, "effect") else f"{pred_label}[{term}]"
        row.update({
            "predictor": label,
            "outcome": out_label,
            # The looped analyte/column. This -- not the row count -- is the
            # Bonferroni denominator, since a decomposed LMM emits two rows
            # per analyte and correcting on rows would silently double it.
            "analyte": pred,
            "term": term or "effect",
            "strata": st["strata_val"] or "ALL",
            "is_proteomic": bool(is_proteomic),
            # Assay panel this predictor belongs to, used as the Bonferroni family.
            "assay_family": st["assay_of"][pred] if is_proteomic else NON_HT_PROTEOMICS_FAMILY,
            # Spike-in / assay controls: kept deliberately so their type-I rate
            # is a visible calibration check, flagged so they are never mistaken
            # for findings.
            "is_control": bool(CONTROL_ANALYTE_RE.search(pred)),
        })
    return {"rows": rows, "status": _FIT_STATUS.get("reason", "ok")}


def _analyte_results(predictors: list[str], n_workers: int, label: str,
                     strata_val: str | None, chunksize: int = 1):
    """Yield (i, result) over `predictors`, serially or across a fork pool.

    `result` is the dict `_fit_and_label` returns: {"rows": [...], "status": ...}.

    `imap` preserves input order, so the emitted sequence — and therefore the row
    order of the output CSV — is identical to the serial path regardless of worker
    count.

    `chunksize` is per-model (CHUNKSIZE_BY_MODEL). An LMM fit takes seconds and varies
    several-fold in cost (1.5-6.8 s observed), so chunks of 1 keep the tail balanced
    and IPC is negligible against the work. A Cox fit takes ~46 ms, where per-task
    round-trips would dominate instead, so its workers take larger batches.
    """
    if n_workers <= 1:
        for i, pred in enumerate(predictors):
            yield i, _fit_and_label(pred)
        return

    log.info(f"  [{label}|{strata_val}] fitting on {n_workers} worker processes "
             f"(chunksize={chunksize})")
    # fork (not spawn): workers must inherit `_LOOP_STATE` — and with it the frame —
    # through copy-on-write. Under spawn they would re-import and re-load it each.
    ctx = mp.get_context("fork")
    with ctx.Pool(n_workers) as pool:
        for i, result in enumerate(pool.imap(_fit_and_label, predictors, chunksize=chunksize)):
            yield i, result


def run_one(
    df: pd.DataFrame,
    spec: dict,
    strata_val: str | None,
    defaults: dict,
    proteomic_set: set[str],
    assay_of: dict[str, str],
    n_workers: int = 1,
    assay_prefixes: list[str] | None = None,
) -> pd.DataFrame:
    """Run one named spec on one strata subset; return results DataFrame."""
    label = spec.get("name", "unnamed")
    model = spec["model"].lower()

    # Stratum subset
    strata_col = defaults.get("strata_col")
    if strata_val is not None and strata_col:
        work = df[df[strata_col] == strata_val].copy()
    else:
        work = df.copy()

    # Sample filter (pandas query)
    sample_filter = spec.get("sample_filter")
    if sample_filter:
        try:
            work = work.query(sample_filter)
        except Exception as e:
            log.warning(f"  [{label}|{strata_val}] sample_filter '{sample_filter}' failed: {e}")
            return pd.DataFrame()

    if len(work) < MIN_N_DEFAULT:
        log.warning(f"  [{label}|{strata_val}] too few rows ({len(work)}) after sample_filter; skip")
        return pd.DataFrame()

    # Visit selection. "earliest_with_predictor" keeps one row per participant: their
    # earliest visit at which THAT predictor was actually measured. Restricting to
    # EVENT_ID == 'BL' instead would discard nearly all of Project 312, which was
    # sampled at V08-V16 (6-16 participants at BL vs 169-266 at their first visit).
    # Same rule as the *_earliest_visit_PC blocks: BL first, then SC, then visit order.
    # Sort once here so the per-predictor step is a mask + drop_duplicates, not a sort.
    visit_mode = spec.get("visit_mode", defaults.get("visit_mode", "as_filtered"))
    if visit_mode == "earliest_with_predictor":
        _pri = work["EVENT_ID"].map(lambda e: {"BL": 0, "SC": 1}.get(str(e), 2))
        work = (work.assign(_visit_priority=_pri)
                    .sort_values(["PATNO", "_visit_priority", "EVENT_ID"])
                    .drop(columns="_visit_priority"))
    elif visit_mode != "as_filtered":
        log.warning(f"  [{label}|{strata_val}] unknown visit_mode '{visit_mode}'; "
                    f"treating as 'as_filtered'")
        visit_mode = "as_filtered"

    # Predictors
    pred_spec = spec.get("predictors", defaults.get("predictors", "both"))
    new_preds = defaults.get("new_predictors", NEW_PREDICTORS)
    predictors = resolve_predictors(work, pred_spec, new_preds, assay_prefixes)
    log.info(f"  [{label}|{strata_val}] {len(predictors)} predictors to fit")

    # Covariates
    base_cov = spec.get("covariates", defaults.get("default_covariates", DEFAULT_COVARIATES))
    add_cov = spec.get("add_covariates", [])
    drop_cov = spec.get("drop_covariates", [])
    covariates = resolve_covariates(base_cov, add_cov, drop_cov)

    # Random group (LMM)
    random_group = spec.get("random_group", defaults.get("random_group", "PATNO"))
    min_n = spec.get("min_n", defaults.get("min_n", MIN_N_DEFAULT))
    min_epp = spec.get("min_events_per_param",
                       defaults.get("min_events_per_param", MIN_EVENTS_PER_PARAM_DEFAULT))
    zscore_binary = spec.get("zscore_binary", defaults.get("zscore_binary", ZSCORE_BINARY_DEFAULT))
    zscore_dosage = spec.get("zscore_dosage", defaults.get("zscore_dosage", ZSCORE_DOSAGE_DEFAULT))

    _patno_values = (work["PATNO"].to_numpy()
                     if visit_mode == "earliest_with_predictor" else None)

    # Publish the loop's read-only inputs before any pool is forked, so workers
    # inherit them rather than receiving a pickled copy per task.
    _LOOP_STATE.update(
        work=work, spec=spec, model=model, covariates=covariates,
        proteomic_set=proteomic_set, assay_of=assay_of, random_group=random_group,
        min_n=min_n, min_events_per_param=min_epp,
        zscore_binary=zscore_binary, zscore_dosage=zscore_dosage,
        visit_mode=visit_mode, patno_values=_patno_values, strata_val=strata_val,
    )

    results: list[dict] = []
    first_logged = False
    tally = {k: 0 for k in FIT_OUTCOMES}

    try:
        for i, result in _analyte_results(predictors, n_workers, label, strata_val,
                                          chunksize=CHUNKSIZE_BY_MODEL.get(model, 1)):
            rows = result["rows"]
            tally[result["status"]] = tally.get(result["status"], 0) + 1
            if not rows:
                continue
            results.extend(rows)

            if not first_logged or (i + 1) % 500 == 0:
                r0 = rows[0]
                log.info(
                    f"    [{i+1}/{len(predictors)}] {rows[0]['analyte']}: β={r0['beta']:+.4g} "
                    f"SE={r0['SE']:.4g} P={r0['P']:.3g} N={r0['N']}"
                )
                first_logged = True
    finally:
        # Drop the frame reference so the next run's `work` is the only copy held.
        _LOOP_STATE.clear()

    # Account for every predictor offered. A fit that never converged and one that was
    # never applicable both produce no row, so without this they are indistinguishable
    # and a run that lost most of its analytes still looks like it succeeded.
    log.info(f"  [{label}|{strata_val}] fit outcomes over {len(predictors):,} analytes: "
             + ", ".join(f"{k}={tally.get(k, 0):,}" for k in FIT_OUTCOMES))
    if tally.get("fit_failed"):
        pct = 100.0 * tally["fit_failed"] / max(1, len(predictors))
        log.warning(f"  [{label}|{strata_val}] {tally['fit_failed']:,} analytes ({pct:.1f}%) "
                    f"failed to fit even without assay PCs — check stratum size vs "
                    f"parameter count")
    if tally.get("ok_no_assay_pcs"):
        log.warning(f"  [{label}|{strata_val}] {tally['ok_no_assay_pcs']:,} analytes fitted "
                    f"only after dropping assay PCs (flagged `assay_pcs_dropped`)")
    if tally.get("thin_events"):
        pct = 100.0 * tally["thin_events"] / max(1, len(predictors))
        log.warning(f"  [{label}|{strata_val}] {tally['thin_events']:,} analytes ({pct:.1f}%) "
                    f"gated by min_events_per_param={min_epp} — too few events to support "
                    f"the design, NOT a fit failure")

    return pd.DataFrame(results)


def _predictor_outcome_labels(spec: dict, looped_col: str) -> tuple[str, str]:
    """Return (predictor_label, outcome_label) for the output row.

    Baseline-style runs (formula LHS is a fixed column, not PROTEIN):
        predictor = looped column, outcome = formula LHS
    Trajectory-style runs (formula LHS is PROTEIN):
        predictor = substituted report_term (e.g. 'YEAR:grp_NSD_vs_HC')
        outcome   = looped column (which IS the LHS after substitution)
    Cox: predictor = looped column, outcome = duration_col.
    """
    model = spec["model"].lower()
    if model == "cox":
        # Cox: predictor column = looped analyte (matches OLS interaction convention).
        # When interaction_with is set the reported β/SE/P is for the interaction
        # term, but the predictor column stays as the analyte name — the run's
        # name suffix (e.g. `_x_PRS`) tells the reader it's an interaction.
        return looped_col, spec.get("duration_col", "")
    formula = spec.get("formula", "")
    if "~" not in formula:
        return looped_col, spec.get("outcome", "")
    lhs = formula.split("~", 1)[0].strip()
    if lhs == PROTEIN_TOKEN:
        report_term = spec.get("report_term", PROTEIN_TOKEN)
        # Substitute PROTEIN in report_term with looped col (unquoted, human-readable label)
        pred_label = re.sub(rf"\b{PROTEIN_TOKEN}\b", looped_col, report_term)
        return pred_label, looped_col
    return looped_col, lhs


def decompose_time_spec(
    report_term: str, time_col: str
) -> tuple[str, list[tuple[str, str]]] | None:
    """Split a `TIME:GROUP` report term into its within- and between-subject parts.

    Returns (rhs_core, report_terms) where rhs_core is the replacement for the
    `TIME * GROUP` block and report_terms pairs a human label with its formula term.

    Rationale: with one coefficient on `TIME:GROUP`, subjects contribute both
    within-subject change and cross-sectional between-subject differences to the
    same number, and the mix depends on panel depth. Measured on the 2026-07-30
    batch, the as-run coefficient correlated 0.98 with the pure within term in
    p293_olink_plasma (94% of subjects have a 2nd visit) but only 0.49 in p314_CSF
    (22%), where it tracked the pure between term at 0.85 instead. Splitting the
    term keeps every subject -- singletons have TIME_dev == 0 so they inform the
    intercept, group main effect, covariates and residual variance but contribute
    nothing to the trajectory term, which is exactly right.
    """
    parts = [p.strip() for p in report_term.split(":")]
    if len(parts) != 2 or time_col not in parts:
        return None
    group = parts[1] if parts[0] == time_col else parts[0]
    dev, mean_c = f"{time_col}_dev", f"{time_col}_mean_c"
    rhs_core = f"{dev} * {group} + {mean_c} * {group}"
    return rhs_core, [("within", f"{dev}:{group}"), ("between", f"{mean_c}:{group}")]


def _add_time_components(sub: pd.DataFrame, time_col: str, random_group: str) -> pd.DataFrame:
    """Attach subject-mean (between) and within-subject-deviation (within) time terms."""
    grp_mean = sub.groupby(random_group)[time_col].transform("mean")
    sub[f"{time_col}_mean_c"] = grp_mean - grp_mean.mean()  # centered: grp term stays interpretable
    sub[f"{time_col}_dev"] = sub[time_col] - grp_mean       # 0 for single-visit subjects
    return sub


def _fit_one_predictor(
    work: pd.DataFrame,
    *,
    spec: dict,
    model: str,
    predictor: str,
    covariates: list[str],
    proteomic_pcs: list[str],
    random_group: str,
    min_n: int,
    min_events_per_param: float = 0.0,
    zscore_binary: bool,
    zscore_dosage: bool,
    row_index: pd.Index | None = None,
) -> list[dict] | None:
    time_col = spec.get("time_col", "YEAR")
    # Figure out all columns the design needs (for dropna + sample-size guard)
    if model == "cox":
        duration_col = spec["duration_col"]
        event_col = spec["event_col"]
        all_cov = covariates + proteomic_pcs
        needed = [duration_col, event_col, predictor] + all_cov
        # Entry time for left truncation: the visit the predictor was measured at.
        if spec.get("delayed_entry", True) and time_col in work.columns:
            needed.append(time_col)
    else:
        formula_tmpl = spec.get("formula")
        if not formula_tmpl:
            sys.exit(f"ERROR: spec '{spec.get('name')}' has model={model} but no 'formula'")
        # LHS column (if not PROTEIN)
        lhs_raw = formula_tmpl.split("~", 1)[0].strip()
        lhs = predictor if lhs_raw == PROTEIN_TOKEN else lhs_raw
        # Terms needed from RHS (strip PROTEIN placeholder -> predictor)
        rhs = formula_tmpl.split("~", 1)[1]
        rhs_cols = _columns_in_formula(rhs, predictor, covariates, proteomic_pcs, random_group, model)
        needed = [lhs] + rhs_cols
    needed = [c for c in needed if c in work.columns]
    # Narrow columns BEFORE selecting rows. `work` can be ~23,800 columns wide, so
    # row-masking it first would copy the whole frame once per predictor.
    sub = work[list(dict.fromkeys(needed))]
    if row_index is not None:
        sub = sub.loc[row_index]
    sub = sub.dropna().copy()
    if len(sub) < min_n:
        _FIT_STATUS["reason"] = "below_min_n"
        return None

    # Z-score the looped predictor (matches LRRK2 style), skipping binary {0,1}
    # and {0,1,2} genotype-dosage predictors unless explicitly enabled.
    if predictor in sub.columns and pd.api.types.is_numeric_dtype(sub[predictor]):
        s = sub[predictor].astype(float)
        kind = _predictor_kind(s)
        if _should_zscore(kind, zscore_binary, zscore_dosage) and s.std(skipna=True) > 0:
            sub[predictor] = zscore(s, nan_policy="omit")

    if model == "cox":
        # Optional Cox PRS-style interaction: report β/SE/P for the
        # `predictor:interaction_with` term. The interaction column is
        # constructed inline (predictor * interaction_with) so it's available
        # in lifelines' design matrix without requiring formula support.
        interaction_with = spec.get("interaction_with")
        entry_col = time_col if (spec.get("delayed_entry", True)
                                 and time_col in sub.columns) else None
        row = fit_cox(
            sub,
            duration_col=spec["duration_col"],
            event_col=spec["event_col"],
            predictor_col=predictor,
            covariates=[c for c in covariates + proteomic_pcs if c in sub.columns],
            skip_zscore=True,  # _fit_one_predictor already z-scored (or not) per switch
            interaction_with=interaction_with,
            entry_col=entry_col,
            min_events_per_param=min_events_per_param,
        )
        return [row] if row else None

    # Build formula
    formula = build_full_formula(spec["formula"], predictor, covariates, proteomic_pcs)
    report_term = substitute_protein_term(spec.get("report_term", PROTEIN_TOKEN), predictor)

    if model in ("ols", "logit"):
        fitter = fit_ols if model == "ols" else fit_logit
        row = fitter(sub, formula, report_term)
        # Assay-PC fallback, mirroring the 3-tier LMM fallback below. The 5 auto-injected
        # assay PCs are 5 of ~16 parameters, and in a small stratum that is enough to push
        # a logit into quasi-complete separation: the IRLS Hessian goes singular even
        # though the design matrix is full rank, and fit_logit's bare `except` turns that
        # into a silently missing analyte. On the 2026-08-15 NSD_vs_notNSD_LRRK2|EUR run
        # it destroyed all four Olink blocks -- 0/30 analytes fitted at N=46, against
        # 29/30 on the SAME rows with the PCs dropped. Refitting without them keeps the
        # analyte instead of losing it; `assay_pcs_dropped` marks which rows are
        # unadjusted so they stay distinguishable from the adjusted ones.
        dropped = False
        if row is None and proteomic_pcs:
            row = fitter(sub, build_full_formula(spec["formula"], predictor, covariates, []),
                         report_term)
            dropped = row is not None
        if row is None:
            _FIT_STATUS["reason"] = "fit_failed"
            return None
        row["assay_pcs_dropped"] = dropped
        _FIT_STATUS["reason"] = "ok_no_assay_pcs" if dropped else "ok"
        return [row]
    if model == "lmm":
        if random_group not in sub.columns:
            log.warning(f"  random_group '{random_group}' missing for {predictor}; skip")
            return None

        report_terms = [("effect", report_term)]
        random_slope = spec.get("random_slope")
        if spec.get("decompose_time") and time_col in sub.columns:
            split = decompose_time_spec(report_term, time_col)
            if split is None:
                log.warning(f"  [{spec.get('name')}] decompose_time set but report_term "
                            f"'{report_term}' is not TIME:GROUP; fitting undecomposed")
            else:
                rhs_core, report_terms = split
                sub = _add_time_components(sub, time_col, random_group)
                # Swap the `TIME * GROUP` block for its within/between expansion.
                lhs, rhs = formula.split("~", 1)
                original_block = report_term.replace(":", " * ")
                flipped = " * ".join(reversed(report_term.split(":")))
                for block in (original_block, flipped):
                    if block in rhs:
                        rhs = rhs.replace(block, rhs_core, 1)
                        break
                else:
                    log.warning(f"  [{spec.get('name')}] could not locate '{original_block}' "
                                f"in formula; fitting undecomposed")
                    report_terms = [("effect", report_term)]
                    rhs = None
                if rhs is not None:
                    formula = f"{lhs}~{rhs}"
                    if random_slope == time_col:
                        random_slope = f"{time_col}_dev"

        min_gpp = spec.get("min_groups_per_param", MIN_GROUPS_PER_PARAM_DEFAULT)
        return fit_lmm(
            sub, formula, report_terms, random_group,
            random_slope=random_slope,
            min_groups_per_param=min_gpp,
        )
    sys.exit(f"ERROR: unknown model '{model}'")


def _columns_in_formula(
    rhs: str,
    predictor: str,
    covariates: list[str],
    proteomic_pcs: list[str],
    random_group: str,
    model: str,
) -> list[str]:
    """Best-effort enumeration of dataframe columns referenced by the RHS."""
    cols: list[str] = [predictor] + list(covariates) + list(proteomic_pcs)
    if model == "lmm":
        cols.append(random_group)
    # Pull out bare identifiers from the RHS (e.g. 'YEAR', 'grp_NSD_vs_HC')
    for tok in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", rhs):
        if tok in (PROTEIN_TOKEN, "Q"):
            continue
        cols.append(tok)
    return list(dict.fromkeys(cols))


# ============================================================================
# Post-fit filtering, Bonferroni, lambda
# ============================================================================

def filter_numerical_failures(df: pd.DataFrame, label: str) -> pd.DataFrame:
    """Drop unusable fits, and FLAG (never drop) numerically-degenerate ones.

    Note there is deliberately no lower bound on SE. SE is scale-dependent: the
    slope_* runs regress a raw clinical slope, so small SEs there are correct, and
    an `SE < 0.1/sqrt(N)` rule would delete 24% of the slope_lowput_ratio fits
    despite a median |Z| of 0.64 and a nominal 4.6% type-I rate. The scale-free
    signal of a collapsed fit is the coefficient itself.
    """
    if df.empty:
        return df
    before = len(df)
    bad = (
        ~np.isfinite(df["beta"])
        | ~np.isfinite(df["SE"])
        | ~np.isfinite(df["P"])
        | (df["beta"].abs() > 20)
        | (df["SE"] > 20)
        | (df["SE"] <= 0)
        | (df["P"] <= 0)
        | (df["P"] > 1)
    )
    n_bad = int(bad.sum())
    if n_bad:
        log.warning(f"  {label}: dropped {n_bad}/{before} numerical failures")
    out = df.loc[~bad].copy()
    if out.empty:
        return out
    out["degenerate"] = out["beta"].abs() < DEGENERATE_BETA_THRESHOLD
    n_deg = int(out["degenerate"].sum())
    if n_deg:
        n_ctrl = int((out["degenerate"] & out.get("is_control", False)).sum())
        log.warning(f"  {label}: flagged {n_deg} numerically-degenerate fits "
                    f"(|β| < {DEGENERATE_BETA_THRESHOLD:g}; {n_ctrl} are control analytes) "
                    f"— retained and flagged, not dropped")
    return out


def apply_bonferroni(df: pd.DataFrame, alpha: float = ALPHA,
                     n_independent: int = N_INDEPENDENT_PROTEINS) -> pd.DataFrame:
    """Attach three Bonferroni lenses to a run's results.

    `bonferroni_threshold` / `significant` -- the run-wide correction over every
    predictor fitted in this run. This is the original behaviour and is what
    meta_analysis.py, plot_results.py and the summary's headline counts use.

    `bonferroni_threshold_family` / `significant_family` -- corrected within each
    assay panel (`assay_family`) instead. Adding Projects 312/314 roughly doubled the
    run-wide denominator, so hits from the pre-existing panels can lose run-wide
    significance purely because more analytes are now being tested alongside them.
    The family lens keeps each panel judged against its own number of tests, which
    stays comparable with earlier cuts.

    `bonferroni_threshold_proteome` / `significant_proteome` -- corrected against a
    FIXED denominator of `n_independent` (default 9,500), the effective number of
    independent proteins across the panel set. This is the lens to use once the
    harmonized blocks are in play: the same protein is measured by several projects and
    again in its harmonized block, so a raw analyte count massively overstates the number
    of independent tests, while the run-wide lens would penalise every panel for the
    duplication. A fixed independent-protein denominator is invariant to how many
    redundant columns a run happens to include.
    """
    if df.empty:
        return df
    # Denominator is the number of ANALYTES tested, not the number of result rows.
    # A decomposed LMM emits two rows per analyte (within + between); dividing by
    # rows would double the correction as a side effect of the output schema.
    # The within and between terms are separate questions, each corrected over the
    # analytes tested — so the family-wise error rate across both is nominally 2·α.
    key = "analyte" if "analyte" in df.columns else None
    n = int(df[key].nunique()) if key else len(df)
    thresh = alpha / n if n > 0 else np.nan
    df["n_analytes_tested"] = n
    df["bonferroni_threshold"] = thresh
    df["significant"] = df["P"] < thresh

    if "assay_family" in df.columns:
        if key:
            fam_n = df.groupby("assay_family")[key].transform("nunique")
        else:
            fam_n = df.groupby("assay_family")["P"].transform("size")
        df["n_analytes_tested_family"] = fam_n
        df["bonferroni_threshold_family"] = alpha / fam_n
        df["significant_family"] = df["P"] < df["bonferroni_threshold_family"]
    else:
        df["n_analytes_tested_family"] = n
        df["bonferroni_threshold_family"] = thresh
        df["significant_family"] = df["significant"]

    # Proteome-wide lens: fixed independent-protein denominator, applied only to
    # high-throughput proteomic predictors. The targeted/clinical biomarkers in
    # NON_HT_PROTEOMICS_FAMILY are not part of the proteome panel set, so they keep the
    # family threshold rather than being judged against 9,500 protein tests.
    df["n_independent_proteins"] = n_independent
    prot_thresh = alpha / n_independent if n_independent > 0 else np.nan
    if "assay_family" in df.columns:
        is_prot = df["assay_family"] != NON_HT_PROTEOMICS_FAMILY
        df["bonferroni_threshold_proteome"] = np.where(
            is_prot, prot_thresh, df["bonferroni_threshold_family"])
    else:
        df["bonferroni_threshold_proteome"] = prot_thresh
    df["significant_proteome"] = df["P"] < df["bonferroni_threshold_proteome"]
    return df


def lambda_gc(pvals: pd.Series) -> tuple[float, int]:
    # Genomic inflation factor. Inputs MUST be 2-sided p-values — this holds
    # by default for statsmodels (OLS/Logit/LMM) and lifelines (Cox). The
    # transform chi2.ppf(1-p, df=1) converts a 2-sided p to a 1-df chi-square.
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

def find_existing_output(cwd: str, output_dir: str, name: str, suffix: str) -> str | None:
    """Newest non-empty CSV already written for this run and stratum, or None.

    Outputs are named f"{name}-{suffix}-{timestamp}.csv", so globbing the literal
    f"{name}-{suffix}-" prefix cannot collide with a longer run name that shares a
    prefix: "slope_moca-EUR-*" does not match "slope_moca_x_PRS-EUR-...".

    Any matching file counts, whatever batch wrote it — the timestamp is not
    compared against the current run.
    """
    pattern = os.path.join(cwd, output_dir, f"{name}-{suffix}-*.csv")
    hits = [p for p in sorted(glob.glob(pattern)) if os.path.getsize(p) > 0]
    return hits[-1] if hits else None


def run_all(config: dict, cwd: str, run_filter: list[str] | None,
            resume: bool = False, n_workers: int = 1) -> None:
    defaults = config.get("defaults", {})
    output_dir = defaults.get("output_dir", "results")
    os.makedirs(os.path.join(cwd, output_dir), exist_ok=True)

    strata_col = defaults.get("strata_col")
    strata_vals = defaults.get("strata") or [None]

    # Resolve the work list before loading the dataset: a resumed batch that has
    # nothing left to do should not spend minutes reading a multi-GB .tab first.
    work: list[tuple[dict, list[tuple[str | None, str]]]] = []
    n_skipped = 0
    for spec in config["runs"]:
        name = spec.get("name", "unnamed")
        if run_filter and name not in run_filter:
            continue

        run_strata = spec.get("strata", strata_vals)
        if not run_strata:
            run_strata = [None]

        pending: list[tuple[str | None, str]] = []
        for sv in run_strata:
            suffix = sv if sv is not None else "ALL"
            done = find_existing_output(cwd, output_dir, name, suffix) if resume else None
            if done:
                log.info(f"  resume: {name}|{suffix} already written "
                         f"({os.path.basename(done)}); skipping")
                n_skipped += 1
            else:
                pending.append((sv, suffix))
        if pending:
            work.append((spec, pending))

    n_pending = sum(len(p) for _, p in work)
    if resume:
        log.info(f"Resume: {n_skipped} run/stratum pair(s) already on disk, "
                 f"{n_pending} to fit.")
    if not work:
        log.info("Nothing to do.")
        return

    input_glob = defaults.get("input_glob", "unified_PPMI-*.tab")
    input_path = resolve_input_path(input_glob, cwd)

    # Resolve the scope BEFORE the read, so the out-of-scope analyte columns are never
    # loaded at all rather than loaded and then filtered. See load_data: the full-width
    # read is what puts peak memory over what a 31 GB machine can hold.
    scope = defaults.get("proteomic_scope", "all")
    assay_prefixes = scoped_assay_prefixes(scope)
    log.info(f"  proteomic_scope={scope} ({len(assay_prefixes)} assay blocks)")
    df = load_data(input_path, assay_prefixes)

    # Precompute proteomic/assay map once, restricted to the configured scope.
    proteomic_cols = list_proteomic_columns(df, assay_prefixes)
    assay_of = proteomic_assay_map(proteomic_cols)
    proteomic_set = set(proteomic_cols)
    log.info(f"  {len(proteomic_cols):,} proteomic columns; {len(NEW_PREDICTORS)} non-proteomic in default list")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    for spec, pending in work:
        name = spec.get("name", "unnamed")
        log.info("")
        log.info("=" * 70)
        log.info(f"RUN: {name}  (model={spec['model']})")
        log.info("=" * 70)

        # Parallelise the LMM runs only. They are the whole cost of a batch — a
        # trajectory stratum takes 13-17 h against minutes for a logit/OLS/Cox run,
        # because each analyte is a full REML optimisation (99.4% of fit time, with
        # design construction at 0.6%). Every other model keeps the serial path
        # untouched, where the pool would only add fork overhead.
        run_workers = n_workers if spec["model"].lower() in PARALLEL_MODELS else 1

        for sv, suffix in pending:
            results = run_one(df, spec, sv, defaults, proteomic_set, assay_of,
                              n_workers=run_workers, assay_prefixes=assay_prefixes)
            results = filter_numerical_failures(results, f"{name}|{sv}")
            results = apply_bonferroni(
                results,
                n_independent=int(defaults.get("n_independent_proteins",
                                               N_INDEPENDENT_PROTEINS)))

            if results.empty:
                log.warning(f"  [{name}|{sv}] no results after filtering; no CSV written")
                continue

            out = os.path.join(cwd, output_dir, f"{name}-{suffix}-{timestamp}.csv")
            results.to_csv(out, index=False)
            n_sig = int(results["significant"].sum())
            thr = results["bonferroni_threshold"].iloc[0]
            n_sig_fam = int(results["significant_family"].sum())
            n_analytes = int(results["n_analytes_tested"].iloc[0])
            log.info(
                f"  [{name}|{sv}] wrote {len(results)} rows over {n_analytes} analytes "
                f"({n_sig} sig @ {thr:.2e} run-wide; "
                f"{n_sig_fam} sig within assay family) -> {os.path.basename(out)}"
            )

            # Lambda GC per reported term. Within- and between-subject terms are
            # differently calibrated (1.14 vs 0.70 in the same panel on the
            # 2026-07-30 batch), so a pooled λ hides both.
            for term, tdf in results.groupby("term", sort=False):
                lam_all, n_all = lambda_gc(tdf["P"])
                log.info(f"    [{term}] λ (overall) = {lam_all:.4f}  N_tests = {n_all}")
                if "is_proteomic" in tdf.columns and tdf["is_proteomic"].any():
                    prot_results = tdf[tdf["is_proteomic"]]
                    for prefix in ASSAY_PREFIXES:
                        mask = prot_results["analyte"].astype(str).str.startswith(prefix + "_")
                        if not mask.any():
                            continue
                        lam_a, n_a = lambda_gc(prot_results.loc[mask, "P"])
                        log.info(f"      {prefix}: λ = {lam_a:.4f}  N_tests = {n_a}")

            # Negative-control calibration: these analytes cannot carry signal, so
            # their type-I rate is the empirical calibration of this run.
            if "is_control" in results.columns and results["is_control"].any():
                for term, tdf in results.groupby("term", sort=False):
                    c = tdf[tdf["is_control"] & ~tdf.get("degenerate", False)]
                    if len(c) < 5:
                        continue
                    rate = 100 * float((c["P"] < 0.05).mean())
                    flag = "  <-- CHECK: >2x nominal" if rate > 10 else ""
                    log.info(f"    [{term}] negative controls: {len(c)} analytes, "
                             f"{rate:.1f}% at P<0.05 (nominal 5%){flag}")


# ============================================================================
# CLI
# ============================================================================

def main() -> None:
    p = argparse.ArgumentParser(description="Batch regressions on unified PPMI")
    p.add_argument("--config", required=True, help="Path to YAML config file")
    p.add_argument("--run", nargs="+", default=None,
                   help="Only execute the named run(s) (by 'name' key). Accepts multiple.")
    p.add_argument("--resume", action="store_true",
                   help="Skip any run/stratum that already has a non-empty CSV in "
                        "output_dir, from any batch. Composes with --run. A stratum "
                        "that legitimately wrote no CSV (every fit a numerical "
                        "failure, as small strata often produce) is indistinguishable "
                        "from one that died before writing, so it is re-fitted.")
    p.add_argument("--workers", type=int, default=1, metavar="N",
                   help="Fit analytes across N worker processes. Applies to LMM "
                        "(trajectory) runs only — every other model is fast enough "
                        "serially. 1 (default) is the serial path. 0 means auto: "
                        "all cores but two, capped at 14. Results are unaffected: "
                        "fits are independent and output order is preserved, so any "
                        "worker count yields a byte-identical CSV.")
    args = p.parse_args()

    n_workers = args.workers
    if n_workers == 0:
        n_workers = max(1, min(14, (os.cpu_count() or 2) - 2))
    if n_workers < 0:
        p.error("--workers must be >= 0")

    cwd = os.path.dirname(os.path.abspath(args.config)) or os.getcwd()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(cwd, f"regressions_{timestamp}.log")
    setup_logging(log_file)

    log.info(f"Config: {args.config}")
    log.info(f"Log:    {log_file}")
    if n_workers > 1:
        log.info(f"Workers: {n_workers} (LMM runs only; 1 BLAS thread each)")
    cfg = load_config(args.config)
    run_all(cfg, cwd, args.run, args.resume, n_workers)
    log.info("\nDone.")


if __name__ == "__main__":
    main()
