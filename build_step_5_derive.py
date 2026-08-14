"""Step 5 — derived analysis variables.

A port of ../recode_analysis_groups.py (group flags, CI flags, DaT lowput ratio,
slopes, time-to-event) plus the two LEDD-step derivations
(time_to_LEDD_years, disease_duration_years) from
../adding_p312_and_P314/unified-PPMI-dataset-proteomics-clinical-PRS-master/LEDD/ledd_analysis.py.

Reads the scaffold rather than the merged file — every derived variable is a function
of clinical columns only, so there is no reason to load 1.2 GB of proteomics to compute
them. Emits a `key`-indexed block that step 7 concatenates.

TWO DELIBERATE DEPARTURES from the reference implementation, both forced by the
refreshed clinical cut. Each is logged at runtime and recorded in the dictionary:

  1. `lowput_ratio` is taken directly from `MIA_LOWPUT_EXPECTED`, which is already the
     observed/expected ratio. This both switches series (the curated cut drops all
     DATSCAN_* columns; MIA is the MIAKAT pipeline, a different quantification, not a
     rename) and FIXES A PRE-EXISTING BUG.

     The reference implementation computed `min(DATSCAN_PUTAMEN_L, R) / lowput_expected`.
     But `lowput_expected` is documented in the source dictionary as the "age-/sex-
     expected lowest putamen *ratio*" — i.e. it is already observed/expected. Dividing
     the observed value by it a second time cancels the imaging term and leaves a
     function of age and sex. Evidence in the previous release:

       - corr(lowput_ratio, observed min-putamen SBR) = 0.055 — no imaging signal
       - lowput_ratio has 0 of 4,312 values below 0.75, while its own dictionary entry
         says "Ratios <0.75 indicate DaT positivity"; `lowput_expected` has 75.8%
       - slope_lowput_ratio is near-constant (IQR -0.0158..-0.0141, sd 0.002),
         tracking the -0.01397/yr age coefficient rather than any participant's DaT

     Consequence: `lowput_ratio` and `slope_lowput_ratio` differ substantially from the
     previous release by design, and any prior result using them should be re-read.

  2. `grp_NMC_*` excludes multi-gene carriers by gene-token test rather than by
     exact-matching the single label "LRRK2 + GBA". The refreshed cut introduces 12
     new `subgroup` values including "LRRK2 + GBA + Normosmic" and
     "LRRK2 + GBA + SNCA + PRKN + Normosmic"; under the original exact-match rule those
     would be counted as clean single-gene carriers in several NMC groups at once.
     The stated intent — "keep singletons clean" — is preserved, and behaviour on the
     old vocabulary is unchanged (gene+phenotype labels like "LRRK2 + RBD" still match).

Emits:
  build_intermediates/derived-<ts>.tab
  build_intermediates/derived_summary-<ts>.tab
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

from build_common import TIMESTAMP, log

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")

VISIT_YEAR_COL = "YEAR"
STAGE_EARLY = {"1a", "1b", "2a", "2b"}
STAGE_LATE = {"3", "4", "5", "6"}

# Tokens in `subgroup` that denote a causal/risk gene. Everything else in the field
# ("Normosmic", "RBD", "Hyposmia", ...) is a prodromal phenotype, not a second gene.
GENE_TOKENS = ["LRRK2", "GBA", "SNCA", "PRKN", "PINK1", "VPS35", "PARK7", "Other GV"]


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} — run the prior step first")
    return hits[-1]


def baseline_mask(df: pd.DataFrame) -> pd.Series:
    return df["EVENT_ID"] == "BL"


def write_patno_value(df: pd.DataFrame, name: str, by_patno: pd.Series) -> None:
    df[name] = df["PATNO"].map(by_patno)


# ---------------------------------------------------------------- group flags

def build_group_flags(df: pd.DataFrame) -> list[str]:
    added: list[str] = []
    bl_idx = df.loc[baseline_mask(df)].drop_duplicates(subset=["PATNO"]).set_index("PATNO")

    nsd = bl_idx["NSD_Status"].astype("float")
    sub = bl_idx["subgroup"].astype("string").fillna("")
    ana = bl_idx["analytic_subgroup"].astype("string").fillna("")
    stage = bl_idx["NSD_STAGE"].astype("string").fillna("")
    cohort = bl_idx["COHORT"].astype("string").fillna("")

    val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
    val[nsd == 1.0] = 1.0
    val[(sub == "Healthy Control") & (nsd == 0.0)] = 0.0
    write_patno_value(df, "grp_NSD_vs_HC", val); added.append("grp_NSD_vs_HC")

    for name, mask in [
        ("grp_NSD_vs_notNSD_sPD",       sub == "Sporadic PD"),
        ("grp_NSD_vs_notNSD_LRRK2",     sub == "LRRK2"),
        ("grp_NSD_vs_notNSD_GBA",       sub == "GBA"),
        ("grp_NSD_vs_notNSD_LRRK2GBA",  sub == "LRRK2 + GBA"),
        ("grp_NSD_vs_notNSD_SNCA",      sub == "SNCA"),
        ("grp_NSD_vs_notNSD_PRKN",      sub == "PRKN"),
        ("grp_NSD_vs_notNSD_prodromal", cohort == "4"),
    ]:
        val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
        in_grp = mask & nsd.notna()
        val[in_grp & (nsd == 1.0)] = 1.0
        val[in_grp & (nsd == 0.0)] = 0.0
        write_patno_value(df, name, val); added.append(name)

    for name, ref, comp in [
        ("grp_sPD_vs_LRRK2", "Sporadic PD", "LRRK2"),
        ("grp_sPD_vs_GBA",   "Sporadic PD", "GBA"),
        ("grp_LRRK2_vs_GBA", "LRRK2",       "GBA"),
    ]:
        val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
        val[sub == ref] = 0.0
        val[sub == comp] = 1.0
        write_patno_value(df, name, val); added.append(name)

    for name, ref, comp in [
        ("grp_HC_vs_Prodromal", "2", "4"),
        ("grp_HC_vs_PD",        "2", "1"),
        ("grp_Prodromal_vs_PD", "4", "1"),
    ]:
        val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
        val[cohort == ref] = 0.0
        val[cohort == comp] = 1.0
        write_patno_value(df, name, val); added.append(name)

    hc_neg = (cohort == "2") & (nsd == 0.0)
    prod_pos = (cohort == "4") & (nsd == 1.0)
    pd_pos = (cohort == "1") & (nsd == 1.0)
    for name, ref_mask, comp_mask in [
        ("grp_HCNSDneg_vs_ProdNSDpos", hc_neg, prod_pos),
        ("grp_HCNSDneg_vs_PDNSDpos",   hc_neg, pd_pos),
        ("grp_ProdNSDpos_vs_PDNSDpos", prod_pos, pd_pos),
    ]:
        val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
        val[ref_mask] = 0.0
        val[comp_mask] = 1.0
        write_patno_value(df, name, val); added.append(name)

    def stage_contrast(name: str, low: set[str], high: set[str]) -> None:
        val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
        val[stage.isin(low)] = 0.0
        val[stage.isin(high)] = 1.0
        write_patno_value(df, name, val); added.append(name)

    stage_contrast("grp_NSD_stage_2A_vs_2B", {"2a"}, {"2b"})
    stage_contrast("grp_NSD_stage_2A_vs_3",  {"2a"}, {"3"})
    stage_contrast("grp_NSD_stage_2B_vs_3",  {"2b"}, {"3"})
    stage_contrast("grp_NSD_stage_early_vs_late", STAGE_EARLY, STAGE_LATE)

    for abc in ["A1", "A2", "B", "C"]:
        val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
        val[ana.isin(["A1", "A2", "B", "C"])] = 0.0
        val[ana == abc] = 1.0
        write_patno_value(df, f"grp_ABC_{abc}", val); added.append(f"grp_ABC_{abc}")

    # --- NMC: single-gene carriers in the prodromal cohort --------------------
    def gene_hits(s: pd.Series) -> pd.DataFrame:
        return pd.DataFrame({g: s.str.contains(g, regex=False) for g in GENE_TOKENS})

    hits = gene_hits(sub)
    n_genes = hits.sum(axis=1)
    in_prod = cohort == "4"

    for gene in ["LRRK2", "GBA", "LRRK2GBA", "SNCA", "PRKN"]:
        val = pd.Series(np.nan, index=bl_idx.index, dtype="float")
        val[in_prod] = 0.0
        if gene == "LRRK2GBA":
            carrier = sub == "LRRK2 + GBA"
        else:
            # target gene present and no other gene alongside it
            carrier = hits[gene] & (n_genes == 1)
        val[in_prod & carrier] = 1.0
        write_patno_value(df, f"grp_NMC_{gene}", val); added.append(f"grp_NMC_{gene}")

    multi = sorted(set(sub[(n_genes > 1) & (sub != "LRRK2 + GBA")]))
    if multi:
        log(f"  NMC multi-gene labels excluded from single-gene groups: {multi}")

    return added


# ------------------------------------------------------------------ CI / DaT

def build_ci_flags(df: pd.DataFrame) -> list[str]:
    cog = pd.to_numeric(df["cogstate"], errors="coerce")
    df["CI_PI"] = np.where(cog.notna(), (cog > 1).astype(float), np.nan)
    moca = pd.to_numeric(df["moca"], errors="coerce")
    df["CI_MOCA"] = np.where(moca.notna(), (moca < 26).astype(float), np.nan)
    return ["CI_PI", "CI_MOCA"]


def build_dat_lowput_ratio(df: pd.DataFrame) -> list[str]:
    """Observed / age-and-sex-expected lowest-putamen SBR, per visit.

    `MIA_LOWPUT_EXPECTED` is, despite its name, *already* that ratio. The curated cut's
    dictionary gives the formula explicitly:

        Min(putamen_r_ref_cwm, putamen_l_ref_cwm)
        / (1.4474 - 0.003780*age_at_DATSCAN + 0.2093*gender)

    so it is taken as-is. See module docstring, departure 1: the reference
    implementation divided the observed putamen SBR by this quantity a second time,
    which cancels the imaging term and leaves a function of age and sex.
    """
    src = pd.to_numeric(df.get("MIA_LOWPUT_EXPECTED"), errors="coerce")
    df["lowput_ratio"] = src
    r = df["lowput_ratio"]
    obs = pd.concat([pd.to_numeric(df.get("MIA_PUTAMEN_L"), errors="coerce"),
                     pd.to_numeric(df.get("MIA_PUTAMEN_R"), errors="coerce")], axis=1).min(axis=1)
    m = r.notna() & obs.notna()
    log(f"  lowput_ratio: {r.notna().sum():,} rows, median {r.median():.3f}, "
        f"{(r < 0.75).sum():,} ({(r < 0.75).sum() / r.notna().sum() * 100:.1f}%) "
        "below the 0.75 DaT-positive mark")
    log(f"    corr with observed min-putamen SBR: {r[m].corr(obs[m]):.3f}  "
        "(previous release: 0.055 — the signal the double division removed)")
    return ["lowput_ratio"]


# --------------------------------------------------------------------- slopes

def _slope(years: np.ndarray, y: np.ndarray) -> float:
    m = np.isfinite(years) & np.isfinite(y)
    if m.sum() < 2:
        return np.nan
    x, yy = years[m], y[m]
    if np.ptp(x) == 0:
        return np.nan
    xm, ym = x.mean(), yy.mean()
    denom = ((x - xm) ** 2).sum()
    if denom == 0:
        return np.nan
    return float(((x - xm) * (yy - ym)).sum() / denom)


def build_slopes(df: pd.DataFrame) -> list[str]:
    added: list[str] = []
    year = pd.to_numeric(df[VISIT_YEAR_COL], errors="coerce")
    for out_name, src in {
        "slope_moca": "moca",
        "slope_updrs3_off": "updrs3_score",
        "slope_lowput_ratio": "lowput_ratio",
    }.items():
        if src not in df.columns:
            log(f"  WARNING: source '{src}' not found; skipping {out_name}")
            continue
        tmp = pd.DataFrame({"PATNO": df["PATNO"],
                            "y": pd.to_numeric(df[src], errors="coerce"), "year": year})
        slopes = {p: _slope(g["year"].to_numpy(float), g["y"].to_numpy(float))
                  for p, g in tmp.groupby("PATNO", sort=False)}
        write_patno_value(df, out_name, pd.Series(slopes))
        log(f"  {out_name}: {df[out_name].notna().sum():,} rows "
            f"({pd.Series(slopes).notna().sum():,} patients)")
        added.append(out_name)
    return added


# ------------------------------------------------------------- time-to-event

def _emit_tte(df: pd.DataFrame, name_prefix: str, tte: dict, ev: dict) -> list[str]:
    tte_name, ev_name = f"tte_{name_prefix}_years", f"event_{name_prefix}"
    write_patno_value(df, tte_name, pd.Series(tte))
    write_patno_value(df, ev_name, pd.Series(ev))
    n_ev = int(sum(v == 1.0 for v in ev.values()))
    log(f"  {name_prefix}: {len(ev):,} at-risk | {n_ev:,} events | {len(ev) - n_ev:,} censored")
    return [tte_name, ev_name]


def build_tte_threshold(df, *, name_prefix, source_col, event_condition, baseline_eligible):
    if source_col not in df.columns:
        log(f"  WARNING: '{source_col}' not in df; skipping {name_prefix}")
        return []
    y = pd.to_numeric(df[source_col], errors="coerce")
    year = pd.to_numeric(df[VISIT_YEAR_COL], errors="coerce")
    tmp = pd.DataFrame({"PATNO": df["PATNO"], "y": y, "year": year, "ev": event_condition(y)})
    tte, ev = {}, {}
    for pat, g in tmp.groupby("PATNO", sort=False):
        g = g.sort_values("year")
        bl = g[g["year"] == 0]
        if bl.empty or pd.isna(bl["y"].iloc[0]) or not baseline_eligible(float(bl["y"].iloc[0])):
            continue
        post = g[g["year"] > 0].dropna(subset=["y"])
        if post.empty:
            continue
        hit = post[post["ev"] == True]  # noqa: E712
        tte[pat], ev[pat] = ((float(hit["year"].iloc[0]), 1.0) if not hit.empty
                             else (float(post["year"].iloc[-1]), 0.0))
    return _emit_tte(df, name_prefix, tte, ev)


def build_tte_indicator(df, *, name_prefix, source_col):
    if source_col not in df.columns:
        log(f"  WARNING: '{source_col}' not in df; skipping {name_prefix}")
        return []
    y = pd.to_numeric(df[source_col], errors="coerce")
    year = pd.to_numeric(df[VISIT_YEAR_COL], errors="coerce")
    tmp = pd.DataFrame({"PATNO": df["PATNO"], "y": y, "year": year})
    tte, ev = {}, {}
    for pat, g in tmp.groupby("PATNO", sort=False):
        g = g.sort_values("year")
        bl = g[g["year"] == 0]
        if bl.empty or pd.isna(bl["y"].iloc[0]) or float(bl["y"].iloc[0]) != 0.0:
            continue
        post = g[g["year"] > 0].dropna(subset=["y"])
        if post.empty:
            continue
        hit = post[post["y"] == 1.0]
        tte[pat], ev[pat] = ((float(hit["year"].iloc[0]), 1.0) if not hit.empty
                             else (float(post["year"].iloc[-1]), 0.0))
    return _emit_tte(df, name_prefix, tte, ev)


def build_tte_nsd_stage(df, *, name_prefix, bl_stages, event_stages):
    year = pd.to_numeric(df[VISIT_YEAR_COL], errors="coerce")
    stage = df["NSD_STAGE"].astype("string").fillna("")
    tmp = pd.DataFrame({"PATNO": df["PATNO"], "stage": stage, "year": year})
    tte, ev = {}, {}
    for pat, g in tmp.groupby("PATNO", sort=False):
        g = g.sort_values("year")
        bl = g[g["year"] == 0]
        if bl.empty or bl["stage"].iloc[0] not in bl_stages:
            continue
        post = g[(g["year"] > 0) & (g["stage"] != "")]
        if post.empty:
            continue
        hit = post[post["stage"].isin(event_stages)]
        tte[pat], ev[pat] = ((float(hit["year"].iloc[0]), 1.0) if not hit.empty
                             else (float(post["year"].iloc[-1]), 0.0))
    return _emit_tte(df, name_prefix, tte, ev)


# ------------------------------------------------------------------ LEDD step

def build_ledd_and_duration(df: pd.DataFrame) -> list[str]:
    """Ported from LEDD/ledd_analysis.py. `time_to_LEDD_years` is written on BL rows
    only — not propagated — matching the reference implementation exactly."""
    dt = pd.to_datetime(df["visit_date"], errors="coerce")
    ledd = pd.to_numeric(df["LEDD"], errors="coerce")
    tmp = pd.DataFrame({"PATNO": df["PATNO"], "EVENT_ID": df["EVENT_ID"],
                        "dt": dt, "LEDD": ledd})

    df["time_to_LEDD_years"] = np.nan
    for _, g in tmp.groupby("PATNO", sort=False):
        bl = g[g["EVENT_ID"] == "BL"]
        if bl.empty or bl["dt"].isna().all():
            continue
        on = g[(g["LEDD"].notna()) & (g["LEDD"] != 0) & (g["dt"].notna())]
        if on.empty:
            continue
        first = on.sort_values("dt")["dt"].iloc[0]
        df.loc[bl.index, "time_to_LEDD_years"] = (first - bl["dt"].iloc[0]).days / 365.25

    age_v = pd.to_numeric(df["age_at_visit"], errors="coerce")
    agediag = pd.to_numeric(df["agediag"], errors="coerce")
    df["disease_duration_years"] = np.where(age_v.notna() & agediag.notna(),
                                            age_v - agediag, np.nan)

    log(f"  time_to_LEDD_years: {df['time_to_LEDD_years'].notna().sum():,} BL rows")
    log(f"  disease_duration_years: {df['disease_duration_years'].notna().sum():,} rows")
    return ["time_to_LEDD_years", "disease_duration_years"]


# ---------------------------------------------------------------------- main

def main() -> None:
    log(f"=== Step 5: derived variables ===  timestamp {TIMESTAMP}")
    df = pd.read_csv(latest("scaffold-*.tab"), sep="\t", low_memory=False,
                     dtype={"key": str, "PATNO": str, "EVENT_ID": str, "COHORT": str,
                            "NSD_STAGE": str, "subgroup": str, "analytic_subgroup": str})
    log(f"scaffold: {df.shape[0]:,} rows x {df.shape[1]} columns")

    log("\n--- LEDD / disease duration ---")
    ledd_cols = build_ledd_and_duration(df)

    log("\n--- Group flags ---")
    group_cols = build_group_flags(df)
    log(f"  added {len(group_cols)} group flags")

    log("\n--- CI flags ---")
    ci_cols = build_ci_flags(df)

    log("\n--- DaT lowput ratio ---")
    ratio_cols = build_dat_lowput_ratio(df)

    log("\n--- Slopes ---")
    slope_cols = build_slopes(df)

    log("\n--- Time-to-event ---")
    tte_cols: list[str] = []
    tte_cols += build_tte_threshold(df, name_prefix="moca_lt26", source_col="moca",
                                   event_condition=lambda y: y < 26,
                                   baseline_eligible=lambda v: v >= 26)
    tte_cols += build_tte_threshold(df, name_prefix="cogstate_worsen", source_col="cogstate",
                                    event_condition=lambda y: y > 1,
                                    baseline_eligible=lambda v: v == 1)
    for pfx, src in [("pm_any", "pm_any"), ("pm_cog_any", "pm_cog_any"),
                     ("pm_mc_any", "pm_mc_any"), ("stage_d", "Stage_D")]:
        tte_cols += build_tte_indicator(df, name_prefix=pfx, source_col=src)
    tte_cols += build_tte_nsd_stage(df, name_prefix="nsd_2a_to_later",
                                    bl_stages={"2a"}, event_stages={"2b", "3", "4", "5", "6"})
    tte_cols += build_tte_nsd_stage(df, name_prefix="nsd_2b_to_later",
                                    bl_stages={"2b"}, event_stages={"3", "4", "5", "6"})

    derived = ledd_cols + group_cols + ci_cols + ratio_cols + slope_cols + tte_cols
    log(f"\ntotal derived columns: {len(derived)}")

    out = df[["key"] + derived]
    f_out = os.path.join(OUT, f"derived-{TIMESTAMP}.tab")
    out.to_csv(f_out, sep="\t", index=False)

    summary = pd.DataFrame([
        dict(variable=c, n_nonnull=int(df[c].notna().sum()),
             n_patients=int(df.loc[df[c].notna(), "PATNO"].nunique()),
             n_positive=int((pd.to_numeric(df[c], errors="coerce") == 1).sum())
             if c.startswith(("grp_", "event_", "CI_")) else -1)
        for c in derived])
    summary.to_csv(os.path.join(OUT, f"derived_summary-{TIMESTAMP}.tab"), sep="\t", index=False)

    log(f"wrote {os.path.basename(f_out)}  ({out.shape[0]:,} x {out.shape[1]})")


if __name__ == "__main__":
    sys.exit(main())
