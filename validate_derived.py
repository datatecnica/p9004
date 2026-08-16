"""Validate the clinical and derived layers against the previous release.

The proteomic layer is checked exhaustively by verify_bit_identical.py. This does the
same job for the two layers that are *expected* to move, and separates legitimate
movement from drift that would indicate a defect.

Compared on the keys present in both releases:

  A. CLINICAL source columns carried over from the previous cut. These come from a
     different data cut (public curated vs internal proteomics), so disagreement is
     possible and is worth quantifying rather than assuming either way.

  A2. HARMONIZED proteomic columns, when the baseline has them. Reported as its own
     layer because it is ~50x the size of the clinical layer: pooling the two gives a
     percentage dominated by proteomics while labelled clinical. Against a
     pre-harmonization baseline this section is skipped, since the columns are not
     shared. NOTE when comparing runs: the clinical rate quoted against a
     pre-harmonization baseline and against a harmonized one were, before 2026-08-15,
     the same label over different column sets.

  B. DERIVED columns. Split by whether the participant's visit set changed, because 569
     of 2,025 previous participants gained a visit in the refreshed cut. Slopes and
     time-to-event are functions of a participant's whole visit history, so for those
     participants a changed value is correct, not a bug. For participants whose visit
     set is unchanged, a per-visit or baseline-determined variable should reproduce
     exactly unless its definition or its clinical inputs changed.

Three variables are expected to differ everywhere and are reported separately:
lowput_ratio and slope_lowput_ratio (bug fix + series change) and grp_NMC_* (definition).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

from build_common import BASELINE_DATASET_STEM, DATASET_STEM, require_build

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = require_build(DATASET_STEM)
OLD = require_build(BASELINE_DATASET_STEM)

DERIVED_PREFIXES = ("grp_", "tte_", "event_", "slope_", "CI_")
DERIVED_EXTRA = {"lowput_ratio", "time_to_LEDD_years", "disease_duration_years"}
EXPECTED_DIFFERENT = {"lowput_ratio", "slope_lowput_ratio"}

VISIT_DEPENDENT = ("slope_", "tte_", "event_")


def compare(a: pd.Series, b: pd.Series) -> tuple[int, int, int, int]:
    """-> (both_null, agree, value_differs, presence_differs)"""
    na, nb = a.isna(), b.isna()
    both_null = int((na & nb).sum())
    presence = int((na ^ nb).sum())
    both = ~na & ~nb
    if both.sum() == 0:
        return both_null, 0, 0, presence
    x, y = a[both], b[both]
    if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
        eq = np.isclose(x.astype(float), y.astype(float), rtol=0, atol=1e-9, equal_nan=True)
    else:
        eq = x.astype(str).str.strip().to_numpy() == y.astype(str).str.strip().to_numpy()
    return both_null, int(eq.sum()), int((~eq).sum()), presence


def main() -> None:
    nh = pd.read_csv(NEW, sep="\t", nrows=0).columns.tolist()
    oh = pd.read_csv(OLD, sep="\t", nrows=0).columns.tolist()
    oset = set(oh)

    derived = [c for c in nh
               if (c.startswith(DERIVED_PREFIXES) or c in DERIVED_EXTRA) and c in oset]
    # `harmonized_*` is proteomic. The project-prefix test below ("p" + digit) matches
    # p277_, p314_ and so on but NOT harmonized_, so before 2026-08-15 the whole
    # harmonized block fell through into `clinical`. That was invisible while the
    # baseline was the pre-harmonization release — those columns were not in `oset`, so
    # they were filtered out as unshared — and only surfaced once the baseline became a
    # release that HAS them, at which point "clinical" silently became 98% proteomics
    # (11,573 of 11,789 columns) and its agreement rate stopped being a clinical figure.
    proteomic = [c for c in nh
                 if (c.split("_")[0].startswith("p") and c[1].isdigit())
                 or c.startswith("harmonized_")]
    harmonized = [c for c in nh if c.startswith("harmonized_") and c in oset
                  and "earliest_visit_PC" not in c]
    clinical = [c for c in nh
                if c in oset and c not in derived and c not in set(proteomic)
                and not c.startswith(("GP2_", "p9005_", "PATNO_", "EVENT_ID_"))
                and "earliest_visit_PC" not in c
                and c not in ("key", "MERGE_INDEX")]

    print(f"clinical columns shared: {len(clinical)}   "
          f"harmonized shared: {len(harmonized)}   derived shared: {len(derived)}")

    a = pd.read_csv(NEW, sep="\t",
                    usecols=["key", "PATNO"] + clinical + harmonized + derived,
                    low_memory=False, dtype={"key": str, "PATNO": str})
    b = pd.read_csv(OLD, sep="\t",
                    usecols=["PATNO", "EVENT_ID"] + clinical + harmonized + derived,
                    low_memory=False, dtype={"PATNO": str, "EVENT_ID": str})
    b["key"] = b["PATNO"].str.strip() + "_" + b["EVENT_ID"].str.strip()

    a = a.set_index("key")
    b = b.set_index("key")
    shared = a.index.intersection(b.index)
    print(f"shared keys: {len(shared):,}\n")

    # Participants whose visit set is unchanged between releases. Compare the actual
    # set of EVENT_IDs, not the count: participants exist with four visits in both
    # releases where V10 was replaced by V12, which correctly moves censoring times.
    o_all = pd.read_csv(OLD, sep="\t", usecols=["PATNO", "EVENT_ID"], dtype=str)
    n_all = pd.read_csv(NEW, sep="\t", usecols=["PATNO", "EVENT_ID"], dtype=str)
    o_sets = o_all.groupby("PATNO")["EVENT_ID"].agg(frozenset)
    n_sets = n_all.groupby("PATNO")["EVENT_ID"].agg(frozenset)
    stable_pats = {p for p in o_sets.index if n_sets.get(p) == o_sets[p]}
    A, B = a.loc[shared], b.loc[shared]
    stable = A["PATNO"].isin(stable_pats).to_numpy()
    print(f"participants with an unchanged visit set: {len(stable_pats):,} of "
          f"{A['PATNO'].nunique():,}  ({stable.sum():,} of {len(shared):,} shared rows)\n")

    # ---------------------------------------------------------------- clinical
    print("=" * 78)
    print("A. CLINICAL source columns, on shared keys")
    print("=" * 78)
    rows = []
    for c in clinical:
        bn, ag, vd, pd_ = compare(A[c], B[c])
        rows.append(dict(column=c, both_null=bn, agree=ag, value_differs=vd,
                         presence_differs=pd_))
    cl = pd.DataFrame(rows)
    tot_cmp = cl.agree.sum() + cl.value_differs.sum()
    print(f"  columns: {len(cl)}")
    print(f"  comparable cells: {tot_cmp:,}   agree: {cl.agree.sum():,} "
          f"({cl.agree.sum() / tot_cmp * 100:.4f}%)   differ: {cl.value_differs.sum():,}")
    print(f"  presence mismatches: {cl.presence_differs.sum():,}")
    bad = cl[(cl.value_differs > 0) | (cl.presence_differs > 0)].sort_values(
        "value_differs", ascending=False)
    print(f"  columns with any disagreement: {len(bad)}")
    if len(bad):
        print(bad.head(20).to_string(index=False))

    # ------------------------------------------------------------- harmonized
    # Reported as its own layer rather than folded into the clinical rate: it is
    # proteomic, it is ~50x the size of the clinical layer, and a pooled percentage
    # would be dominated by it while being labelled clinical.
    if harmonized:
        print()
        print("=" * 78)
        print("A2. HARMONIZED proteomic columns, on shared keys")
        print("=" * 78)
        rows = []
        for c in harmonized:
            bn, ag, vd, pd_ = compare(A[c], B[c])
            rows.append(dict(column=c, both_null=bn, agree=ag, value_differs=vd,
                             presence_differs=pd_))
        hm = pd.DataFrame(rows)
        tot_h = hm.agree.sum() + hm.value_differs.sum()
        print(f"  columns: {len(hm)}")
        print(f"  comparable cells: {tot_h:,}   agree: {hm.agree.sum():,} "
              f"({hm.agree.sum() / tot_h * 100:.4f}%)   differ: {hm.value_differs.sum():,}")
        print(f"  presence mismatches: {hm.presence_differs.sum():,}")
        badh = hm[(hm.value_differs > 0) | (hm.presence_differs > 0)]
        print(f"  columns with any disagreement: {len(badh)}")
        if len(badh):
            print(badh.sort_values("value_differs", ascending=False)
                  .head(20).to_string(index=False))
    else:
        print("\nA2. HARMONIZED — baseline carries no harmonized_* columns; not compared.")

    # ---------------------------------------------------------------- derived
    print()
    print("=" * 78)
    print("B. DERIVED columns, participants with an UNCHANGED visit set")
    print("=" * 78)
    rows = []
    for c in derived:
        bn, ag, vd, pd_ = compare(A.loc[stable, c], B.loc[stable, c])
        n_cmp = ag + vd
        rows.append(dict(column=c, comparable=n_cmp, agree=ag, differs=vd,
                         presence_differs=pd_,
                         pct=round(ag / n_cmp * 100, 3) if n_cmp else np.nan))
    dv = pd.DataFrame(rows)

    expected = dv[dv.column.isin(EXPECTED_DIFFERENT) | dv.column.str.startswith("grp_NMC_")]
    other = dv[~dv.index.isin(expected.index)]

    clean = other[(other.differs == 0) & (other.presence_differs == 0)]
    dirty = other[(other.differs > 0) | (other.presence_differs > 0)]
    print(f"  columns expected to reproduce: {len(other)}")
    print(f"    reproduce exactly: {len(clean)}")
    print(f"    disagree:          {len(dirty)}")
    if len(dirty):
        print()
        print(dirty.sort_values("differs", ascending=False).to_string(index=False))

    print()
    print("  deliberately changed (reported, not a failure):")
    print(expected.to_string(index=False))

    print()
    print("=" * 78)
    print("C. DERIVED columns, participants who GAINED visits")
    print("=" * 78)
    rows = []
    for c in derived:
        if c in EXPECTED_DIFFERENT or c.startswith("grp_NMC_"):
            continue
        bn, ag, vd, pd_ = compare(A.loc[~stable, c], B.loc[~stable, c])
        n_cmp = ag + vd
        rows.append(dict(column=c, visit_dependent=c.startswith(VISIT_DEPENDENT),
                         comparable=n_cmp, agree=ag, differs=vd,
                         pct=round(ag / n_cmp * 100, 3) if n_cmp else np.nan))
    gv = pd.DataFrame(rows)
    vd_cols = gv[gv.visit_dependent]
    nv_cols = gv[~gv.visit_dependent]
    print(f"  visit-dependent (slope_/tte_/event_): {len(vd_cols)} columns, "
          f"{vd_cols.differs.sum():,} differing cells — expected, these are functions of "
          "the whole visit history")
    print(f"  not visit-dependent (grp_/CI_/...):   {len(nv_cols)} columns, "
          f"{nv_cols.differs.sum():,} differing cells")
    if nv_cols.differs.sum():
        print(nv_cols[nv_cols.differs > 0].sort_values("differs", ascending=False).to_string(index=False))

    # ------------------------------------------------------------ attribution
    # A derived column that disagrees is only a defect if its clinical INPUTS agree.
    # Where the inputs changed between the two cuts, a changed output is correct.
    print()
    print("=" * 78)
    print("D. Attribution — are derived disagreements explained by clinical input changes?")
    print("=" * 78)

    SOURCES: dict[str, list[str]] = {
        "slope_moca": ["moca"], "slope_updrs3_off": ["updrs3_score"],
        "CI_MOCA": ["moca"], "CI_PI": ["cogstate"],
        "tte_moca_lt26_years": ["moca"], "event_moca_lt26": ["moca"],
        "tte_cogstate_worsen_years": ["cogstate"], "event_cogstate_worsen": ["cogstate"],
        "tte_pm_any_years": ["pm_any"], "event_pm_any": ["pm_any"],
        "tte_pm_cog_any_years": ["pm_cog_any"], "event_pm_cog_any": ["pm_cog_any"],
        "tte_pm_mc_any_years": ["pm_mc_any"], "event_pm_mc_any": ["pm_mc_any"],
        "tte_stage_d_years": ["Stage_D"], "event_stage_d": ["Stage_D"],
        "tte_nsd_2a_to_later_years": ["NSD_STAGE"], "event_nsd_2a_to_later": ["NSD_STAGE"],
        "tte_nsd_2b_to_later_years": ["NSD_STAGE"], "event_nsd_2b_to_later": ["NSD_STAGE"],
        "disease_duration_years": ["age_at_visit", "agediag"],
        "time_to_LEDD_years": ["LEDD", "visit_date"],
    }
    for c in derived:
        if c.startswith("grp_ABC_"):
            SOURCES[c] = ["analytic_subgroup"]
        elif c.startswith("grp_"):
            SOURCES[c] = ["NSD_Status", "subgroup", "NSD_STAGE", "COHORT"]

    Astab, Bstab = A.loc[stable], B.loc[stable]
    # participants with any clinical-input disagreement, per source column
    src_bad: dict[str, set] = {}
    for s in {x for v in SOURCES.values() for x in v} & set(clinical):
        _, _, vd, pdif = compare(Astab[s], Bstab[s])
        if vd or pdif:
            na, nb = Astab[s].isna(), Bstab[s].isna()
            both = ~na & ~nb
            neq = pd.Series(False, index=Astab.index)
            if both.sum():
                x, y = Astab.loc[both, s], Bstab.loc[both, s]
                if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
                    e = ~np.isclose(x.astype(float), y.astype(float), rtol=0, atol=1e-9)
                else:
                    e = x.astype(str).str.strip().to_numpy() != y.astype(str).str.strip().to_numpy()
                neq.loc[both] = e
            neq |= (na ^ nb)
            src_bad[s] = set(Astab.loc[neq, "PATNO"])

    print(f"{'column':30s}{'bad PATNOs':>12s}{'explained':>11s}{'unexplained':>13s}")
    unexplained_total = 0
    for _, r in dirty.iterrows():
        c = r["column"]
        na, nb = Astab[c].isna(), Bstab[c].isna()
        both = ~na & ~nb
        neq = pd.Series(False, index=Astab.index)
        if both.sum():
            x, y = Astab.loc[both, c], Bstab.loc[both, c]
            if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
                e = ~np.isclose(x.astype(float), y.astype(float), rtol=0, atol=1e-9)
            else:
                e = x.astype(str).str.strip().to_numpy() != y.astype(str).str.strip().to_numpy()
            neq.loc[both] = e
        neq |= (na ^ nb)
        pats = set(Astab.loc[neq, "PATNO"])
        explained_by = set().union(*[src_bad.get(s, set()) for s in SOURCES.get(c, [])]) \
            if SOURCES.get(c) else set()
        exp = len(pats & explained_by)
        unexp = len(pats - explained_by)
        unexplained_total += unexp
        print(f"{c:30s}{len(pats):12d}{exp:11d}{unexp:13d}")

    print()
    if unexplained_total == 0:
        print("PASS — every derived disagreement traces to a participant whose clinical "
              "inputs also changed between the two cuts. The derivation logic reproduces.")
    else:
        print(f"{unexplained_total} participant-level disagreement(s) NOT explained by a "
              "clinical input change — inspect these.")

    out = os.path.join(HERE, "build_intermediates", "validation_derived.tab")
    pd.concat([cl.assign(layer="clinical"),
               dv.assign(layer="derived_stable_visits")]).to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    sys.exit(main())
