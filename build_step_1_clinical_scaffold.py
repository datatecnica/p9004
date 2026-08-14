"""Step 1 — build the clinical scaffold from the refreshed curated data cut.

Replaces PPMI_Proteomics_Data_Cut_INTERNAL_20251215 (10,354 rows / 2,025 participants)
with PPMI_Curated_Data_Cut_Public_20260511 sheet `20260511`
(19,450 rows / 4,788 participants).

Emits:
  build_intermediates/scaffold-<ts>.tab            the clinical base, with `key`
  build_intermediates/scaffold_dict-<ts>.tab       the cut's own dictionary, collapsed 1:1
  build_intermediates/key_audit-<ts>.tab           hygiene findings

The row count written here is the invariant every later step asserts against.
"""

from __future__ import annotations

import os
import sys

import pandas as pd

from build_common import TIMESTAMP, assert_no_fatal, log, make_key, report_findings

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "clinical_and_p9001_updates",
                   "PPMI_Curated_Data_Cut_Public_20260511 (2).xlsx")
SHEET = "20260511"
OUT = os.path.join(HERE, "build_intermediates")

EXPECTED_ROWS = 19450
EXPECTED_PARTICIPANTS = 4788


def collapse_dictionary(dd: pd.DataFrame) -> pd.DataFrame:
    """The cut ships its dictionary sparsely: a variable's Code/Decode pairs occupy
    extra rows with every other field blank. Collapse to one row per variable, folding
    the code list into a single cell, so it can be asserted 1:1 against the dataset."""
    dd = dd.copy()
    # A new variable starts wherever `Variable` is populated.
    dd["_grp"] = dd["Variable"].notna().cumsum()
    dd = dd[dd["_grp"] > 0]

    def _codes(g: pd.DataFrame, col: str) -> str:
        vals = [str(v).strip() for v in g[col].dropna() if str(v).strip() not in ("", "-")]
        return "; ".join(dict.fromkeys(vals)) if vals else "-"

    rows = []
    for _, g in dd.groupby("_grp", sort=True):
        head = g.iloc[0]
        rows.append({
            "Category": head["Category"],
            "Variable": head["Variable"],
            "Description": head["Description"],
            "Code": _codes(g, "Code"),
            "Decode": _codes(g, "Decode"),
            "Derived Variable": head["Derived Variable"],
            "Original Variable(s)": head["Original Variable(s)"],
            "Original Dataset(s)": head["Original Dataset(s)"],
            "Derivation Notes": head["Derivation Notes"],
        })
    out = pd.DataFrame(rows)
    out["Variable"] = out["Variable"].astype(str).str.strip()
    return out


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    log(f"=== Step 1: clinical scaffold ===  timestamp {TIMESTAMP}")
    log(f"source: {os.path.basename(SRC)}  sheet {SHEET}")

    df = pd.read_excel(SRC, sheet_name=SHEET)
    log(f"  read {df.shape[0]:,} rows x {df.shape[1]} columns")

    if df.shape[0] != EXPECTED_ROWS:
        log(f"  WARNING: expected {EXPECTED_ROWS:,} rows, got {df.shape[0]:,}. "
            "The cut may have been reissued — check before trusting downstream asserts.")

    # --- key ---------------------------------------------------------------
    key, findings = make_key(df, source="clinical_scaffold")
    report_findings(findings, "clinical_scaffold")
    assert_no_fatal(findings)

    df.insert(0, "key", key)

    # MERGE_INDEX is retained for continuity with the existing *-to_merge.tab files,
    # which are all keyed on it. `key` is built independently; step 4 asserts they agree.
    df.insert(1, "MERGE_INDEX", key)

    n_pat = df["PATNO"].nunique()
    log(f"  participants: {n_pat:,}")
    if n_pat != EXPECTED_PARTICIPANTS:
        log(f"  WARNING: expected {EXPECTED_PARTICIPANTS:,} participants")

    assert df["key"].is_unique, "duplicate key in scaffold — cannot be a merge base"
    log(f"  key is unique across {len(df):,} rows")
    log(f"  EVENT_ID vocabulary: {sorted(df['EVENT_ID'].astype(str).unique())}")

    # --- dictionary --------------------------------------------------------
    dd_raw = pd.read_excel(SRC, sheet_name="Data dictionary")
    dd = collapse_dictionary(dd_raw)
    log(f"  dictionary: {len(dd_raw)} raw rows -> {len(dd)} variables")

    cols = set(df.columns) - {"key", "MERGE_INDEX"}

    # The workbook disagrees with itself on capitalisation for a handful of variables
    # (`age_datscan` on the data sheet vs `age_DATSCAN` in the dictionary, and four
    # others). The data sheet is authoritative — it is what becomes a dataset column —
    # so realign the dictionary onto it rather than the reverse.
    by_lower = {c.lower(): c for c in cols}
    drift = []
    fixed = []
    for v in dd["Variable"]:
        if v not in cols and v.lower() in by_lower:
            drift.append((v, by_lower[v.lower()]))
            fixed.append(by_lower[v.lower()])
        else:
            fixed.append(v)
    dd["Variable"] = fixed
    if drift:
        log(f"  dictionary case drift realigned onto data-sheet spelling: {len(drift)}")
        for a, b in drift:
            log(f"    {a}  ->  {b}")

    documented = set(dd["Variable"])
    undocumented = sorted(cols - documented)
    orphaned = sorted(documented - cols)
    log(f"  scaffold columns without a dictionary entry: {len(undocumented)}")
    if undocumented:
        log(f"    {undocumented[:12]}{' ...' if len(undocumented) > 12 else ''}")
    log(f"  dictionary entries with no scaffold column:  {len(orphaned)}")
    if orphaned:
        log(f"    {orphaned[:12]}{' ...' if len(orphaned) > 12 else ''}")

    assert not undocumented and not orphaned, \
        "scaffold dictionary is not 1:1 with the scaffold columns"
    log("  dictionary is 1:1 with scaffold columns")

    # --- write -------------------------------------------------------------
    f_sc = os.path.join(OUT, f"scaffold-{TIMESTAMP}.tab")
    f_dd = os.path.join(OUT, f"scaffold_dict-{TIMESTAMP}.tab")
    f_au = os.path.join(OUT, f"key_audit-{TIMESTAMP}.tab")

    df.to_csv(f_sc, sep="\t", index=False)
    dd.to_csv(f_dd, sep="\t", index=False)
    findings.to_csv(f_au, sep="\t", index=False)

    log(f"\nwrote {f_sc}  ({df.shape[0]:,} x {df.shape[1]})")
    log(f"wrote {f_dd}  ({dd.shape[0]:,} entries)")
    log(f"wrote {f_au}  ({len(findings)} finding rows)")


if __name__ == "__main__":
    sys.exit(main())
