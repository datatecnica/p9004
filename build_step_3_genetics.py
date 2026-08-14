"""Step 3 — prepare both genetics blocks.

Two blocks, both participant-level (one row per PATNO, propagated to every visit):

  GP2_*    retained unchanged from GP2_NBA_PPMI_PRS.txt via GP2_PGS-to_merge.tab.
           `GP2_clinical_id` becomes the PATNO join key and `GP2_GP2ID` is dropped,
           matching merge_datasets.py and the previous release's 14 GP2_ columns.

  p9001_*  new, from GP2_R12_PPMI_PRS_PCs.csv. `p9001_` is prepended to each source
           column name verbatim, with no other transformation. EVENT_ID is `SC` on
           every row (genetics is static) so it is dropped rather than joined on.

Emits:
  build_intermediates/genetics_gp2-<ts>.tab
  build_intermediates/genetics_p9001-<ts>.tab
  build_intermediates/key_audit_genetics-<ts>.tab
"""

from __future__ import annotations

import glob
import os
import sys

import pandas as pd

from build_common import TIMESTAMP, assert_no_fatal, log, make_patno, report_findings

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")

GP2_SRC = os.path.join(HERE, "initial_build_assets", "MJFF_proteomics-EDA",
                       "GP2_PGS-to_merge.tab")
P9001_SRC = os.path.join(HERE, "clinical_and_p9001_updates", "GP2_R12_PPMI_PRS_PCs.csv")


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} — run the prior step first")
    return hits[-1]


def main() -> None:
    log(f"=== Step 3: genetics ===  timestamp {TIMESTAMP}")

    scaffold = pd.read_csv(latest("scaffold-*.tab"), sep="\t", usecols=["PATNO"], dtype=str)
    participants = set(scaffold["PATNO"])
    log(f"scaffold participants: {len(participants):,}")

    findings_all = []

    # --- GP2 (retained) ----------------------------------------------------
    gp2 = pd.read_csv(GP2_SRC, sep="\t", dtype={"GP2_clinical_id": str})
    log(f"\nGP2  source {os.path.basename(GP2_SRC)}  {gp2.shape[0]:,} x {gp2.shape[1]}")

    gp2 = gp2.rename(columns={"GP2_clinical_id": "PATNO"})
    dropped = [c for c in ("GP2_GP2ID",) if c in gp2.columns]
    gp2 = gp2.drop(columns=dropped)
    if dropped:
        log(f"  dropped external identifier(s): {dropped}  (not in the previous release)")

    patno, f = make_patno(gp2, source="GP2")
    report_findings(f, "GP2")
    assert_no_fatal(f)
    if not f.empty:
        findings_all.append(f)
    gp2["PATNO"] = patno

    n_match = gp2["PATNO"].isin(participants).sum()
    log(f"  {gp2.shape[1] - 1} GP2_ columns, {len(gp2):,} participants, "
        f"{n_match:,} ({n_match / len(gp2) * 100:.1f}%) in scaffold")
    log(f"  columns: {[c for c in gp2.columns if c != 'PATNO']}")

    # --- p9001 (new) -------------------------------------------------------
    p9 = pd.read_csv(P9001_SRC, dtype={"PATNO": str})
    log(f"\np9001  source {os.path.basename(P9001_SRC)}  {p9.shape[0]:,} x {p9.shape[1]}")

    ev = p9["EVENT_ID"].astype(str).unique()
    log(f"  EVENT_ID values: {list(ev)}  -> genetics is static, dropping EVENT_ID")
    assert len(ev) == 1, "p9001 carries more than one EVENT_ID; the static assumption fails"
    p9 = p9.drop(columns=["EVENT_ID"])

    patno, f = make_patno(p9, source="p9001")
    report_findings(f, "p9001")
    assert_no_fatal(f)
    if not f.empty:
        findings_all.append(f)
    p9["PATNO"] = patno

    # Prefix every source column verbatim; PATNO stays the join key.
    p9 = p9.rename(columns={c: f"p9001_{c}" for c in p9.columns if c != "PATNO"})

    n_match = p9["PATNO"].isin(participants).sum()
    log(f"  {p9.shape[1] - 1} p9001_ columns, {len(p9):,} participants, "
        f"{n_match:,} ({n_match / len(p9) * 100:.1f}%) in scaffold")

    prs = [c for c in p9.columns if "_PRS_PRS" in c]
    pcs = [c for c in p9.columns if "_PRS_PC" in c]
    snps = [c for c in p9.columns if c.startswith("p9001_rs")]
    other = [c for c in p9.columns if c not in prs + pcs + snps and c != "PATNO"]
    log(f"    {len(prs)} PRS: {prs}")
    log(f"    {len(pcs)} PCs, {len(snps)} SNP dosages, other: {other}")

    # Overlap between the two genetics blocks, since both are retained.
    both = set(gp2["PATNO"]) & set(p9["PATNO"])
    log(f"\n  participants in both genetics blocks: {len(both):,}  "
        f"(GP2 only {len(set(gp2['PATNO']) - set(p9['PATNO'])):,}, "
        f"p9001 only {len(set(p9['PATNO']) - set(gp2['PATNO'])):,})")

    # --- write -------------------------------------------------------------
    f_gp2 = os.path.join(OUT, f"genetics_gp2-{TIMESTAMP}.tab")
    f_p9 = os.path.join(OUT, f"genetics_p9001-{TIMESTAMP}.tab")
    f_au = os.path.join(OUT, f"key_audit_genetics-{TIMESTAMP}.tab")

    gp2.to_csv(f_gp2, sep="\t", index=False)
    p9.to_csv(f_p9, sep="\t", index=False)
    fnd = pd.concat(findings_all, ignore_index=True) if findings_all else pd.DataFrame(
        columns=["source", "column", "defect", "n", "action", "example"])
    fnd.to_csv(f_au, sep="\t", index=False)

    log(f"\nwrote {os.path.basename(f_gp2)}  ({gp2.shape[0]:,} x {gp2.shape[1]})")
    log(f"wrote {os.path.basename(f_p9)}  ({p9.shape[0]:,} x {p9.shape[1]})")


if __name__ == "__main__":
    sys.exit(main())
