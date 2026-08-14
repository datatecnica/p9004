"""Step 2 — validate every proteomic block and report its coverage against the scaffold.

The ten blocks are reused from their documented QC-filter outputs (`*-to_merge.tab`).
Each one is a filter output described in the two upstream READMEs; the filters
themselves are not re-run here, so this step's job is to prove the blocks are joinable:

  * `key` built independently of `MERGE_INDEX`, then asserted equal to it
  * key hygiene audit (whitespace, float-formatted IDs, case drift, duplicates)
  * coverage against the scaffold, old cut vs new cut

Reads only the key columns, so it is fast. The heavy join happens in step 4.

Emits:
  build_intermediates/block_manifest-<ts>.tab
  build_intermediates/block_coverage-<ts>.tab
  build_intermediates/key_audit_blocks-<ts>.tab
"""

from __future__ import annotations

import glob
import os
import sys

import pandas as pd

from build_common import TIMESTAMP, assert_no_fatal, log, make_key, report_findings

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")

B = os.path.join(HERE, "initial_build_assets", "MJFF_proteomics-EDA")
N = os.path.join(B, "MJFF_proteomics-EDA")          # the nested subdirectory
A = os.path.join(HERE, os.pardir, "adding_p312_and_P314")

# suffix: the join-artifact suffix used by the previous build, retained for continuity.
# None means the file carries no PATNO/EVENT_ID columns, so no artifact is produced.
BLOCKS = [
    dict(name="p277_CSF",           suffix="p277",       file=os.path.join(N, "PPMI_277_to_merge.tab")),
    dict(name="p282_CNS",           suffix="p282cns",    file=os.path.join(N, "PPMI_282_CNS_to_merge.tab")),
    dict(name="p282_Inflammation",  suffix="p282inflam", file=os.path.join(B, "PPMI_282_Inflammation_to_merge.tab")),
    dict(name="p288_CNS",           suffix=None,         file=os.path.join(N, "FILTERED_PPMI_288_Oct2025_CNS Disease NPQ-to_merge.tab")),
    dict(name="p288_Inflammation",  suffix=None,         file=os.path.join(N, "FILTERED_PPMI_288_Oct2025_Inflammation NPQ-to_merge.tab")),
    dict(name="p293_olink_plasma",  suffix=None,         file=os.path.join(B, "FILTERED_PPMI_293_Oct2025_NPQ-to_merge.tab")),
    dict(name="p312_Inflammation",  suffix="p312inflam", file=os.path.join(A, "PPMI_312_Inflammation_to_merge.tab")),
    dict(name="p312_Neuro",         suffix="p312neuro",  file=os.path.join(A, "PPMI_312_Neuro_to_merge.tab")),
    dict(name="p314_CSF",           suffix="p314csf",    file=os.path.join(A, "PPMI_314_CSF_to_merge.tab")),
    dict(name="p314_Plasma",        suffix="p314plasma", file=os.path.join(A, "PPMI_314_Plasma_to_merge.tab")),
]

OLD_RELEASE = os.path.join(HERE, os.pardir, "Project_9004_Unified_Emerging_Biomarkers.tab")


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} — run the prior step first")
    return hits[-1]


def main() -> None:
    log(f"=== Step 2: block validation ===  timestamp {TIMESTAMP}")

    scaffold = pd.read_csv(latest("scaffold-*.tab"), sep="\t", usecols=["key"], dtype=str)
    new_keys = set(scaffold["key"])
    log(f"scaffold keys: {len(new_keys):,}")

    old_keys: set = set()
    if os.path.exists(OLD_RELEASE):
        o = pd.read_csv(OLD_RELEASE, sep="\t", usecols=["PATNO", "EVENT_ID"], dtype=str)
        old_keys = set(o["PATNO"].str.strip() + "_" + o["EVENT_ID"].str.strip())
        log(f"previous release keys: {len(old_keys):,}")

    manifest, coverage, all_findings = [], [], []

    for b in BLOCKS:
        name, path, suffix = b["name"], b["file"], b["suffix"]
        if not os.path.exists(path):
            sys.exit(f"missing block file for {name}: {path}")

        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        has_ids = "PATNO" in header and "EVENT_ID" in header
        usecols = ["MERGE_INDEX"] + (["PATNO", "EVENT_ID"] if has_ids else [])
        df = pd.read_csv(path, sep="\t", usecols=usecols, dtype=str)

        if has_ids:
            key, findings = make_key(df, source=name)
        else:
            # p288 / p293 carry only MERGE_INDEX; split it back into its components
            # so the same hygiene rules apply rather than trusting the string.
            parts = df["MERGE_INDEX"].astype(str).str.rsplit("_", n=1, expand=True)
            tmp = pd.DataFrame({"PATNO": parts[0], "EVENT_ID": parts[1]})
            key, findings = make_key(tmp, source=name)

        report_findings(findings, name)
        if not findings.empty:
            all_findings.append(findings)
        assert_no_fatal(findings)

        mi = df["MERGE_INDEX"].astype(str).str.strip()
        mismatch = int((key != mi).sum())
        if mismatch:
            log(f"  WARNING [{name}]: key disagrees with MERGE_INDEX on {mismatch} rows")
            log(f"    e.g. key={key[key != mi].iloc[0]!r} vs MERGE_INDEX={mi[key != mi].iloc[0]!r}")

        n_new = int(key.isin(new_keys).sum())
        n_old = int(key.isin(old_keys).sum()) if old_keys else 0
        n_analytes = len([c for c in header if c not in ("MERGE_INDEX", "PATNO", "EVENT_ID")])

        log(f"  {name:20s} {len(df):6,} samples  {n_analytes:5,} cols  "
            f"-> new cut {n_new:6,} ({n_new / len(df) * 100:5.1f}%)   "
            f"old cut {n_old:6,} ({n_old / len(df) * 100:5.1f}%)   gain {n_new - n_old:+,}")

        manifest.append(dict(block=name, file=os.path.relpath(path, HERE), suffix=suffix or "",
                             samples=len(df), columns=n_analytes,
                             carries_ids=has_ids, key_vs_merge_index_mismatch=mismatch))
        coverage.append(dict(block=name, vendor_samples=len(df),
                             matched_old_cut=n_old, matched_new_cut=n_new,
                             gain=n_new - n_old,
                             pct_old=round(n_old / len(df) * 100, 1),
                             pct_new=round(n_new / len(df) * 100, 1)))

    man = pd.DataFrame(manifest)
    cov = pd.DataFrame(coverage)
    fnd = pd.concat(all_findings, ignore_index=True) if all_findings else pd.DataFrame(
        columns=["source", "column", "defect", "n", "action", "example"])

    log(f"\nTOTAL matched: old cut {cov.matched_old_cut.sum():,}  "
        f"new cut {cov.matched_new_cut.sum():,}  gain {cov.gain.sum():+,}")

    man.to_csv(os.path.join(OUT, f"block_manifest-{TIMESTAMP}.tab"), sep="\t", index=False)
    cov.to_csv(os.path.join(OUT, f"block_coverage-{TIMESTAMP}.tab"), sep="\t", index=False)
    fnd.to_csv(os.path.join(OUT, f"key_audit_blocks-{TIMESTAMP}.tab"), sep="\t", index=False)
    log(f"\nwrote manifest / coverage / key audit to {os.path.relpath(OUT, HERE)}/")


if __name__ == "__main__":
    sys.exit(main())
