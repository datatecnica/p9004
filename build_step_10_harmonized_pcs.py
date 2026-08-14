"""Step 10 — principal components on the six harmonized blocks.

Six blocks x 10 PCs = 60 columns, named `harmonized_<block>_earliest_visit_PC1..PC10`.
Method is identical to the per-project PC blocks in step 6:

  1. each participant's EARLIEST VISIT carrying data for that block
  2. drop analytes missing in >=20% of those samples
  3. mean-impute the residual gaps
  4. StandardScaler, then PCA(n_components=10)
  5. join scores back on PATNO, so baseline PCs carry to all of a participant's visits

These must be computed fresh rather than inherited from the per-project blocks: a
participant's earliest measurement for a harmonized block can now come from a different
project than before, which moves the origin for slopes, time-to-event and
disease-duration covariates.

Step 1 takes the earliest qualifying ROW, not `GroupBy.first()` — the legacy pipeline's
`.first()` took the first non-null value *per column*, stitching a composite "baseline"
from several visits while EVENT_ID still read BL. That defect is not inherited here.

Proteomic analytes only: no clinical, genetic, provenance or PLATE_ID column enters
the PCA.

Emits:
  build_intermediates/harmonized_pcs-<ts>.tab
  build_intermediates/harmonized_pc_variance-<ts>.tab
  EDA_PCA_plots/EDA_PCA-harmonized_<block>-{eigenvectors,loadings}-<ts>.tab
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from build_common import TIMESTAMP, log
from build_step_6_block_pcs import MAX_MISSING, N_COMPONENTS, VISIT_ORDER, run_block

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")

BLOCKS = ["olink_plasma", "olink_csf", "nulisa_cns_plasma", "nulisa_cns_csf",
          "nulisa_inf_plasma", "nulisa_inf_csf"]


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} — run the prior step first")
    return hits[-1]


def main() -> None:
    log(f"=== Step 10: harmonized block PCs ===  timestamp {TIMESTAMP}")
    f_harm = latest("harmonized-*.tab")
    header = pd.read_csv(f_harm, sep="\t", nrows=0).columns.tolist()

    block_feats = {b: [c for c in header
                       if c.startswith(f"harmonized_{b}_")
                       and not c.endswith("_src")
                       and "earliest_visit_PC" not in c]
                   for b in BLOCKS}
    feats = sorted({c for v in block_feats.values() for c in v})
    log(f"harmonized analyte columns: {len(feats):,}")

    log("reading harmonized matrix (float32) ...")
    mat = pd.read_csv(f_harm, sep="\t", usecols=feats,
                      dtype={c: np.float32 for c in feats})
    meta = pd.read_csv(latest("merged-*.tab"), sep="\t",
                       usecols=["key", "PATNO", "EVENT_ID", "YEAR"],
                       dtype={"key": str, "PATNO": str, "EVENT_ID": str})
    log(f"  matrix {mat.shape[0]:,} x {mat.shape[1]:,} "
        f"({mat.memory_usage(deep=True).sum() / 1e9:.2f} GB)")

    order = {v: i for i, v in enumerate(VISIT_ORDER)}
    meta["_order"] = (pd.to_numeric(meta["YEAR"], errors="coerce").fillna(9e3) * 100
                      + meta["EVENT_ID"].map(order).fillna(99))

    variance: list = []
    scores: list[pd.DataFrame] = []
    log("")
    for b in BLOCKS:
        s = run_block(f"harmonized_{b}", block_feats[b], meta, mat, variance)
        if s is not None:
            scores.append(s)
    del mat

    pcs = pd.concat(scores, axis=1)
    joined = meta[["PATNO"]].join(pcs, on="PATNO").drop(columns=["PATNO"])
    out = pd.concat([meta[["key"]].reset_index(drop=True),
                     joined.reset_index(drop=True)], axis=1)

    assert len(out) == len(meta), "PC join changed the row count"
    log(f"\nPC columns: {out.shape[1] - 1} across {len(variance)} blocks")

    f_out = os.path.join(OUT, f"harmonized_pcs-{TIMESTAMP}.tab")
    out.to_csv(f_out, sep="\t", index=False)
    pd.DataFrame(variance).to_csv(
        os.path.join(OUT, f"harmonized_pc_variance-{TIMESTAMP}.tab"), sep="\t", index=False)
    log(f"wrote {os.path.basename(f_out)}  ({out.shape[0]:,} x {out.shape[1]})")


if __name__ == "__main__":
    sys.exit(main())
