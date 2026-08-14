"""Step 6 — per-block proteomic principal components.

Twelve proteomic panel blocks x 10 PCs = 120 columns, named
`<prefix>_earliest_visit_PC1..PC10`. PCs are computed from proteomic analytes only —
no clinical, genetic or PLATE_ID column enters the PCA.

`All_Combined` / `All_Combined_no_p293` are NOT built here. The previous versions
required a sample to carry data in every project, which across the six current projects
leaves 100 participants — too thin to be worth a block, and the pooled matrix is not on
a common scale anyway (NPQ vs NPX, plate-control vs intensity normalization, differing
panel versions). A cross-project block belongs downstream of harmonizing projects by
assay and analyte type, where the analytes are z-scored within project first.

Method, per block (unchanged from the p312/p314 build):
  1. take each participant's EARLIEST VISIT THAT CARRIES DATA for that block
  2. drop assays missing in >=20% of those samples
  3. mean-impute the residual gaps
  4. StandardScaler, then PCA(n_components=10)
  5. join scores back on PATNO, so a participant's baseline PCs carry to all visits

Step 1 takes the earliest qualifying ROW, not `GroupBy.first()`. The legacy pipeline
used `.first()`, which takes the first non-null value *per column* — so a participant's
"baseline" row was a composite stitched from several visits, while EVENT_ID (never null)
still read `BL`. That defect is documented in the p312/p314 README and is not inherited
here.

PC values are not comparable to the previous release even for untouched blocks: every
block's PCA sample set grew with the refreshed scaffold, so expect sign flips and
component reordering.

Emits:
  build_intermediates/block_pcs-<ts>.tab
  build_intermediates/pc_variance-<ts>.tab
  EDA_PCA_plots/EDA_PCA-<block>-{eigenvectors,loadings}-<ts>.tab
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

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "build_intermediates")
PLOTS = os.path.join(HERE, "EDA_PCA_plots")

N_COMPONENTS = 10
MAX_MISSING = 0.20          # drop assays missing in >= this fraction of block samples

PANEL_PREFIXES = [
    "p277_CSF", "p282_CNS_CSF", "p282_Inflammation_CSF",
    "p288_CNS_plasma", "p288_Inflammation_plasma", "p293_olink_plasma",
    "p312_Inflammation_CSF", "p312_Inflammation_Plasma",
    "p312_Neuro_CSF", "p312_Neuro_Plasma",
    "p314_CSF", "p314_Plasma",
]

# Visit ordering for the "earliest visit" rule. YEAR is the primary key; EVENT_ID
# breaks ties within a year deterministically.
VISIT_ORDER = ["BL", "PW", "ST", "V04", "V06", "V08", "V10", "V12", "V13", "V14",
               "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22"]


def latest(pattern: str) -> str:
    hits = sorted(glob.glob(os.path.join(OUT, pattern)))
    if not hits:
        sys.exit(f"missing {pattern} — run the prior step first")
    return hits[-1]


def run_block(name: str, feats: list[str], meta: pd.DataFrame,
              mat: pd.DataFrame, variance: list) -> pd.DataFrame | None:
    """PCA one block; return a PATNO-indexed frame of PC scores."""
    sub = mat[feats]
    has_data = sub.notna().any(axis=1)
    if not has_data.any():
        log(f"  {name:26s} SKIPPED — no rows carry data")
        return None

    # earliest qualifying ROW per participant
    idx = meta.loc[has_data].sort_values(["PATNO", "_order"]).groupby("PATNO", sort=False).head(1).index
    X = sub.loc[idx]
    who = meta.loc[idx, "PATNO"]

    # drop sparse assays, then mean-impute what is left
    miss = X.isna().mean()
    keep = miss[miss < MAX_MISSING].index.tolist()
    if len(keep) < 2:
        log(f"  {name:26s} SKIPPED — only {len(keep)} feature(s) survive the <20% missing cut")
        return None
    X = X[keep]
    n_imputed_cells = int(X.isna().sum().sum())
    pct_cells = n_imputed_cells / (X.shape[0] * X.shape[1]) * 100
    X = X.fillna(X.mean())
    # a feature that is all-NaN across these samples leaves NaN after mean-impute
    X = X.loc[:, X.notna().all()]

    Z = StandardScaler().fit_transform(X.to_numpy(dtype=np.float64))
    n_comp = min(N_COMPONENTS, *Z.shape)
    pca = PCA(n_components=n_comp, random_state=0)
    scores = pca.fit_transform(Z)

    cols = [f"{name}_earliest_visit_PC{i + 1}" for i in range(n_comp)]
    sc = pd.DataFrame(scores, columns=cols)
    sc.insert(0, "PATNO", who.to_numpy())
    sc = sc.drop_duplicates(subset=["PATNO"]).set_index("PATNO")

    evr = pca.explained_variance_ratio_
    at_bl = float((meta.loc[idx, "EVENT_ID"] == "BL").mean() * 100)
    log(f"  {name:26s} {X.shape[0]:5,} samples x {X.shape[1]:5,} feats  "
        f"PC1 {evr[0] * 100:5.1f}%  cum {evr.sum() * 100:5.1f}%  "
        f"imputed {pct_cells:4.2f}%  at BL {at_bl:5.1f}%")

    variance.append(dict(block=name, samples=X.shape[0], features=X.shape[1],
                         features_dropped=len(feats) - X.shape[1],
                         pct_cells_imputed=round(pct_cells, 3),
                         pct_samples_at_BL=round(at_bl, 1),
                         pc1_var=round(float(evr[0]) * 100, 2),
                         cum_var_pc1_10=round(float(evr.sum()) * 100, 2)))

    # eigenvectors (participant scores) and loadings, per the existing convention
    os.makedirs(PLOTS, exist_ok=True)
    eig = sc.reset_index()
    eig.to_csv(os.path.join(PLOTS, f"EDA_PCA-{name}-eigenvectors-{TIMESTAMP}.tab"),
               sep="\t", index=False)
    pd.DataFrame(pca.components_.T, index=X.columns, columns=cols).reset_index(
        names="feature").to_csv(
        os.path.join(PLOTS, f"EDA_PCA-{name}-loadings-{TIMESTAMP}.tab"), sep="\t", index=False)
    return sc


def main() -> None:
    log(f"=== Step 6: block PCs ===  timestamp {TIMESTAMP}")
    src = latest("merged-*.tab")

    header = pd.read_csv(src, sep="\t", nrows=0).columns.tolist()
    block_feats = {p: [c for c in header
                       if c.startswith(p + "_") and "PLATE_ID" not in c]
                   for p in PANEL_PREFIXES}
    all_feats = sorted({c for f in block_feats.values() for c in f})
    log(f"analyte columns across {len(PANEL_PREFIXES)} panels: {len(all_feats):,}")

    meta_cols = ["key", "PATNO", "EVENT_ID", "YEAR"]
    log("reading analyte matrix (float32)...")
    mat = pd.read_csv(src, sep="\t", usecols=all_feats,
                      dtype={c: np.float32 for c in all_feats})
    meta = pd.read_csv(src, sep="\t", usecols=meta_cols,
                       dtype={"key": str, "PATNO": str, "EVENT_ID": str})
    log(f"  matrix {mat.shape[0]:,} x {mat.shape[1]:,}  "
        f"({mat.memory_usage(deep=True).sum() / 1e9:.2f} GB)")

    order = {v: i for i, v in enumerate(VISIT_ORDER)}
    meta["_order"] = (pd.to_numeric(meta["YEAR"], errors="coerce").fillna(9e3) * 100
                      + meta["EVENT_ID"].map(order).fillna(99))

    variance: list = []
    scores: list[pd.DataFrame] = []

    log("\n--- panel blocks ---")
    for p in PANEL_PREFIXES:
        s = run_block(p, block_feats[p], meta, mat, variance)
        if s is not None:
            scores.append(s)

    del mat

    pcs = pd.concat(scores, axis=1)
    out = pd.DataFrame({"key": meta["key"]})
    joined = meta["PATNO"].map(lambda p: p).to_frame("PATNO").join(pcs, on="PATNO")
    out = pd.concat([out, joined.drop(columns=["PATNO"])], axis=1)

    assert len(out) == len(meta), "PC join changed the row count"
    log(f"\nPC columns: {out.shape[1] - 1} across {len(variance)} blocks")

    f_out = os.path.join(OUT, f"block_pcs-{TIMESTAMP}.tab")
    out.to_csv(f_out, sep="\t", index=False)
    pd.DataFrame(variance).to_csv(
        os.path.join(OUT, f"pc_variance-{TIMESTAMP}.tab"), sep="\t", index=False)
    log(f"wrote {os.path.basename(f_out)}  ({out.shape[0]:,} x {out.shape[1]})")


if __name__ == "__main__":
    sys.exit(main())
