# Phase 1 wrap — rebuild plan for Project 9004

Rebuild `Project_9004_Unified_Emerging_Biomarkers.tab` from scratch on a refreshed
clinical scaffold, following the filters established in
the MJFF proteomics EDA README and the p312/p314 build README. Both are working
documents in the restricted data directory and are not part of this repository.

Status: **BUILT.** This is the planning and source-audit document; for what was
actually produced and how it validates, see [README.md](README.md) and
[comparison_to_previous_release.md](comparison_to_previous_release.md).

Deltas from this plan, decided during the build:
- **All_Combined PC blocks dropped** — 12 blocks / 120 PC columns, not 14 / 140. The
  six-project intersection is 100 participants, and a cross-project block belongs
  downstream of harmonizing by assay and analyte type.
- **54 derived columns**, not the ~66 estimated below.
- **`lowput_ratio` was a pre-existing bug**, not just a series change — see README.
- **p277 and p282_CNS were not blocked**; their filtered files were in the nested
  `MJFF_proteomics-EDA/MJFF_proteomics-EDA/` subdirectory.

---

## Headline: the scaffold refresh is the whole point

The current release is built on `PPMI_Proteomics_Data_Cut_INTERNAL_20251215`
(10,354 rows, 2,025 participants). The replacement,
`PPMI_Curated_Data_Cut_Public_20260511 (2).xlsx` sheet `20260511`, carries
**19,450 rows across 4,788 participants** — 2.4× the participants.

The p312/p314 README named this exactly: *"Refreshing the clinical cut is the single
change that would most improve these projects' usable N."* Measured against every
QC-passing vendor file, joining on `PATNO_EVENT_ID`:

| Vendor block | Samples | Matched to **old** cut | Matched to **new** cut | Gain |
|---|---|---|---|---|
| p277 CSF | 1,200 | 1,167 (97.2%) | 1,167 (97.2%) | +0 |
| p282 CNS | 2,301 | 2,282 (99.2%) | 2,282 (99.2%) | +0 |
| p282 Inflammation | 1,915 | 1,899 (99.2%) | 1,899 (99.2%) | +0 |
| p288 CNS | 2,754 | 2,665 (96.8%) | 2,665 (96.8%) | +0 |
| p288 Inflammation | 2,758 | 2,669 (96.8%) | 2,669 (96.8%) | +0 |
| p293 plasma | 850 | 784 (92.2%) | 784 (92.2%) | +0 |
| **p312 Inflammation** | 2,204 | 478 (21.7%) | **1,995 (90.5%)** | **+1,517** |
| **p312 Neuro** | 2,187 | 471 (21.5%) | **1,978 (90.4%)** | **+1,507** |
| **p314 CSF** | 2,128 | 1,236 (58.1%) | **2,098 (98.6%)** | **+862** |
| **p314 Plasma** | 3,399 | 1,561 (45.9%) | **2,942 (86.6%)** | **+1,381** |
| **TOTAL** | | 15,212 | **20,479** | **+5,267** |

**The entire gain is p312 and p314.** The six older blocks were already matched at
92–99% and do not move — expected, since those projects were assayed on the cohort the
INTERNAL proteomics cut was built around, whereas p312/p314 reach a much wider
participant set.

Project 312 goes from a 21.7% salvage to essentially complete. That converts it from
"effectively cross-sectional, 382 usable participants" into a real block, and it is the
single result that justifies the rebuild.

The residual non-matches are stable across both cuts and are the `V01/V02/V09/V11`
visit codes that neither clinical cut carries for anyone, plus a small number of
participants absent from both.

**The refresh is non-destructive at the participant level.** Every one of the 2,025
old PATNOs is present in the new cut; only 5 of 10,354 old `PATNO_EVENT_ID` keys are
absent. The `EVENT_ID` vocabulary is identical (`BL, PW, ST, V04, V06, V08, V10,
V12–V22`), so no visit remapping is needed and the residual p312/p314 losses are still
the `V01/V02/V09/V11` codes that neither cut carries.

---

## Source availability

**All ten blocks have usable inputs.** Note the nested
`initial_build_assets/MJFF_proteomics-EDA/MJFF_proteomics-EDA/` subdirectory — it holds
the p277, p282_CNS and both p288 filtered files plus the p282 CNS raw xlsx, and is easy
to miss because the parent directory has near-identical contents.

| Block | Raw vendor source | Filtered `*-to_merge` | Rebuild path |
|---|---|---|---|
| Clinical scaffold | ✅ curated xlsx (new) | — | from raw |
| p277 Olink CSF | ❌ (methods PDFs only) | ✅ `PPMI_277_to_merge.tab` † | reuse filtered |
| p282 CNS | ✅ `…CNSDiseasePanel_NPQCounts_20260120.xlsx` † | ✅ `PPMI_282_CNS_to_merge.tab` † | either |
| p282 Inflammation | ❌ (report docx only) | ✅ `PPMI_282_Inflammation_to_merge.tab` | reuse filtered |
| p288 CNS + Inflammation | ✅ `PPMI-288-P-000671-NULISA-20250828.xlsx` (both panels as sheets) | ✅ both † | either |
| p293 Olink plasma | ⚠️ stale `…20250328.parquet` | ✅ verified current | **reuse filtered** |
| p312 Inflammation + Neuro | ✅ 4 xlsx | ✅ both | either |
| p314 CSF + Plasma | ✅ 2 parquet | ✅ both | either |
| GP2 (existing block) | ✅ `GP2_NBA_PPMI_PRS.txt` | ✅ | either |
| p9001 (new block) | ✅ `GP2_R12_PPMI_PRS_PCs.csv` | — | from raw |

† in the nested `MJFF_proteomics-EDA/MJFF_proteomics-EDA/` subdirectory.

**Default: reuse the filtered `*-to_merge.tab` files.** Every one is a documented
QC-filter output keyed on `MERGE_INDEX`, none contains duplicate keys, and reusing them
makes the proteomic layer identical-by-construction to the current release — which
turns validation check 4 (value identity on the unchanged path) into a real test rather
than a tautology. Re-deriving from raw where possible is a useful independent check but
is not required for the rebuild, and for p293 it is actively wrong.

### On p293 — reuse the filtered file, do not re-derive

`data_from_sources/` holds the **20250328** parquet, but the 2026-01-02 refresh moved
p293 to **20251121**. Rebuilding p293 from the parquet on disk would silently regress
it by one data release.

The filtered file is the current one. Verified directly: 60 sampled analytes in
`FILTERED_PPMI_293_Oct2025_NPQ-to_merge.tab` are **value-identical** (`atol=1e-9`) to
the shipped release across all matched rows.

---

## Blocks

Ten data blocks plus a derived layer. Current release column counts, for reference:

| # | Block | Platform / metric | Fluid | Cols in current release |
|---|---|---|---|---|
| 0 | Clinical scaffold | — | — | 254 source columns |
| 1 | `p277_CSF` | Olink Explore HT / NPX | CSF | 5,426 |
| 2 | `p282_CNS_CSF` | NULISA / NPQ | CSF | 139 |
| 3 | `p282_Inflammation_CSF` | NULISA / NPQ | CSF | 257 |
| 4 | `p288_CNS_plasma` | NULISA / NPQ | plasma | 142 |
| 5 | `p288_Inflammation_plasma` | NULISA / NPQ | plasma | 260 |
| 6 | `p293_olink_plasma` | Olink Explore HT / PCNormalizedNPX | plasma | 5,426 |
| 7 | `p312_{Inflammation,Neuro}_{CSF,Plasma}` | NULISA / NPQ | both | 975 |
| 8 | `p314_{CSF,Plasma}` | Olink Explore HT / NPX | both | 10,839 |
| 9 | `GP2_*` | genetics (existing) | — | 16 |
| 10 | `p9001_*` | genetics (**new**) | — | — |

> **Scales are not comparable across blocks.** NPQ vs NPX, plate-control vs intensity
> normalization, differing panel versions, no formal bridging. Any cross-block
> analysis must z-score within block first. This carries over unchanged; see the
> p312/p314 README's warning section.

---

## Filtering spec — carried over verbatim

No filter changes. Reproduced here so the build scripts can be checked against one page.

**p288 (NULISA)** — `SampleQC == "PASS"`, `HEMOLYSIS == "Passed"`; both panels.

**p293 (Olink)** — `AssayQC == "PASS"`, `SampleQC == "PASS"`, and
`SampleBlockQCWarn / SampleBlockQCFail / BlockQCFail / AssayQCWarn == 1`
(1 = baseline "no issues"; >1 counts warnings/failures). Uses `PCNormalizedNPX`.

**p282 (NULISA)** — `SampleQC == "passed"` (lowercase vocabulary, unlike p312),
`Biofluid == "CSF"`, non-missing `PATNO` / `CLINICAL_EVENT`.

**p277 (Olink)** — non-missing `PATNO` / `EVENT_ID`, `AssayQC == "PASS"`,
`SampleQC == "PASS"`.

**p312 (NULISA)** — header read with `skiprows=1` (plate-specific LOD note on row 1);
`SampleQC == "PASS"`, `SampleType == "Sample"` (drops NC/IPC/SC wells),
`Biofluid == "CSF"|"Plasma"` per file, non-missing `PatNo` / `ClinicalEvent`.
Asymmetric CSF/plasma column sets are expected — 10 high-abundance plasma proteins
report no NPQ in CSF.

**p314 (Olink)** — nine filters per the methods document: `AssayType == "assay"`;
`SampleType == "SAMPLE"`; non-missing `PATNO`/`EVENT_ID`; `SampleQC == "PASS"`;
`AssayQC == "PASS"` and `AssayQCWarn != 2`; `Normalization != "EXCLUDED"`;
`SampleBlockQCWarn/SampleBlockQCFail/BlockQCFail <= 1`; `LAB_QC_FLAG` carries no code
6 or 7 (code 2 = below LOD is **retained** by design); NPX beyond ±5 MAD-scaled SD of
the per-`OlinkID` median masked to `NaN`. Uses `NPX`, **not** `PCNormalizedNPX`.

---

## What changes, beyond the scaffold

### 1. `lowput_ratio` must be redefined — MIA is not a drop-in

The curated cut drops all 15 `DATSCAN_*` / `con_*` / `ips_*` / `mean_*` columns and
keeps only the `MIA_*` DaTscan series. `build_dat_lowput_ratio()` currently computes
`min(DATSCAN_PUTAMEN_L, DATSCAN_PUTAMEN_R) / lowput_expected`.

Tested on the current release, which carries **both** series:

| Pair | Non-null (old / MIA) | Identical? | Pearson r |
|---|---|---|---|
| `DATSCAN_PUTAMEN_L` vs `MIA_PUTAMEN_L` | 4,312 / 5,178 | no (max Δ 2.36) | 0.896 |
| `DATSCAN_PUTAMEN_R` vs `MIA_PUTAMEN_R` | 4,312 / 5,178 | no (max Δ 1.87) | 0.886 |
| `lowput_expected` vs `MIA_LOWPUT_EXPECTED` | 4,312 / 5,178 | no (max Δ 0.78) | 0.905 |

These are **different quantification pipelines, not a rename.** The MIA series has
~20% more coverage, so the switch is a net gain — but `lowput_ratio` and
`slope_lowput_ratio` will not reproduce their old values, and the 0.75 DaT-positivity
convention should be re-checked against the MIA distribution rather than assumed.
This must be recorded in the dictionary's `Derivation Notes`.

### 2. 59 source-clinical columns drop, 17 arrive

Dropped, in three coherent groups:

- **12 assay-availability / near-visit indicators** (`NULISA_CSF`, `NULISA_CSF_VIS`,
  `NULISA_CSF_VISCM`, … `Olink_CSF_VISCM`). These powered the `'1*'` near-visit rescue
  mechanism, which the p312/p314 audit reported on but deliberately never applied. No
  derived variable depends on them. Loss is cosmetic.
- **32 SAA replicate/detail columns** (`Fmax_*`, `TTT_*`, `AUC_*`, `TSmax_*`,
  `SLOPEMax_*`, `T50_*`, `SLOPE_*`, `CSFSAARundate`, `CSFSAATRANCHE`, `CSFSAA_150h`).
  Headline `CSFSAA` and `CSFSAA_assay` are retained.
- **15 DaTscan columns** — see §1.

Arriving: `Clinical_Stage`, `Death_Status`, `age_death`, `Death_Date`, `CGI`, `PGI`,
`PDAQ27`, `nqol_{cogns,comms,lefs,uefs}`, `NFL_CSF_ULOD`, `Roche_IL_1B_CSF`,
`Roche_IL_6_CSF(_LLOD)`, `Roche_YKL40_CSF(_LLOD)`.

**No group-defining variable is lost.** Every derived flag keys off `subgroup`,
`analytic_subgroup`, `NSD_Status`, `NSD_STAGE`, `COHORT` — all present.

### 3. New `p9001_*` genetics block, alongside the retained `GP2_*` block

Both blocks are kept, per instruction. `GP2_*` rebuilds unchanged from
`GP2_NBA_PPMI_PRS.txt` (3,359 participants, 16 cols).

`GP2_R12_PPMI_PRS_PCs.csv` adds **4,985 participants** (+48% over GP2) × 172 data
columns: 4 PRS variants (`PRS157/154/152/149`), `InfPop` ancestry label, `PC1–PC10`,
and **157 individual SNP dosages** (PD GWAS risk loci, 0/1/2 coded).

`EVENT_ID` is `SC` for every row and `PATNO` is unique — genetics is static, so merge
**on `PATNO` only** and drop `EVENT_ID`, matching how `GP2_*` was joined.

Naming: `p9001_` is prepended to the source column name verbatim, with no other
transformation, and no other block's names change — so
`p9001_Genetic_PRS_PRS157`, `p9001_Genetic_PRS_InfPop`, `p9001_Genetic_PRS_PC1..PC10`,
`p9001_rs356182_A`.

Ancestry distribution: EUR 3,727 · AJ 1,008 · AMR 70 · AFR 39 · AAC 36 · EAS 35 ·
MDE 23 · CAH 21 · CAS 13 · SAS 9 · FIN 4.

> **Downstream note (phase 2):** `batch.yaml` strata are defined on `GP2_nba_label`.
> `p9001_Genetic_PRS_InfPop` is the larger and more current label but is *not*
> automatically adopted — switching strata is an analysis decision, deliberately out of
> scope here. Worth noting the two labels will disagree for some participants; a
> crosstab belongs in the step 8 comparison.

### 4. `key` — new hard-keyed merge column with a hygiene audit

Per the notes, a `key = PATNO + "_" + VISIT` string column becomes the merge key, with
an explicit audit that **fails loudly** rather than coercing. Checks on every input
before any join:

- leading/trailing whitespace, non-breaking spaces (` `), zero-width characters
- `PATNO` arriving as float (`"3000.0"`) or with Excel scientific notation
- `EVENT_ID` case/spacing variants (`bl` / `BL ` / `V04`)
- empty, `NaN`, or literal `"nan"` / `"None"` components
- duplicate keys within any single block

The new scaffold is already clean on all of these — 0 whitespace defects, `PATNO` is
`int64`, 0 duplicate keys across 19,450 rows — but the vendor files are the historical
offenders, and this is where the safety belongs.

`MERGE_INDEX` is retained as-is for continuity with the existing `*-to_merge.tab`
files; `key` is built independently and the two are asserted equal where both exist.

### 5. Per-block PCs — 14 blocks, regenerated

Per instruction, per-block PCs as before: **12 panel blocks** (`p277_CSF`,
`p282_{CNS,Inflammation}_CSF`, `p288_{CNS,Inflammation}_plasma`, `p293_olink_plasma`,
`p312_{Inflammation,Neuro}_{CSF,Plasma}`, `p314_{CSF,Plasma}`) plus
`All_Combined` and `All_Combined_no_p293` = **140 columns**.

> **Superseded:** the All_Combined blocks were dropped. 12 blocks / 120 columns shipped.

Method unchanged: earliest visit carrying data for that prefix → drop assays missing
in ≥20% → mean-impute residual → `StandardScaler` → `PCA(n_components=10)`, scores
joined back on `PATNO` so baseline PCs propagate to all visits. Eigenvectors and
loadings written to `EDA_PCA_plots/` under the run timestamp.

Two things change by necessity:

- **`All_Combined*` now includes p312/p314.** The p312/p314 merge deliberately left
  these stale to avoid perturbing existing analysis columns. A from-scratch rebuild
  has no such constraint, and leaving them stale would mean shipping PCs that ignore
  half the dataset. Regenerating is the coherent choice — flagged because it *will*
  move every result derived from those columns.
- **Every PC block shifts anyway**, because the PCA sample sets grow with the new
  scaffold (p312's inputs roughly quadruple). PC values are not comparable to the
  current release even for untouched blocks; expect sign flips and reordering.

Naming stays `*_earliest_visit_PC*`. Do **not** reintroduce `*_baseline_PC*` — and
note the legacy `GroupBy.first()` defect documented in the p312/p314 README (it took
the first non-null value *per column*, not the first row, so "baseline" rows were
composites). The rebuild recomputes from source and does not inherit it.

---

## Derived variables — port, don't rewrite

`recode_analysis_groups.py` (31 KB) is the reference implementation. Inventory to
reproduce, **54 columns** (an earlier estimate of ~66 in this document was wrong):

| Family | Count | Inputs |
|---|---|---|
| `grp_*` contrast flags | 30 | `subgroup`, `analytic_subgroup`, `NSD_Status`, `NSD_STAGE`, `COHORT` |
| `tte_*` / `event_*` pairs | 16 | `moca`, `cogstate`, `pm_*`, `Stage_D`, `NSD_STAGE`, `YEAR` |
| `slope_*` | 3 | `moca`, `updrs3_score`, `lowput_ratio` |
| `CI_PI`, `CI_MOCA` | 2 | cognitive impairment flags |
| `lowput_ratio` | 1 | **MIA_\* now** — see §1 |
| `time_to_LEDD_years`, `disease_duration_years` | 2 | `LEDD`, `ageonset`, `age_at_visit` |

All group definitions are computed **at baseline index** and written per-PATNO, so the
larger scaffold changes group sizes but not definitions. Expect every `grp_*` N to grow.

---

## Build workflow

Scripts land in this directory (`phase_1_wrap/`) for the eventual GitHub push.

```bash
python3 build_step_1_clinical_scaffold.py   # curated xlsx -> scaffold + key + hygiene audit
python3 build_step_2_format_blocks.py       # per-project QC filters -> *-to_merge.tab
python3 build_step_3_genetics.py            # GP2_* (retained) + p9001_* (new)
python3 build_step_4_merge.py               # left-join all blocks onto scaffold on key
python3 build_step_5_derive.py              # port of recode_analysis_groups.py
python3 build_step_6_block_pcs.py           # 14 PC blocks -> 140 columns
python3 build_step_7_dictionary.py          # .tab + .xlsx, 1:1 assertion
python3 build_step_8_compare.py             # vs current release
```

**Step 1** reads sheet `20260511`, builds `key`, runs the §4 audit, and writes a
scaffold + audit report. Row count here is the invariant every later step asserts
against — 19,450 in, 19,450 out.

**Step 2** dispatches per project, reusing the verified `*-to_merge.tab` files as the
default path (see [Source availability](#source-availability)). It re-asserts each
block's filter provenance, checks for duplicate `MERGE_INDEX`, and builds `key`
independently of `MERGE_INDEX` so the two can be cross-checked.

**Step 4** left-joins onto the scaffold so the clinical baseline is preserved exactly,
matching every previous merge. Join artifacts (`PATNO_p312inflam`, …) are retained and
documented as `Technical` stubs, as before.

**Step 7** writes both dictionary formats. The `.xlsx` is not a convenience copy — 44
cells contain embedded newlines that Excel's text import breaks across rows, so the
`.tab` is for programmatic use and the `.xlsx` for hand inspection. One row per
column, bold frozen header, autofilter, wrap-text widths.

Helpful: the curated cut ships its **own data dictionary sheet in the identical
9-column schema** (`Category, Variable, Description, Code, Decode, Derived Variable,
Original Variable(s), Original Dataset(s), Derivation Notes`) — 448 rows. It drops
straight into the clinical layer. Note it is stored sparsely: multi-row `Code`/`Decode`
blocks carry `NaN` in `Variable` and must be collapsed to one row per variable.

---

## Validation plan

Ordered by what would actually catch a defect.

1. **Row invariant** — 19,450 scaffold rows in, 19,450 out of every step. Assert, don't print.
2. **Key integrity** — zero duplicate `key`; every block's key set ⊆ scaffold key set; hygiene audit clean.
3. **Dictionary 1:1** — every dataset column has exactly one dictionary entry and vice versa, asserted on write. This caught a real defect before (`p288_*_PLATE_ID` documented as analytes because a `startswith` test shadowed the plate-ID branch — verify both are `Technical`/`PlateID`).
4. **Value identity on the unchanged path** — for the blocks reused from verified `*-to_merge.tab`, values joined to shared keys must be bit-identical to the current release. Any drift means the join changed, not the data.
5. **Coverage delta table** — per block, samples matched old vs new, reproducing the headline table from the actual build rather than from this pre-flight estimate.
6. **Per-block sparsity + participant-level flags** — port `flag_missing_ids_312_314.py` to all blocks. Analyte columns only: `*_PLATE_ID` is technical and `*_earliest_visit_PC*` propagates, so both would report coverage where no measurement exists. `matched_but_all_missing == 0` is the check that matters.
7. **Derived-variable reconciliation** — every `grp_*` / `tte_*` N, old vs new, with growth explained by the participant expansion. `lowput_ratio` gets its own before/after distribution comparison given §1.
8. **Dictionary diff** — added/removed/recategorized entries vs `Project_9004_Data_Dictionary.tab` (23,821 entries, current).

---

## Decisions taken

1. **Both genetics blocks are kept.** `GP2_*` is retained as-is and `p9001_*` is added
   alongside it, rather than replacing it.
2. **`p9001_*` naming: prefix only, no other transformation.** Columns from
   `GP2_R12_PPMI_PRS_PCs.csv` get `p9001_` prepended to the source name verbatim —
   `p9001_Genetic_PRS_PRS157`, `p9001_Genetic_PRS_InfPop`, `p9001_Genetic_PRS_PC1`,
   `p9001_rs356182_A`. No other block's names change. (The doubled `p9001_Genetic_PRS_`
   is verbose; it is a one-line change in step 3 if it grates once seen in the
   dictionary.)
3. **Per-block PCs are regenerated**, all 14 blocks / 140 columns — see §5.
4. **`lowput_ratio` switches to the MIA series** — no alternative exists in the new cut.
   Documented in `Derivation Notes`; the 0.75 DaT-positivity threshold is re-checked
   against the MIA distribution rather than assumed. See §1.
5. **`All_Combined*` PCs are regenerated including p312/p314.** Consistent with a
   from-scratch build; flagged because it moves every downstream result derived from
   those columns.

### Resolved during the audit

- **p277 and p282_CNS are not blocked** — their filtered files were in the nested
  subdirectory. An earlier draft of this plan flagged them as blocked and warned that
  the p277 fallback would forfeit a large coverage gain. Both claims were wrong: p277's
  vendor file holds 1,200 samples and already matches at 97.2% on *both* cuts, so there
  was no gain to forfeit.
- **p293 does not gain 66 samples.** An earlier draft inferred that from the release
  join; measured against the new cut it is 784 (92.2%) on both.

### Still open

- **Only p282 Inflammation lacks a raw source** (report `.docx` only). Its filtered file
  is present and used. Worth sourcing
  `PPMI_Project_282_NULISAseq_Inflamation_Panel_NPQCounts_20260120.xlsx` for
  completeness of the archive, but it does not block anything.

---

## Out of scope (phase 2)

Per the notes, deferred until the data foundation is asserted solid:

- Combining projects by assay and analyte type
- Re-running the regression / meta-analysis batch and the app
- Adding p312/p314 (and p9001) prefixes to `PROTEOMIC_PREFIXES` in `regressions.py` —
  still hard-coded to `p277_CSF`, `p282_*`, `p288_*`, `p293_olink_plasma`, so the new
  projects are invisible to the batch until edited
- Whether to move `batch.yaml` strata from `GP2_nba_label` to `p9001_InfPop`
