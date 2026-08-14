# Project 9004 — unified dataset rebuild (phase 1 wrap)

Rebuilds `Project_9004_Unified_Emerging_Biomarkers.tab` and its data dictionary from
scratch on a refreshed clinical scaffold, following the filters established in
the MJFF proteomics EDA README and the p312/p314 build README. Both are working
documents in the restricted data directory and are not part of this repository.

| Output | Shape |
|---|---|
| `Project_9004_Unified_Emerging_Biomarkers.tab` | 19,450 rows × 23,932 columns |
| `Project_9004_Data_Dictionary.tab` | 23,932 entries × 9 columns, 1:1 with the dataset |
| `Project_9004_Data_Dictionary.xlsx` | same content, one row per entry |
| `comparison_to_previous_release.md` | full diff against the previous release |

Planning document: [PHASE1_REBUILD_PLAN.md](PHASE1_REBUILD_PLAN.md).

## What actually changed

**Only two projects gained data. Every pre-existing proteomic value is unchanged.**

The proteomic blocks are reused from their documented QC-filter outputs
(`*-to_merge.tab`), so the analyte layer is identical-by-construction. This is verified
exhaustively, not sampled — [verify_bit_identical.py](verify_bit_identical.py) checks
**every** analyte column on **every** shared key, and checks the null patterns too
(a value present in one release and absent in the other is a discrepancy even though no
value "differs"):

```
TOTAL   23,335 columns   27,507,875 value-carrying cells   0 value diffs   0 presence diffs

PASS — all 23,335 analyte columns are identical across 27,507,875 value-carrying cells
and 10,349 shared keys, with matching null patterns.
```

What moved is the clinical scaffold underneath them —
`PPMI_Proteomics_Data_Cut_INTERNAL_20251215` (10,354 rows / 2,025 participants) replaced
by `PPMI_Curated_Data_Cut_Public_20260511` sheet `20260511` (19,450 / 4,788):

| Block | Vendor samples | Matched, previous cut | Matched, new cut | Gain |
|---|---|---|---|---|
| p312 Inflammation | 2,204 | 478 (21.7%) | **1,995 (90.5%)** | **+1,517** |
| p312 Neuro | 2,187 | 471 (21.5%) | **1,978 (90.4%)** | **+1,507** |
| p314 CSF | 2,128 | 1,236 (58.1%) | **2,098 (98.6%)** | **+862** |
| p314 Plasma | 3,399 | 1,561 (45.9%) | **2,942 (86.6%)** | **+1,381** |
| p277, p282 ×2, p288 ×2, p293 | — | — | unchanged | +0 |
| **total** | | 15,212 | **20,479** | **+5,267** |

The six older blocks were already matched at 92–99% and do not move — they were assayed
on the cohort the internal cut was built around. Project 312 goes from a 21.7% salvage
to essentially complete, which is the result that justifies the rebuild.

The refresh is non-destructive: all 2,025 previous participants are present, only 5 of
10,354 previous keys drop, and the `EVENT_ID` vocabulary is identical so no visit
remapping is needed.

## Pipeline

```bash
pip install -r requirements.txt

python3 build_step_1_clinical_scaffold.py   # curated xlsx -> scaffold + key + hygiene audit
python3 build_step_2_format_blocks.py       # validate 10 blocks, coverage vs scaffold
python3 build_step_3_genetics.py            # GP2_* (retained) + p9001_* (new)
python3 build_step_4_merge.py               # left-join everything onto the scaffold
python3 build_step_5_derive.py              # 54 derived analysis variables
python3 build_step_6_block_pcs.py           # 12 PC blocks -> 120 columns
python3 build_step_7_dictionary.py          # dataset + dictionary (.tab and .xlsx)
python3 build_step_8_compare.py             # diff vs the previous release
```

Intermediates land in `build_intermediates/` under one shared timestamp; PCA
eigenvectors and loadings in `EDA_PCA_plots/`. Each step asserts the scaffold row count
(19,450 in, 19,450 out) and aborts rather than emitting a silently wrong join.

## Validation

Three layers, three different expectations. Run after a build:

```bash
python3 verify_bit_identical.py   # proteomic layer — must be identical
python3 validate_derived.py       # clinical + derived layers — must be explainable
python3 build_step_8_compare.py   # full narrative diff -> comparison_to_previous_release.md
```

| Layer | Expectation | Result |
|---|---|---|
| Proteomic (23,335 cols) | identical on shared keys | **0 value diffs, 0 presence diffs** across 27.5M cells |
| Clinical (195 shared cols) | small revisions — different cut | **99.9573%** agree (606 of 1,417,951 cells) |
| Derived (54 cols) | moves only where inputs moved | every disagreement **attributed** |

The derived layer is the subtle one, so it is checked by attribution rather than by
equality. Two things legitimately move a derived value: the participant's visit history
changed, or the clinical inputs changed. `validate_derived.py` separates both:

- **Participants who gained or changed visits** (569 of 2,025): the visit-dependent
  columns (`slope_*`, `tte_*`, `event_*`) differ, which is correct — they are functions
  of the whole visit history. The 29 columns that are *not* visit-dependent
  (`grp_*`, `CI_*`) show **0 differing cells**, which is the check that matters.
- **Participants with an identical visit set**: every remaining disagreement traces to a
  participant whose clinical inputs also changed between the two cuts. **0 unexplained.**

Comparing visit *sets*, not visit counts, matters here: participants exist with four
visits in both releases where a year-4 visit (`pm_any` null) was replaced by a year-5
visit (populated), which correctly moves three censoring times. Counting visits would
have mislabelled those as unexplained defects.

`lowput_ratio`, `slope_lowput_ratio` and `grp_NMC_*` are excluded from the reproduce
test and reported separately — they changed by design, as described below.

### `key` — the merge safety belt

Every join goes through `key = PATNO + '_' + EVENT_ID`, built by `make_key()` in
[build_common.py](build_common.py) and asserted equal to the legacy `MERGE_INDEX`.

The failure mode being guarded against is silent: a stray space or an Excel float turns
`####_BL` into `#### _BL` or `####.0_BL`, the left join misses, and the block arrives
100% missing without anything erroring. So defects that are unambiguous typos are
repaired **and reported** — surrounding whitespace, non-breaking and zero-width
characters, float-formatted IDs, lowercase visit codes. Defects that would require
guessing — non-numeric PATNO, malformed visit code, duplicate key — abort the build.

All ten blocks, both genetics files and the scaffold came through clean.

## Source availability

All blocks have usable inputs. Note the nested
`initial_build_assets/MJFF_proteomics-EDA/MJFF_proteomics-EDA/` subdirectory — it holds
the p277, p282_CNS and both p288 filtered files plus the p282 CNS raw xlsx, and is easy
to miss because the parent has near-identical contents.

**p293 is a trap.** `data_from_sources/` holds the 20250328 parquet, but the pipeline
moved to 20251121 in Jan 2026. Re-deriving p293 from the parquet on disk would silently
regress it one data release. The filtered file is current (verified value-identical to
the shipped release) and is what the build uses.

Only p282 Inflammation lacks a raw source (report `.docx` only); its filtered file is
present and used.

## Changes beyond the scaffold

### `lowput_ratio` — a pre-existing bug, fixed

The reference implementation computed
`min(DATSCAN_PUTAMEN_L, R) / lowput_expected`. But `lowput_expected` is documented in
the source dictionary as the *"age-/sex-expected lowest putamen **ratio**"* — it is
already observed/expected. Dividing by it a second time cancels the imaging term and
leaves a function of age and sex.

Evidence in the previous release:

- `corr(lowput_ratio, observed min-putamen SBR) = 0.055` — no imaging signal
- 0 of 4,312 values below 0.75, despite the dictionary entry stating *"Ratios <0.75
  indicate DaT positivity"*; `lowput_expected` itself had 75.8% below 0.75
- `slope_lowput_ratio` was near-constant across participants (sd 0.002), tracking the
  −0.01397/yr age coefficient rather than anyone's dopaminergic decline

The rebuild takes `MIA_LOWPUT_EXPECTED` directly, which the curated cut defines
explicitly as the ratio. Result: correlation with observed min-putamen SBR **0.986**,
and 62.5% of values below 0.75.

**Any prior result using `lowput_ratio` or `slope_lowput_ratio` should be re-read.**
The full derivation and evidence are in the dictionary entry for both variables.

Note the DaTscan series also changed: the curated public cut carries only the MIAKAT
(`MIA_*`) reconstruction, not the `DATSCAN_*` series, and the two are different
quantifications (r ≈ 0.89–0.90) — so values would not have reproduced regardless.

### New `p9001_*` genetics block, alongside the retained `GP2_*`

Both are kept. `GP2_*` rebuilds unchanged (3,359 participants, 14 columns).
`GP2_R12_PPMI_PRS_PCs.csv` adds 4,985 participants × 172 columns: 4 PRS variants, an
inferred-ancestry label, 10 within-population PCs, and 157 SNP dosages. `p9001_` is
prepended to each source name verbatim; genetics is static (`EVENT_ID` is `SC`
throughout) so it joins on `PATNO`.

p9001 covers **3,676** scaffold participants against GP2's 2,620. Where both exist the
ancestry labels agree on **100.00%** of 2,574 participants. Analysis strata remain
defined on `GP2_nba_label` — switching is a phase-2 decision.

### `grp_NMC_*` definition, and the enriched carrier cohort

The refreshed cut introduces 12 new `subgroup` values, including
`LRRK2 + GBA + Normosmic`. The original rule excluded double carriers by exact-matching
the single label `LRRK2 + GBA`, so the new labels would have been counted as clean
single-gene carriers in several NMC groups at once. Exclusion is now by gene-token test,
which preserves the stated intent and leaves behaviour on the old vocabulary unchanged
(gene+phenotype labels like `LRRK2 + RBD` still match).

Separately, the curated cut brings in the **enriched LRRK2/GBA cohort**, which the
proteomics-selected internal cut excluded. LRRK2-carrying prodromal participants at
baseline go from **6 to 239**. The NMC contrasts existed as columns before but had too
few participants to fit; they are now analyzable. Two cautions carried into phase 2:
the enrichment cohort is not proteomically assayed at the same rate, and adding it
shifts the genetic composition of every contrast that does not condition on carrier
status.

### PC blocks — 12, not 14

`<prefix>_earliest_visit_PC1..PC10` for the 12 proteomic panels = 120 columns. PCs are
computed from **proteomic analytes only** — no clinical, genetic or `PLATE_ID` column
enters the PCA.

`All_Combined` / `All_Combined_no_p293` are **not** built. The previous construction
required a sample to carry data in every project, which across six projects leaves 100
participants; and the pooled matrix is not on a common scale anyway (NPQ vs NPX,
plate-control vs intensity normalization, differing panel versions). A cross-project
block belongs **downstream of harmonizing projects by assay and analyte type**, where
analytes are z-scored within project first — see phase 2.

The "earliest visit" step takes the earliest qualifying **row**, not `GroupBy.first()`.
The legacy pipeline used `.first()`, which takes the first non-null value *per column*,
so a participant's "baseline" was a composite stitched from several visits while
`EVENT_ID` still read `BL`. That defect is not inherited.

PC values are not comparable to the previous release even for untouched blocks — every
block's PCA sample set grew, so expect sign flips and component reordering.

### Clinical columns: 59 out, 17 in

Dropped: 12 assay-availability / near-visit indicators (`NULISA_*`, `Olink_*` — they
powered a near-visit rescue mechanism the p312/p314 audit reported on but never
applied), 32 SAA replicate/detail columns (headline `CSFSAA` retained), and the 15
DaTscan columns above.

Added: `Clinical_Stage`, `Death_Status`, `age_death`, `Death_Date`, `CGI`, `PGI`,
`PDAQ27`, `nqol_*`, `NFL_CSF_ULOD`, `Roche_IL_1B_CSF`, `Roche_IL_6_CSF(_LLOD)`,
`Roche_YKL40_CSF(_LLOD)`.

No group-defining variable is lost — every derived flag keys off `subgroup`,
`analytic_subgroup`, `NSD_Status`, `NSD_STAGE`, `COHORT`, all present.

## Reading the outputs

**Values are on native, project-specific scales and are not comparable across
projects.** NPQ vs NPX, different normalizations, different panel versions, no formal
bridging. Any analysis putting analytes from two projects side by side must z-score
each analyte within its own project first. This carries over unchanged — see the
warning section of the p312/p314 README.

**For analysis, read the N column, not the % column.** The pipeline agrees: every gate
in `../regressions.py` and `batch.yaml` is an absolute row count (`min_n: 20` rows,
`min_groups_per_param: 2`) applied after `dropna()` on the design columns. There is no
fraction-based coverage floor anywhere in the batch, so a block whose N is unchanged
fits exactly as it did before — the higher percentage is cosmetic.

The one place to be careful is any *ad-hoc* analysis that filters on a missingness
fraction. A threshold tuned on the previous release will be materially stricter here,
and it will bite the older blocks hardest — precisely the blocks that did not change.

**Use the `.tab` for programmatic work and the `.xlsx` for anything opened by hand.**
42 dictionary cells contain embedded newlines. They are correctly quoted in the `.tab`,
so any standards-compliant TSV reader handles them — but Excel's text import does not
honour the quoting and splits those entries across rows. The `.xlsx` writes them as real
cells with the newline inside, one row per entry.

## Harmonized proteomic blocks

Six blocks combine the two projects that ran the same platform × panel × biofluid, so an
analyte measured twice under the same assay becomes one column spanning both projects.
Project-prefixed columns are **retained** — harmonized columns are added alongside, never
in place of, so every result stays checkable against the unpooled data.

| Block | Reference (larger) | rows | Mapped (smaller) | rows | Core analytes |
|---|---|---|---|---|---|
| `olink_plasma` | p314_Plasma | 2,942 | p293_olink_plasma | 784 | 5,400 |
| `olink_csf` | p314_CSF | 2,098 | p277_CSF | 1,167 | 5,416 |
| `nulisa_cns_plasma` | p288_CNS_plasma | 2,665 | p312_Neuro_Plasma | 1,618 | 127 |
| `nulisa_cns_csf` | p282_CNS_CSF | 2,282 | p312_Neuro_CSF | 1,197 | 128 |
| `nulisa_inf_plasma` | p288_Inflammation_plasma | 2,669 | p312_Inflammation_Plasma | 1,616 | 249 |
| `nulisa_inf_csf` | p282_Inflammation_CSF | 1,899 | p312_Inflammation_CSF | 1,234 | 247 |

**The reference is whichever project has more rows**, not whichever is older. That
maximises the share of each harmonized column that is uncorrected native value and
minimises exposure to correction error. On these data it makes p314 the reference for
both Olink blocks; the four NULISA blocks keep the older project, which is also larger.

`nulisa_cns_*` pairs p282/p288's **CNS Disease** panel with p312's **Neuro 220** panel —
confirmed as counterparts by analyte overlap (128 of p282's 129, 127 of p288's 131).

### Naming

```
harmonized_<block>_<ANALYTE>_<NPX|NPQ>     the harmonized value
harmonized_<block>_src                      which project supplied this row's value
harmonized_<block>_earliest_visit_PC1..PC10 PCs on the harmonized block
collection_era                              categorical collection period (see p277 below)
```

Harmonized names carry the **true** metric — `NPQ` for NULISA, `NPX` for Olink. Note
p288's 380 legacy columns are suffixed `_NPX` despite being NULISA NPQ values
(`prepare_for_merge.py:19` builds the name while line 25 pivots `values='NPQ'`). Those
legacy names are left untouched for continuity, so `p288_CNS_plasma_ACHE_NPX` and
`harmonized_nulisa_cns_plasma_ACHE_NPQ` refer to the same metric under different
suffixes. The dictionary records this for every affected entry.

### The correction

The mapped project is placed on the reference's scale per analyte:

```
mapped_corrected = slope × mapped + intercept
```

fitted robustly on **(EVENT_ID × COHORT) matched cells** carrying ≥10 samples per
project — not on raw pooled values. That matters because cohort composition differs
sharply between projects at shared visits: at V10 in `nulisa_cns_csf`, p282 is 92% PD /
0% HC against p312's 12% PD / 40% HC. Within each cell the median and IQR of each
project are taken, then combined across cells weighted by the smaller cell count:

```
slope     = IQR_reference / IQR_mapped
intercept = median_reference − slope × median_mapped
```

The reference value always wins; the corrected mapped value fills **only** rows the
reference did not cover. Values are not averaged, so every cell traces to one assay run.

Values stay on native NPQ/NPX scales. **They are deliberately not z-scored** — see
[Why not z-score](#why-not-z-score) below.

### Validation — six analyses, and what each established

Full outputs: [harmonization_input_validation.md](harmonization_input_validation.md),
[harmonization_block_descriptives.md](harmonization_block_descriptives.md),
[harmonization_corrections.md](harmonization_corrections.md),
[harmonization_positive_controls.md](harmonization_positive_controls.md).

**1. Are the projects on the same scale?** Yes, all six. Per-analyte median differences
run −0.08 to +0.31 log2 against a within-analyte SD of ~0.8–1.0, and IQR ratios 0.85 to
1.28. For reference, the NPQ↔NPX gap between *platforms* is ~11.6 log2 — two orders of
magnitude larger, which is why platforms are never combined.

**2. Is the difference a batch offset or a visit effect?** A batch offset. Project is
partly collinear with visit (p282/p288 are the baseline arm, p312 the follow-up arm),
so the offset was re-estimated at every visit with ≥10 samples in both projects. It is
flat across 3–8 visits per block — spread 0.035 to 0.119 log2. A progression artifact
would drift with visit; these don't. Within-project level drift is ≤0.11 log2 at every
visit, which is what licenses extrapolating the correction to visits only one project
covers (the three CSF blocks have no shared visits past V10–V12).

**3. Is it cohort composition rather than assay?** No. Direct standardization within
(visit × cohort) cells moves the offsets by at most **0.018 log2**. The offset is a
median over hundreds-to-thousands of analytes, and disease shifts a minority of them, so
the median analyte is composition-insensitive.

**4. Do the corrections generalize?** Yes. Fitted on half the matched cells and evaluated
on the held-out half, the location gap falls **64–95%** and the dispersion ratio moves
from 0.88–1.27 to **1.00–1.03**. Intercept alone fixes location but leaves dispersion
untouched by construction — that residual spread mismatch is exactly the
heteroscedasticity-tracking-project artifact that corrupts longitudinal models, which is
why slope is applied and not just intercept.

**5. Are the overlapping samples genuine replicates?** Only in one block. Where the same
PATNO_EVENT_ID appears in both projects, a sample should match itself across projects
better than it matches other samples once each analyte is centred:

| Block | overlapping visits | self-match r | self ranked #1 | verdict |
|---|---|---|---|---|
| `olink_plasma` | 56 | **+0.647** | **54 / 56** | genuine replicates |
| `olink_csf` | 30 | +0.074 | **4 / 30** | **not replicate-grade** |
| four NULISA blocks | 0 | — | — | none exist by design |

Restricting to analytes with SD ≥ 1.0 in both projects — the attenuation test — raises
`olink_plasma`'s median r² from 0.303 to **0.493**, exactly as real replicate pairs
diluted by flat analytes should behave. `olink_csf` moves the *wrong* way (0.032 →
0.020), which no attenuation model produces, so its 30 overlaps are separate draws
sharing a visit code rather than split aliquots and are excluded from fitting. The four
NULISA blocks have no overlap at all: p282/p288 and p312 share 149–277 participants but
never the same visit, so their corrections come from matched groups instead.

**6. Does harmonization preserve known biology?** Yes — this is the test that matters.
Contrast `grp_NSD_vs_HC`, OLS with PATNO-clustered SEs over all visits, outcome
z-scored. A random-intercept LMM is *not* usable here: these panels average ~1.2
observations per participant, so the variance component is unestimable and the fit is
singular.

DDC is the anchor — the strongest established signal in this dataset, 29 robust hits
across 7 panel columns in the previous batch — and it behaves exactly as a working
harmonization should, holding beta while N and precision rise:

| Block | reference beta (P) | harmonized beta (P) | observations |
|---|---|---|---|
| `nulisa_cns_csf` | +1.307 (4.7e-74) | +1.249 (**4.5e-146**) | 1,774 → 2,382 |
| `olink_csf` | +1.175 (2.2e-54) | +1.009 (**1.3e-74**) | 1,199 → 2,312 |
| `olink_plasma` | +0.545 (1.4e-15) | +0.554 (**3.3e-26**) | 1,624 → 2,326 |
| `nulisa_cns_plasma` | +0.373 (7.9e-6) | +0.351 (**1.1e-8**) | 2,120 → 2,744 |

Across all 32 positive-control fits: median SE fell 27% (0.083 → 0.061), p-values
improved in 18 of 32, and **no beta shifted further than sampling noise explains**
(|z| > 1.96 in 1 of 32). Naive pooling without correction lands 2.1× further from the
reference estimate than harmonization does.

**Harmonization is not a free p-value win.** Among the 15 fits where the reference was
already significant, only 8 improved — a small beta shift toward the mapped estimate
partly cancels the precision gain. It buys generalizability and precision, not uniformly
smaller p-values.

### p277 carries a collection-era confound — found by the negative controls

The spike-in controls are the pipeline's built-in null, and one failed:

| Control | beta | P |
|---|---|---|
| **`p277_CSF_CTRL`** | **−0.262** | **3.9e-04** |
| `p293_olink_plasma_CTRL` | +0.055 | 0.62 |
| `p314_CSF_CTRL` | −0.010 | 0.90 |
| `p314_Plasma_CTRL` | +0.053 | 0.53 |
| `p312_*_mCherry` ×4 | — | *zero variance, unusable* |

A spike-in cannot associate with disease status, so p277 has internal structure tracking
diagnosis. The cause is collection era: p277's post-2020 samples contain **2 healthy
controls against 372 NSD+**. It is panel-wide, not confined to the control — under an
identical model p277 shows **λ = 2.486** against p314's 1.113.

| p277 subset | n (HC / NSD+) | CTRL P | λ |
|---|---|---|---|
| all samples | 1,142 (293 / 849) | 3.9e-4 | 2.486 |
| pre-2020 | 768 (291 / 477) | 0.097 | 1.506 |
| 2020+ | 374 (**2** / 372) | 1.7e-2 | **57.9** |

**This is a p277 data problem, not a harmonization problem** — it affects existing
p277-based results independently of anything done here, and is worth raising separately.

**Adjust, do not exclude.** This is a positivity violation: only 6 of 15 collection years
carry both arms at n≥10, so in the rest the within-year contrast is not estimable. With a
dummy per year those samples' means are absorbed and they contribute nothing to the group
coefficient — categorical adjustment is algebraically equivalent to restriction for the
contrast, while keeping the samples for nuisance terms and adapting per contrast. The 374
post-2020 samples are unusable for an HC-referenced contrast but perfectly good for
within-PD or NSD-stage analyses, so deleting them would forfeit that.

| Strategy | n | CTRL P | λ |
|---|---|---|---|
| no adjustment | 1,142 | 3.9e-4 | 2.486 |
| **+ categorical era** | 1,142 | 5.2e-2 | **1.459** |
| + *linear* year | 1,142 | 3.7e-2 | 1.731 |
| restrict to pre-2020 | 768 | 9.7e-2 | 1.506 |

**Era must be categorical.** A linear year term is worse than the categorical version
(λ 1.731 vs 1.459) because it extrapolates a trend through the single-arm years instead
of neutralising them.

`collection_era` is therefore added to the dataset as a categorical variable, and no
samples are dropped. Applying it as a covariate is an analysis-time decision for
`batch.yaml`, not a build-time one.

Two caveats: neither strategy fully resolves it — λ stays 1.43–1.51 against p314's 1.11 —
and era is only a proxy for what is presumably assay batch or run date, which p277 does
not carry (it has no `PLATE_ID`). p277-driven findings warrant extra scrutiny regardless.

### Why not z-score

Z-scoring each project before combining would force offset = 0 and IQR ratio = 1 for
every analyte, making all six blocks "harmonize" perfectly and unverifiably — every
diagnostic above would become vacuous, since the transform guarantees the result rather
than measuring it.

It would also destroy real signal. NEFL differs by **−2.87 log2** between p282 and p312
in CSF; p312's samples sit 2–8 years later in follow-up and NEFL is the canonical
axonal-damage marker. That difference is either genuine progression or a real
assay-version difference, and z-scoring erases it either way with no way to tell which.
A +0.06 offset correction leaves it intact, which is correct — that difference is not
batch. Between 43% and 50% of analytes in the NULISA and `olink_plasma` blocks carry
per-analyte structure beyond the common offset that z-scoring would flatten.

Z-scoring does belong inside the PCA (`StandardScaler`, already there) and in any
cross-platform comparison, where NPQ and NPX genuinely cannot be compared raw.
`regressions.py` also z-scores each predictor per run and stratum, so dataset-level
z-scoring would be redundant for the batch regardless.

### Known limitations

- **`olink_csf`'s correction is the least evidenced** — no replicate-grade overlaps, and
  its mapped project (p277) carries the era confound above. Its positive controls behave
  worst of the six: NEFL loses significance (8.9e-3 → 0.77) and GFAP/MAPT/CHI3L1/SNCA
  move toward p277's estimates.
- **p312's negative controls are unusable** (zero variance), so the four NULISA blocks
  cannot be checked this way.
- **Slope estimates rest on 4–12 matched cells.** Bootstrapping shows 32–62% of
  nominally-flagged slopes have CIs spanning the tolerance band, so slope is applied only
  where the bootstrap CI supports it and set to 1 otherwise.
- **Harmonized columns are a union, not an average.** The two projects overlap in only
  3–11% of participants, so `harmonized` = reference participants + mapped-only
  participants. A harmonized effect is therefore not guaranteed to lie between the two
  project-specific estimates.
- **Blocks differ in correction dependence**: `olink_plasma` draws 87% of its harmonized
  values from corrected data, `nulisa_cns_csf` only 24%. The per-block share is recorded
  in the dictionary.

## Re-running the analyses

[batch.yaml](batch.yaml) is the phase-1-wrap analysis config, derived from the
p312/p314 batch config. 52 runs.

**Status — the config is validated on one run; 51 remain to be run.** `NSD_vs_HC` was
executed over both strata on 2026-08-13 as a smoke test of the new config against the
rebuilt dataset. It completed clean, so what is left is compute time rather than
configuration risk:

| Stratum | Predictors offered | Rows written | Numerical failures | λ overall | Run-wide sig | Family sig |
|---|---|---|---|---|---|---|
| EUR | 34,916 | 34,910 | 2 | 1.232 | 14 @ 1.43e-06 | 25 |
| AJ | 34,916 | 23,900 | 33 | 1.274 | 4 @ 2.09e-06 | 8 |

AJ writes fewer rows because analytes that cannot meet `min_n: 20` in the smaller
stratum are dropped before fitting, not because anything failed — the same reason
`p277_CSF` and `p293_olink_plasma` produce no AJ family at all.

Inflation is concentrated in the small panels: every 5,400-analyte panel lands between
1.07 and 1.38 in both strata, while every family above 1.9 carries 99–249 tests, where λ
is estimated from too few points to be stable. The extremes are 2.24
`harmonized_nulisa_cns_plasma` (127 tests) in EUR and 2.57 `p282_CNS_CSF` (99 tests) in AJ.

**`p277_CSF` comes in at λ = 1.16 in EUR**, against the 2.49 recorded for that panel in
the unadjusted era diagnostic above. Not a controlled before/after, but it is the
direction `collection_era` was added to produce.

Both strata together took 13.1 minutes wall-clock, which is the basis for the ~11 hour
serial estimate below.

The analysis chain is **copied from the parent directory** so this directory is
self-contained and the pushed repository is runnable end to end:

| Script | Status |
|---|---|
| [regressions.py](regressions.py) | **modified** — harmonized prefixes, PRS157, proteome-wide lens |
| [meta_analysis.py](meta_analysis.py) | copied unchanged |
| [build_results_summary.py](build_results_summary.py) | copied unchanged |
| [calibration_check.py](calibration_check.py) | copied unchanged |

The parent's copies remain canonical; if they change later, re-sync. `plot_results.py`
was not copied and still lives only in the parent.

```bash
python3 regressions.py --config batch.yaml     # 52 runs -> results/
python3 meta_analysis.py                       # pool EUR + AJ -> meta/
python3 calibration_check.py --results-dir results
python3 build_results_summary.py
```

`results/` holds one CSV per run × stratum — 7–10 MB each, one row per analyte carrying
`beta`, `SE`, `P`, `N` and the three significance lenses. They are summary statistics,
not participant-level records, but they are large and fully regenerable, so they are not
committed.

### What changed from the previous batch

| | Previous | Now |
|---|---|---|
| Input | `../unified_PPMI-*.tab` | `Project_9004_Unified_Emerging_Biomarkers.tab` |
| Ancestry strata | `GP2_nba_label` | **`p9001_Genetic_PRS_InfPop`** |
| Genetic PCs | `GP2_PC1–5` | **`p9001_Genetic_PRS_PC1–5`** |
| PRS | `GP2_PRS_zscore` | **`p9001_Genetic_PRS_PRS157`** |
| Assay prefixes | 12 project blocks | 12 project + **6 harmonized** |
| Multiple testing | run-wide + per-family | + **proteome-wide over 9,500 proteins** |
| Covariates | — | + **`collection_era`** |

**Ancestry and PRS now come from Project 9001**, which raises stratum coverage from
2,620 to **3,676** participants — EUR 1,920 → 2,825, AJ 582 → 669.

`PRS157` is the right counterpart to the retired `GP2_PRS_zscore`: it correlates with it
at **r = 0.992**, so PRS-interaction results stay comparable to earlier batches while
gaining ~1,000 participants. The same file ships PRS154 (r=0.963), and **PRS152 / PRS149
which are a different construct entirely** (r ≈ 0.46–0.49 against the legacy score) —
they are not interchangeable, and the Project 9001 documentation describes the variant
sets in detail.

**Genetics are static and propagate across all timepoints.** The p9001 source file is
`EVENT_ID == SC` only; it is joined on `PATNO` alone, so every one of the 18 visit codes
carries the values (76.8–100% of rows) with zero participants holding conflicting
values. `SC` itself does not appear in the dataset.

### Multiple testing — three lenses

Adding the harmonized blocks alongside the project-specific ones means the same protein
is tested several times: `p282_CNS_CSF_DDC_NPQ` and `harmonized_nulisa_cns_csf_DDC_NPQ`
are heavily correlated, because the harmonized column *contains* the reference project's
values. A raw analyte count is therefore the wrong denominator — it would penalise every
panel for that redundancy.

| Column | Denominator | Use for |
|---|---|---|
| `significant` | analytes in the run | continuity with earlier batches |
| `significant_family` | analytes in that assay block | per-block context |
| **`significant_proteome`** | **fixed 9,500 independent proteins** | **the headline lens** |

`n_independent_proteins: 9500` is set in `batch.yaml` and is invariant to how many
redundant columns a run happens to include. It applies only to high-throughput proteomic
predictors; the targeted/clinical biomarkers keep the family threshold, since they are
not part of the proteome panel set.

### `collection_era` is applied globally, not targeted

The intent was to adjust only where Project 277 contributes. That turned out not to be
expressible: **every run loops all proteomic prefixes**, so `p277_CSF` and
`harmonized_olink_csf` analytes appear in every run, and there is no run-level scope that
isolates them. It is therefore in `default_covariates`.

The cost elsewhere is one degree of freedom; the benefit is that p277's era confound —
which makes its spike-in control significant at P=3.9e-4 and inflates its whole panel to
λ=2.49 against p314's 1.11 — is adjusted wherever it can bite. It **must** stay
categorical; a linear year term is worse than useless (λ 1.731 vs 1.459) because it
extrapolates through the single-arm years. Any run can drop it via `drop_covariates` as
a sensitivity check.

This does mean every result differs from the previous batch by one covariate, so
run-to-run comparisons with the 2026-08-01 batch are not like-for-like.

## Phase 2 (not done here)

- Re-run the regression / meta-analysis batch and the app on the harmonized blocks,
  adding `collection_era` as a categorical covariate to `batch.yaml` where p277 contributes
- Add the p312/p314 and `harmonized_*` prefixes to `PROTEOMIC_PREFIXES` in
  `../regressions.py` — still hard-coded to `p277_CSF`, `p282_*`, `p288_*`,
  `p293_olink_plasma`, so neither the new projects nor the harmonized blocks are visible
  to the batch until edited
- Decide whether to move `batch.yaml` strata from `GP2_nba_label` to the p9001 label
- Raise the p277 collection-era confound with the data providers — it affects existing
  p277 results independently of this build
- Re-read any prior result that used `lowput_ratio` / `slope_lowput_ratio`

## Data size and coverage

`Project_9004_Unified_Emerging_Biomarkers.tab` — **19,450 rows × 35,566 columns**,
2.19 GB. One row per participant-visit, keyed on `key` = `PATNO_EVENT_ID`.
4,788 participants across 18 visit codes (`BL`, `PW`, `ST`, `V04`–`V22`).

| Layer | Columns |
|---|---|
| Clinical scaffold (curated cut) | 212 |
| Project-specific proteomic analytes | 23,335 |
| **Harmonized proteomic analytes** | **11,567** |
| Per-project PCs (12 blocks × 10) | 120 |
| **Harmonized PCs (6 blocks × 10)** | **60** |
| Genetics — `GP2_*` 14 + `p9001_*` 172 | 186 |
| Derived analysis variables | 54 |
| Technical (`key`, `MERGE_INDEX`, plate IDs, join artifacts, `_src`, `collection_era`) | 32 |
| **Total** | **35,566** |

### Coverage per block

Rows and participants carrying at least one analyte value. Harmonized blocks are listed
against the two projects they pool.

| Block | Analytes | Rows | Participants |
|---|---|---|---|
| p277_CSF | 5,416 | 1,167 | 746 |
| p282_CNS_CSF | 129 | 2,282 | 1,909 |
| p282_Inflammation_CSF | 247 | 1,899 | 1,650 |
| p288_CNS_plasma | 131 | 2,665 | 1,861 |
| p288_Inflammation_plasma | 249 | 2,669 | 1,864 |
| p293_olink_plasma | 5,415 | 784 | 304 |
| p312_Inflammation_CSF | 248 | 1,234 | 916 |
| p312_Inflammation_Plasma | 250 | 1,616 | 1,021 |
| p312_Neuro_CSF | 212 | 1,197 | 888 |
| p312_Neuro_Plasma | 221 | 1,618 | 1,020 |
| p314_CSF | 5,416 | 2,098 | 1,648 |
| p314_Plasma | 5,401 | 2,942 | 2,087 |
| **harmonized_olink_plasma** | 5,400 | **3,670** | **2,323** |
| **harmonized_olink_csf** | 5,416 | **3,235** | **2,239** |
| **harmonized_nulisa_cns_plasma** | 127 | **4,283** | **2,732** |
| **harmonized_nulisa_cns_csf** | 128 | **3,479** | **2,520** |
| **harmonized_nulisa_inf_plasma** | 249 | **4,285** | **2,735** |
| **harmonized_nulisa_inf_csf** | 247 | **3,133** | **2,323** |

Every harmonized block covers more participants than either project alone — the largest
gain is `harmonized_olink_plasma` at 2,323 against p293's 304 and p314_Plasma's 2,087.

Genetics: **2,620** participants carry `GP2_*`, **3,676** carry `p9001_*`. Both are
static and propagate to every visit.

### Reading the coverage numbers

**Percent-missing is the wrong statistic here.** The scaffold grew from 10,354 to 19,450
rows, so a block whose absolute N is unchanged shows a higher missing percentage purely
from the larger denominator. Read the rows and participants columns above. The analysis
pipeline agrees — every gate in `regressions.py` and `batch.yaml` is an absolute row
count (`min_n: 20`) applied after `dropna()`, so no fit is affected by the denominator
change. Only an ad-hoc filter written as a missingness *fraction* would tighten.

**Analyte columns are not independent tests.** 34,902 proteomic columns represent far
fewer independent proteins: the same protein appears in several projects and again in
its harmonized block, where the harmonized column literally contains the reference
project's values. This is why the batch reports `significant_proteome` against a fixed
9,500-protein denominator — see [Multiple testing](#multiple-testing--three-lenses).

**Batch runtime.** One run over both strata takes ~13 minutes at 34,916 predictors per
stratum, so the 52-run config is roughly **11 hours** serially. The parent directory's
`run_batch_sharded.sh` shards it across processes if that matters.

## Data availability

**No data is distributed with this repository — code and documentation only.** Every
input and output named here lives behind its originating study's access process, and
nothing in the pipeline will run without them.

| Input | Source | Access |
|---|---|---|
| Clinical scaffold — `PPMI_Curated_Data_Cut_Public_20260511` | Parkinson's Progression Markers Initiative (PPMI), sponsored by the Michael J. Fox Foundation | PPMI data access application |
| Proteomic blocks — Projects 277, 282, 288, 293, 312, 314 | PPMI project data, vendor QC-filtered releases | PPMI data access application |
| Genetics — `GP2_*`, `p9001_*` (release 12 PRS and PCs) | Global Parkinson's Genetics Program (GP2), funded by the Aligning Science Across Parkinson's initiative | GP2 data access application |

The rebuilt `Project_9004_Unified_Emerging_Biomarkers.tab`, its data dictionary and the
`results/` regression outputs are **derivatives of restricted data** and are likewise not
included; they are obtained through the same access path, not from this repository.

Anyone using these sources is bound by the PPMI and GP2 data use agreements, including
their citation and acknowledgment requirements. The counts, coverage tables and
diagnostics reported throughout this README are aggregate summaries published to
document the build — they are not a data release.

## Note on sharing

This directory holds restricted PPMI data, so [.gitignore](.gitignore) is an
**allowlist, not a blocklist**: it ignores everything, then re-admits only the build
code and documentation at the top level.

```gitignore
*
!/.gitignore
!/*.py
!/*.md
!/*.yaml
!/requirements.txt
results/
meta/
plots/
```

A blocklist leaks — it only excludes the extensions someone remembered to name, and any
new data format or stray file in a subdirectory ships by default. With the allowlist,
git does not descend into `build_intermediates/`, `EDA_PCA_plots/`,
`initial_build_assets/` or `clinical_and_p9001_updates/` at all, and anything added
later is excluded until explicitly admitted.

### What ships, and what was checked

Verified against the real tree (498 files) by pointing a scratch git index at this
directory and running `git add -A -n`: **33 files stage, 440 KB total** — 22 `.py`,
8 `.md`, `batch.yaml`, `requirements.txt` and `.gitignore`. Nothing from any
subdirectory, no `.tab`, no `.log`, no `results/`.

Because the repository is **public**, the staged set was also read for disclosure, not
just filtered by extension:

| Check | Result |
|---|---|
| Absolute paths, usernames, home or drive paths | none |
| Emails, URLs, hostnames | none |
| API keys, tokens, credentials | none |
| Participant identifiers (`PATNO`, `PATNO_EVENT_ID` values) | none — one study ID in a code comment and one in this README were removed |
| Participant-level rows or measured values | none — every table is block- or cohort-level aggregate |
| Vendor / source data files | none staged; referenced by filename only |

The documentation does report aggregate counts down to single digits in two places —
`harmonization_corrections.md` gives 12 and 4 matched cells for the two Olink blocks,
being the number of participant-visits assayed by both projects in a pair. These are
assay-overlap counts, not attributes of identifiable people, and they are load-bearing
for the correction's credibility, so they stay.

Note this means `EDA_PCA_plots/*-eigenvectors-*.tab` are excluded, which is intended —
they are participant-level (PATNO + PC scores). The dataset and dictionary themselves
are excluded too; they go through the usual access path, not a git push.
