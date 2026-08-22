# AD/PD abstract — submission draft

Drafted 2026-08-20 against the AD/PD formatting rules in `abstract_outline.txt`;
RESULTS restructured 2026-08-21 onto the four questions the stakeholders asked.

Every result-level number comes from the newest `meta/META_*.csv` per run — the phase 1
wrap meta-analyses — read by `abstract_results_buckets.py`; concordance comes from
`concordance.csv`. Nothing is typed in by hand: `python3 validate_abstract_numbers.py`
re-derives every figure in this document and fails loudly on a mismatch (65/65 pass,
including a row-by-row check that the browser app's snapshot carries the identical
robust set).

**Rules being met**

| Requirement | Limit | This draft |
|---|---|---|
| Title, ALL CAPS | ≤ 20 words | **16** |
| Structured sections | Objectives / Methods / Results / Conclusions | all four present |
| Body word count | ≤ 280 words | **277** |
| Figure | JPG, < 500 KB, ≤ 600 × 800 px | `adpd_figure.pptx` — one slide, four panels, sized to export at exactly 600 × 800 |
| Tables/graphs in text | count toward the limit | none used, so none charged |
| Original data | unpublished | yes — 2026-08-20 batch |

Still to supply at submission: author full names, emails, institutional affiliations
(department, hospital/institution, city, country), and the Theme / Topic / Sub-Topic
assignment.

---

## Title

> **A HARMONIZED MULTI-PLATFORM PROTEOMIC RESOURCE FOR PARKINSON'S DISEASE WITH OPEN
> CODE, BROWSABLE RESULTS AND AI-ASSISTED INTERPRETATION**

16 words.

---

## Body

### OBJECTIVES

PPMI proteomics span six projects, two platforms and two biofluids, confining most
analyses to one panel. We harmonized them with GP2 genetics and asked the four
questions stakeholders posed, releasing data, code, a results browser and AI-assisted
interpretation.

### METHODS

4,788 participants, 19,450 visits, 35,566 variables. Six harmonized blocks pool project
pairs sharing platform, panel and biofluid (11,567 analytes; 2,906 participants),
retaining source columns. Fifty-two prespecified logistic, linear, Cox and mixed-effects
analyses ran in European and Ashkenazi Jewish strata, then meta-analysed.
Robust = Bonferroni-significant, replicated in both strata, direction-concordant.

### RESULTS

Of 711,927 associations tested, 296 were robust; 280 answered the four questions, across
28 analyses and 136 proteins.

Baseline differences by NSD status, genotype or NSD-ISS stage gave 88 robust results (53
proteins). DDC was raised wherever neuronal synuclein disease was (NSD versus HC P=2×10⁻⁴⁰;
CSF and plasma, both platforms) and rose with stage alongside CD276; within
GBA carriers, CSF SPAG7 and LIMD1 separated NSD-positive participants. Polygenic risk,
carrier contrasts excluded, rose with NSD status (P=6×10⁻¹³) and stage.

Genetic versus sporadic PD was thinnest — 3 proteins, led by lower CSF ENO2 in LRRK2-PD
(P=3×10⁻⁵), surviving NSD-status adjustment.

Thirty-five results predicted progression: low CSF NPTX2 and high phospho-tau preceded
MoCA decline, MoCA<26 and cognitive milestones; DDC predicted milestones and stage D. OFF
motor slope and DaT gave 0 and 1.

Longitudinally, 146 results across 9 contrasts included 75 within-subject rate
differences (32 proteins), DDC recurring in six, then FLT1, AOC3, NEFL, VEGFD.

### CONCLUSIONS

All four questions returned replicated protein signals and DDC answered every one, while
genetic-versus-sporadic contrasts stayed thin. Genomic inflation was controlled (median
λ=1.10). Harmonized data, code, browser and AI interpretation are openly available.

**Word budget** — Objectives 39 · Methods 50 · Results 154 · Conclusions 34 · **total 277 / 280**. Section headings are counted.

---

## Results by stakeholder question

The batch is organised by model — logistic, linear, Cox, mixed-effects. The stakeholders
asked four questions that cut across those models, so `abstract_results_buckets.py`
assigns every one of the 52 analyses to exactly one bucket, *including* the ones that fall
outside the four. The assignment asserts exhaustiveness: an analysis added to `batch.yaml`
without a bucket fails the build rather than quietly vanishing from the abstract.

| Question | Analyses | With a robust result | Robust | Proteins | Named in the abstract |
|---|---:|---:|---:|---:|---|
| 1 · Baseline: NSD status, genotype, NSD-ISS stage | 12 | 8 | 88 | 53 | DDC, CD276, SPAG7, LIMD1 |
| 2 · Genetic PD vs sporadic PD | 4 | 4 | 11 | 3 | ENO2 |
| 3 · Baseline predictors of progression | 12 | 7 | 35 | 20 | NPTX2, phospho-tau, DDC |
| 4 · Longitudinal change | 10 | 9 | 146 | 72 | DDC, FLT1, AOC3, NEFL, VEGFD |
| *outside the four* | 14 | 5 | 16 | 10 | *not quoted* |

**280 of the 296 robust results answer one of the four**, across 28 analyses and 136
proteins. Full per-question hit tables are in `abstract_results_buckets.md`; every robust
row with its bucket, both ancestry betas and I² is in `abstract_results_buckets.csv`.

**What a count is.** Rows are analyte-level, so one protein measured in CSF and plasma on
two platforms can contribute four robust results to a single analysis — which is the
replication that makes it credible, not double-counting. Protein counts are therefore
always quoted alongside.

**Question 1.** DDC is robust in six of the eight baseline analyses that returned anything
at all, in both fluids and on both platforms. CD276 rises with stage on both platforms in
three stage analyses. SPAG7 (P=1.3×10⁻⁵⁶) and LIMD1 separate NSD-positive from
NSD-negative GBA carriers. Four analyses returned nothing robust — sporadic PD's NSD
contrast, 2A vs 2B in either sample, and PD-only 2A vs 3.

**Question 2.** The thinnest by a wide margin. Of its 11 robust results only **4 are
proteins** (ENO2 twice, TLR3, DDC); the other 7 are PRS157, two BMP lipids and UPSIT. ENO2
is the one that carries: lower in LRRK2-PD than sporadic PD, and still there once NSD
status is in the model.

**Question 3.** The negatives are as informative as the hits. OFF-state motor slope
returned nothing robust with or without LEDD adjustment, and neither did conversion from
stage 2A or 2B, nor cognitive-state worsening. The DaT slope returned exactly one (IL4);
`lowput_ratio` also changed definition in this build (see the header of
`build_step_5_derive.py`), so any DaT result should be read against that note.

**Question 4.** The 146 splits into **75 within-subject rate terms** (32 proteins) and 71
between-subject terms. Only the within term is a protein moving faster in one group; the
between term carries who was observed when. The abstract quotes 75 explicitly for that
reason, and the figure names proteins from the within term only.

## Polygenic risk — the associations that are NOT carrier-driven

`PRS157` contains LRRK2 and GBA variants, so it partly proxies carrier status. Any PRS
result read off a *carrier contrast* is therefore circular. Filtering the meta-analysis to
`robust_FE` (replicated in both ancestry strata, family-wise Bonferroni significant,
direction-concordant, not stratum-dominated) gives **8** PRS associations, of which **2**
sit in carrier contrasts and are dropped. The remaining **6** are below — all are PRS main
effects on disease status, stage or progression, none is a carrier contrast, and none is
an `_x_PRS` protein interaction.

| Analysis | β (FE) | P (FE) | N | β EUR | β AJ | I² |
|---|---|---|---|---|---|---|
| `NSD_vs_HC` | +0.67 | 5.9×10⁻¹³ | 2,312 | +0.63 | +2.87 | 90.7 |
| `NSD_stage_early_vs_late` | +0.37 | 2.4×10⁻¹¹ | 2,007 | +0.37 | +0.31 | 0.0 |
| `NSD_stage_2A_vs_3` | +0.53 | 4.9×10⁻¹¹ | 1,235 | +0.53 | +0.63 | 0.0 |
| `NSD_stage_2B_vs_3` | +0.28 | 4.6×10⁻⁵ | 1,194 | +0.30 | +0.09 | 0.0 |
| `NSD_vs_notNSD_prodromal` | −0.25 | 1.5×10⁻⁴ | 1,805 | −0.29 | −0.02 | 47.9 |
| `cox_stage_d` | +0.23 | 4.9×10⁻⁴ | 1,334 | +0.18 | +0.43 | 55.7 |

Written to `prs_robust_noncarrier.csv`.

**The one-line reading:** PD polygenic risk rises with NSD status and with NSD-ISS stage,
and predicts progression to stage D — replicated across both ancestries.

**Residual caveat, worth stating if asked.** Excluding carrier *contrasts* removes the
circularity in the comparison, but LRRK2 and GBA carriers are still present in these
general-cohort samples and still contribute to PRS157. A fully clean version would refit
with carriers excluded from the sample; that has not been run. `NSD_vs_HC` also carries
I²=90.7, so its pooled β is heterogeneous even though both strata are individually
significant and agree in direction.

## Figure 1 — one figure, four panels

Built: **`adpd_figure.pptx`** (33 KB) — one slide, 6.25 × 8.33 in, which exports to exactly
600 × 800 px at 96 DPI, the AD/PD image cap. Native PowerPoint shapes throughout, so every
box and label stays editable; rebuild with `python3 build_adpd_figure.py`. Panels stack top
to bottom, and **panel B is now the centre of the figure**: one band per stakeholder
question, each naming the proteins that question returned.

Every count and every protein name is read at build time from `meta/META_*.csv` through
`abstract_results_buckets.py` — the same source as the abstract text — so the two cannot
drift. Nothing about the results is typed into the figure script.

The figure is pitched at a non-specialist: no axis to decode, one line of plain language
per question, and colour used for exactly one thing throughout — **teal = CSF, amber =
plasma, dark = a protein robust in both**. The build prints each question's hits and the
lowest shape edge (7.65 in of 8.33) so overflow is caught without opening PowerPoint.

### Panel A — how the resource is built

Four stacked steps, each a full-width band with an arrow to the next:

1. **Six assay projects** — Olink and NULISA, each run in CSF or plasma or both.
   *34,902 analytes, held apart.*
2. **Pooled in pairs** — projects sharing a platform, panel and biofluid go onto one
   scale. Olink is never mixed with NULISA.
   *11,567 analytes · 2,906 participants · original columns kept.*
3. **52 prespecified analyses** — 18 logistic, 16 survival, 10 trajectory, 8 linear, each
   run in the European and Ashkenazi Jewish strata and then combined.
4. **296 robust associations** — replicated in both ancestries and concordant in
   direction; 296 of 711,927 tested, across 33 of the 52 analyses.

### Panel B — the four questions, and the proteins each returned

The panel opens with the definition in words a non-specialist can hold — **robust =
significant after correction + found in both ancestry groups + pointing the same way in
each** — and the line that 280 of the 296 robust results answer these four questions.

Then one band per question. Each band carries the question in plain words, a bar scaled
across the four so 11 next to 146 is visible without reading, the counts, and a full-width
row of named protein chips coloured by fluid:

| # | Question | Robust | Proteins | Named on the slide |
|---|---|---:|---:|---|
| 1 | Who differs at baseline | 88 | 53 | DDC 16 · CD276 5 · pTau-181 3 · ITGAV 2 · CA3 2 · SPAG7 |
| 2 | Genetic PD vs sporadic PD | 11 | 3 | ENO2 2 · TLR3 · DDC |
| 3 | What predicts later progression | 35 | 20 | DDC 3 · NPTX2 3 · BD-pTau-231 2 · IL23 2 · pTau-181 2 · pTau-217 2 |
| 4 | What changes over time | 146 | 72 | DDC 15 · FLT1 5 · AOC3 4 · NEFL 4 · VEGFD 4 · ABHD14B 3 |

The number on a chip is how many robust results that protein carries in that question —
so DDC 16 means sixteen, across fluids, platforms and analyses, not one finding repeated.

Two selection rules, both deliberate. Chips rank by recurrence first, because replication
across blocks is what `robust` is buying, **but the single strongest result is forced onto
the slide if that ranking drops it** — which is what keeps SPAG7 (P=1.3×10⁻⁵⁶, seen once)
beside DDC (seen sixteen times). And for question 4 the chips are drawn from the
within-subject term only, because that is the term that means *changing*; the band's count
still covers everything those analyses returned, with "75 are rate-of-change" printed
next to it.

A single muted line closes the panel: the genetic risk score PRS157 gives 6 robust
results, rising with NSD status (β=0.67) and with NSD-ISS stage and predicting stage D,
with carrier contrasts excluded because PRS157 contains LRRK2 and GBA variants.

### Panel C — one protein answers all four

**DDC, 51 robust results.** Dopa decarboxylase is the enzyme that makes dopamine, so it
needs no explaining to this audience: raised where disease is at baseline, higher at later
stages, predicting milestones and stage D, and rising faster afterwards — spinal fluid and
blood, both technologies, both ancestry groups.

Underneath, its count under each of the four questions — 16 · 1 · 3 · 31 — so the panel
ties straight back to panel B rather than restating it.

### Panel D — open at every stage

A single band, four items, no arrows — a statement of availability, not a process:
controlled-access harmonized data · public codebase · results browser covering all 52
analyses rather than only the hits · AI interpretation in plain language.

Cross-platform concordance stays cut from the figure — it is a methods point, and the
figure has to land for a non-specialist in one pass. It is kept in the appendix below.

## Appendix — cross-platform concordance (cut from the figure, kept for questions)

Pairing one participant at one visit; analytes at P>0.05 are reported NS and excluded from
the median, with their R² still recorded in `concordance.csv`.

| Comparison | Analytes | NS (P>0.05) | Median R² among significant |
|---|---:|---:|---:|
| Olink vs NULISA, plasma | 319 | 21 (6.6%) | **0.54** |
| Olink vs NULISA, CSF | 320 | 47 (14.7%) | **0.14** |
| CSF vs plasma, same platform | 5,772 | 4,168 (72%) | 0.009 |

Two technologies agree on about half the variance in plasma, far less in CSF, and the two
biofluids share almost nothing. Note the sample sizes: with 2,000–3,600 matched
participant-visits, nominal P<0.05 is cleared by R²≈0.001, so the NS counts and the
Bonferroni column in `concordance.csv` carry more weight than the nominal threshold.
Reproduce with `python3 concordance.py`.

---

## Open items before submission

1. **λ figures are parsed from the meta-analysis log of the run that produced `meta/`**
   (`meta_analysis_20260820_120837.log`, 62 run × term cells): median 1.101, IQR
   1.02–1.19, full range 0.49–1.90, 76% within 1.2 — the same 62 cells the app's
   reliability table carries. By analysis type — baseline 1.19 (0.99–1.37), time-to-event
   1.14 (0.92–1.41), trajectory 1.04 (0.69–1.90), progression slope 1.01 (0.49–1.12). The
   separate permutation sweep is still running and is **not** needed for this abstract; it
   feeds the app's false-positive-rate tile instead.

2. **PRS157 proxies carrier status** — it contains LRRK2 and GBA variants. The abstract
   and figure now quote it in one clause each, and both say *carrier contrasts excluded*:
   the 6 reportable results are PRS main effects on status, stage and stage-D progression,
   never a carrier contrast and never an `_x_PRS` protein interaction. The residual caveat
   still stands — carriers remain inside the general-cohort samples — and the 15 robust
   `_x_PRS` results are excluded from every count in this document.

3. **The negative-control check is inert under harmonized scope.** Only 2 of 10 spike-in
   controls survive into the harmonized blocks and `control_report` requires 4 per cell,
   so it reports nothing — while the project README still describes it as a live readout.
   Either lower the floor or reword the claim.

4. **The figure needs exporting to JPG.** The PowerPoint is built and sized correctly;
   PowerPoint's *File → Export → JPEG* at 96 DPI gives the 600 × 800 the submission
   wants, well under the 500 KB cap. **Panel B has not been eyeballed at final size** —
   this box has no PowerPoint or renderer, so the layout was checked geometrically
   (no shape out of bounds, no chip overlap, lowest edge 7.65 in of 8.33). Six chips per
   band at 5.8 pt is the tightest thing on the slide; check it survives the export.

5. **Two source scripts are new** and belong in the phase 1 wrap commit:
   `abstract_results_buckets.py` (assigns all 52 analyses to the four questions, writes
   the per-question tables) and `validate_abstract_numbers.py` (re-derives every number
   in this document from `meta/`, and checks the app snapshot matches row for row).
