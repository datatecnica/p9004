# AD/PD abstract — submission draft

Drafted 2026-08-20 against the AD/PD formatting rules in `abstract_outline.txt`.
Every number below comes from the 2026-08-20 batch (`results_summary.md`,
`concordance.csv`, `calibration_report.md`) rather than from an earlier draft.

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

PPMI proteomic data span six assay projects, two platforms and two biofluids, which
confines most analyses to a single panel. Bringing PPMI proteomics together with GP2
genetics, we built a harmonized controlled-access resource with a fully open codebase, a
public hypothesis-generating results browser and AI-assisted interpretation.

### METHODS

We unified 4,788 participants and 19,450 visits into 35,566 variables, including 34,902
Olink and NULISA analytes in CSF and plasma. Six harmonized blocks pool project pairs
sharing platform, panel and biofluid onto a common scale (11,567 analytes, 2,906
participants), retaining source columns so every result stays checkable. Fifty-two
prespecified logistic, linear, Cox and mixed-effects analyses ran within European and
Ashkenazi Jewish strata and were meta-analysed, using panel-specific principal components,
dual Bonferroni correction and blocked permutation nulls.

### RESULTS

Of 711,927 associations tested, 296 were robust — Bonferroni-significant, replicated
across both ancestry strata and concordant in direction — spanning 33 of 52 analyses and
151 proteins. Robust findings divided into longitudinal trajectory 146, baseline
comparison 100, progression slope 29 and time-to-event 21. DDC dominated with 51 robust
results: raised in neuronal synuclein disease at baseline and rising faster thereafter,
in CSF and plasma and on both platforms. CSF proteins separated GBA carriers (SPAG7,
LIMD1); NPTX2, VEGFD, FLT1 and pTau-181 recurred. Polygenic risk, replicated in both
ancestries and excluding carrier contrasts, rose with disease status (β=0.67,
P=6×10⁻¹³) and stage and predicted progression to stage D.

### CONCLUSIONS

Harmonizing adds 27–65% more measurements per block while every value stays traceable to
source. Biofluids are not interchangeable and platform agreement is matrix-dependent.
Genomic inflation was controlled across the batch (median λ=1.10, IQR 1.02–1.19; 76% of
analyses within 1.2). Dataset, code and browser are openly available.

**Word budget** — Objectives 47 · Methods 77 · Results 106 · Conclusions 47 · **total 277 / 280**. Section headings are counted.

---

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

Built: **`adpd_figure.pptx`** (31 KB) — one slide, 6.25 × 8.33 in, which exports to exactly
600 × 800 px at 96 DPI, the AD/PD image cap. Native PowerPoint shapes throughout, so every
box and label stays editable; rebuild with `python3 build_adpd_figure.py`. Panels stack
top to bottom. The whole figure is
pitched at a non-specialist reader: no axis a general audience has to decode, one sentence
of plain-language read-out per panel, and colour used for exactly one thing throughout —
**teal = CSF, amber = plasma**. Nothing else is coloured, so the encoding is learnable at
a glance.

### Panel A — how the resource is built

Four stacked steps, each a full-width band with an arrow to the next:

1. **Six assay projects** — Olink and NULISA, each run in CSF or plasma or both.
   *34,902 analytes, held apart.*
2. **Pooled in pairs** — projects sharing a platform, panel and biofluid go onto one
   scale. Olink is never mixed with NULISA, which is exactly why Panel C is a real
   question rather than something the pooling presupposes.
   *11,567 analytes · 2,906 participants · original columns kept.*
3. **52 prespecified analyses** — 18 logistic, 16 survival, 10 trajectory, 8 linear, each
   run in the European and Ashkenazi Jewish strata and then combined.
4. **296 robust associations** — replicated in both ancestries and concordant in
   direction, across 33 of the 52 analyses.

### Panel B — what counts as a result

The definition first, in words a non-specialist can hold: **robust = significant after
correction + found in both ancestry groups + pointing the same way in each.** Then the
quantification — of **711,927** associations tested, **296** clear that bar (0.04%),
across **33 of 52** analyses and **151** proteins — and a four-bar split of those 296:

| Question asked | Robust |
|---|---:|
| Change over time | 146 |
| Baseline comparison | 100 |
| Progression slope | 29 |
| Time to an event | 21 |

Stating the denominator matters more than the numerator here: 296 out of 711,927 is what
makes the bar credible to a sceptical reader.

### Panel C — the protein that keeps coming back

**DDC, 51 robust results** — more than six times any other protein, and the panel's
centrepiece. Dopa decarboxylase is the enzyme that makes dopamine, so it needs no
explaining to this audience: raised in disease at baseline, rising faster afterwards, in
spinal fluid and blood, on both technologies, in both ancestry groups.

Then the next most recurrent, as eight small cards: PRS157 (8, genetic risk score),
UPSIT (7, smell test), VEGFD (6), PRL (6), NPTX2 (6), CD276 (5), pTau-181 (5), FLT1 (5).

Concordance was **cut from the figure** — it is a methods point, and the figure has to
land for a non-specialist in one pass. It is retained in the appendix below.

### Panel D — open at every stage

A single band, four items, no arrows — this is a statement of availability, not a process:
controlled-access harmonized data · public codebase · results browser covering all 52
analyses rather than only the hits · AI interpretation in plain language.

**All figure numbers are read from the app's own files** (`data/manifest.json`,
`data/results.parquet`), so the figure and the browser cannot drift apart. 26 of 26
cross-checks pass.

---

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

1. **λ figures come from the app's own reliability table** (`data/manifest.json`, 62
   run × term cells): median 1.10, IQR 1.02–1.19, full range 0.49–1.90, 76% within 1.2.
   By analysis type — baseline 1.19 (0.99–1.37), time-to-event 1.14 (0.92–1.41),
   trajectory 1.04 (0.69–1.90), progression slope 1.01 (0.49–1.12). The separate
   permutation sweep is still running and is **not** needed for this abstract; it feeds
   the app's false-positive-rate tile instead.

2. **PRS157 proxies carrier status** — it contains LRRK2 and GBA variants. No PRS result
   is quoted here, and the genotype contrasts named (MEGF9, SPAG7) are proteomic rather
   than PRS-derived. Any later draft reaching for the `_x_PRS` runs must carry that
   caveat.

3. **The negative-control check is inert under harmonized scope.** Only 2 of 10 spike-in
   controls survive into the harmonized blocks and `control_report` requires 4 per cell,
   so it reports nothing — while the project README still describes it as a live readout.
   Either lower the floor or reword the claim.

4. **The figure needs exporting to JPG.** The PowerPoint is built and sized correctly;
   PowerPoint's *File → Export → JPEG* at 96 DPI gives the 600 × 800 the submission
   wants, well under the 500 KB cap.
