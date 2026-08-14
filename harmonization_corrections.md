# Paired-project comparison and proposed corrections

Per-analyte table: `harmonization_corrections.tab` (11,567 rows). The older project is the reference; the newer is mapped onto it via `B_corrected = slope * B + intercept`, fitted on (EVENT_ID x COHORT) matched cells so cohort composition cannot drive the estimate.

## Gross comparison

| block | projects | core | matched cells | overlapping visits | raw median gap | fitted slope (median) | slope IQR | % analytes slope-corrected | intercept (median) |
|---|---|---|---|---|---|---|---|---|---|
| `olink_plasma` | p293_olink_plasma -> p314_Plasma | 5,400 | 12 | 56 | +0.326 | **1.247** | [1.10, 1.50] | 77% | +0.312 |
| `olink_csf` | p277_CSF -> p314_CSF | 5,416 | 4 | 30 | +0.032 | **0.945** | [0.84, 1.03] | 50% | +0.021 |
| `nulisa_cns_plasma` | p288_CNS_plasma -> p312_Neuro_Plasma | 127 | 10 | 0 | -0.057 | **0.967** | [0.85, 1.09] | 57% | +0.202 |
| `nulisa_cns_csf` | p282_CNS_CSF -> p312_Neuro_CSF | 128 | 6 | 0 | +0.017 | **0.859** | [0.77, 1.00] | 80% | +1.610 |
| `nulisa_inf_plasma` | p288_Inflammation_plasma -> p312_Inflammation_Plasma | 249 | 10 | 0 | -0.050 | **0.978** | [0.93, 1.02] | 20% | +0.212 |
| `nulisa_inf_csf` | p282_Inflammation_CSF -> p312_Inflammation_CSF | 247 | 6 | 0 | +0.194 | **0.891** | [0.77, 1.02] | 71% | +1.094 |

## Does the correction hold out?

Fitted on half the matched cells, residual measured on the other half. `|gap| before` is the same quantity with no correction applied.

| block | median \|gap\| before | median \|residual\| after | reduction |
|---|---|---|---|
| `olink_plasma` | 0.567 | **0.032** | 94% |
| `olink_csf` | 0.050 | **0.011** | 79% |
| `nulisa_cns_plasma` | 0.370 | **0.020** | 95% |
| `nulisa_cns_csf` | 0.486 | **0.034** | 93% |
| `nulisa_inf_plasma` | 0.063 | **0.023** | 64% |
| `nulisa_inf_csf` | 0.300 | **0.037** | 88% |

## Per-analyte distributions

| block | slope <0.9 | slope 0.9-1.1 | slope >1.1 | \|intercept\| >0.5 | analytes with overlap r2 > 0.5 |
|---|---|---|---|---|---|
| `olink_plasma` | 121 | 1,261 | 4,018 | 2,923 | 579 / 1,162 |
| `olink_csf` | 2,033 | 2,714 | 669 | 45 | 26 / 1,335 |
| `nulisa_cns_plasma` | 42 | 55 | 30 | 91 | — no usable overlaps |
| `nulisa_cns_csf` | 81 | 25 | 22 | 111 | — no usable overlaps |
| `nulisa_inf_plasma` | 30 | 199 | 20 | 53 | — no usable overlaps |
| `nulisa_inf_csf` | 130 | 72 | 45 | 182 | — no usable overlaps |

Columns in `harmonization_corrections.tab`: block, platform, fluid, analyte, col_reference, col_mapped, n_ref, n_map, median_ref, median_map, iqr_ref, iqr_map, raw_median_gap, fitted_slope, fitted_intercept, applied_slope, applied_intercept, heldout_residual, n_matched_cells, n_overlap, overlap_r, overlap_r2, overlap_deming, overlap_sd_ref, overlap_sd_map

