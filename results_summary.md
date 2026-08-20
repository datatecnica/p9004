# Proteomics Data Mine — Results Summary

Generated: 2026-08-20 12:20

Inputs: `results/` (1351623 per-stratum rows across 97 CSVs), `meta/` (717596 META rows across 52 CSVs).

Total runs surveyed: **52**.  Bonferroni-significant unique hits: **698** `(run, predictor, outcome)` tuples across **45** runs (deduplicated from 1247 per-source rows).

Two Bonferroni lenses are reported. **Run-wide** divides 0.05 by every predictor fitted in the run (374 hits); **family** divides by the predictors in the same assay panel (698 hits), which keeps each panel comparable with earlier cuts now that Projects 312/314 have roughly doubled the run-wide denominator. 374 hits pass both. The `lens` column marks which.


## Coverage table

`λ` (genomic inflation) is recomputed here from the P column of each CSV (validates the regression-script log values). LMM runs decompose time into a **within**-subject (trajectory) and a **between**-subject (cross-sectional) term, reported as `within=… / between=…`; the two are separate questions and are calibrated separately.

| Run | Model | EUR N | AJ N | Meta? | λ EUR | λ AJ | λ FE | λ RE |
|---|---|---|---|---|---|---|---|---|
| `CI_MOCA_baseline` | LOGIT | [68, 987] | [64, 216] | ✓ | 1.019 | 1.160 | 1.029 | 0.858 |
| `CI_PI_baseline` | LOGIT | [51, 972] | [61, 215] | ✓ | 1.118 | 1.268 | 1.179 | 0.935 |
| `NSD_stage_2A_vs_2B` | LOGIT | [45, 992] | [25, 99] | ✓ | 1.025 | 1.404 | 1.023 | 0.870 |
| `NSD_stage_2A_vs_2B_PDonly` | LOGIT | [69, 180] | — | ✓ | 1.254 | — | 1.254 | 1.254 |
| `NSD_stage_2A_vs_3` | LOGIT | [127, 1,101] | [30, 134] | ✓ | 1.263 | 1.259 | 1.260 | 0.977 |
| `NSD_stage_2A_vs_3_PDonly` | LOGIT | [125, 519] | — | ✓ | 1.089 | — | 1.089 | 1.089 |
| `NSD_stage_2B_vs_3` | LOGIT | [37, 1,073] | [35, 121] | ✓ | 1.309 | 1.257 | 1.284 | 0.955 |
| `NSD_stage_2B_vs_3_PDonly` | LOGIT | [37, 659] | [34, 96] | ✓ | 1.210 | 1.415 | 1.211 | 0.940 |
| `NSD_stage_early_vs_late` | LOGIT | [45, 1,789] | [25, 218] | ✓ | 1.395 | 1.235 | 1.369 | 0.991 |
| `NSD_vs_HC` | LOGIT | [68, 2,051] | [20, 261] | ✓ | 1.259 | 1.275 | 1.282 | 1.035 |
| `NSD_vs_notNSD_GBA` | LOGIT | [29, 66] | [84, 178] | ✓ | 1.866 | 1.245 | 1.226 | 0.961 |
| `NSD_vs_notNSD_LRRK2` | LOGIT | [39, 107] | [74, 232] | ✓ | 1.389 | 1.229 | 1.200 | 0.930 |
| `NSD_vs_notNSD_prodromal` | LOGIT | [83, 1,418] | [44, 387] | ✓ | 1.412 | 1.080 | 1.373 | 1.007 |
| `NSD_vs_notNSD_sPD` | LOGIT | [150, 792] | — | ✓ | 1.146 | — | 1.146 | 1.146 |
| `cox_cogstate_worsen` | COX | [108, 590] | [43, 158] | ✓ | 1.159 | 1.416 | 1.177 | 0.882 |
| `cox_cogstate_worsen_x_PRS` | COX | [108, 576] | [63, 157] | ✓ | 0.998 | 1.654 | 1.219 | 0.776 |
| `cox_moca_lt26` | COX | [46, 703] | [51, 133] | ✓ | 1.088 | 1.382 | 1.144 | 0.937 |
| `cox_moca_lt26_x_PRS` | COX | [95, 692] | [51, 132] | ✓ | 0.970 | 1.240 | 0.983 | 0.782 |
| `cox_nsd_2a_to_later` | COX | [82, 392] | — | ✓ | 1.097 | — | 1.097 | 1.097 |
| `cox_nsd_2a_to_later_x_PRS` | COX | [120, 390] | — | ✓ | 1.409 | — | 1.409 | 1.409 |
| `cox_nsd_2b_to_later` | COX | [35, 409] | — | ✓ | 1.061 | — | 1.061 | 1.061 |
| `cox_nsd_2b_to_later_x_PRS` | COX | [35, 407] | — | ✓ | 1.090 | — | 1.090 | 1.090 |
| `cox_pm_any` | COX | [51, 826] | [52, 155] | ✓ | 1.092 | 1.284 | 1.169 | 0.960 |
| `cox_pm_any_x_PRS` | COX | [51, 811] | [52, 154] | ✓ | 0.983 | 1.119 | 0.939 | 0.731 |
| `cox_pm_cog_any` | COX | [128, 891] | [62, 185] | ✓ | 1.073 | 1.265 | 1.168 | 0.935 |
| `cox_pm_cog_any_x_PRS` | COX | [128, 875] | [62, 184] | ✓ | 1.171 | 1.139 | 1.317 | 0.988 |
| `cox_pm_mc_any` | COX | [61, 918] | [61, 187] | ✓ | 1.103 | 1.217 | 1.168 | 0.952 |
| `cox_pm_mc_any_x_PRS` | COX | [124, 900] | [61, 186] | ✓ | 1.134 | 1.362 | 1.130 | 0.854 |
| `cox_stage_d` | COX | [224, 1,001] | [174, 333] | ✓ | 1.083 | 1.099 | 1.143 | 0.897 |
| `cox_stage_d_x_PRS` | COX | [228, 996] | [179, 331] | ✓ | 1.563 | 1.091 | 0.918 | 0.611 |
| `sPD_vs_GBA` | LOGIT | [168, 876] | [45, 107] | ✓ | 0.974 | 1.520 | 1.016 | 0.900 |
| `sPD_vs_GBA_SAAadj` | LOGIT | [164, 826] | [41, 86] | ✓ | 0.956 | 1.575 | 0.986 | 0.896 |
| `sPD_vs_LRRK2` | LOGIT | [171, 886] | [51, 147] | ✓ | 1.107 | 1.622 | 1.065 | 0.874 |
| `sPD_vs_LRRK2_SAAadj` | LOGIT | [167, 844] | [50, 132] | ✓ | 1.191 | 1.815 | 1.103 | 0.886 |
| `slope_lowput_ratio` | OLS | [58, 851] | [47, 143] | ✓ | 0.925 | 0.962 | 0.962 | 0.830 |
| `slope_lowput_ratio_x_PRS` | OLS | [58, 835] | [47, 142] | ✓ | 0.799 | 0.956 | 0.798 | 0.612 |
| `slope_moca` | OLS | [65, 932] | [23, 196] | ✓ | 1.021 | 1.045 | 1.050 | 0.876 |
| `slope_moca_x_PRS` | OLS | [65, 914] | [23, 195] | ✓ | 0.603 | 0.790 | 0.493 | 0.416 |
| `slope_updrs3_off` | OLS | [62, 863] | [20, 156] | ✓ | 0.972 | 0.986 | 1.009 | 0.842 |
| `slope_updrs3_off_LEDD_adj` | OLS | [62, 863] | [20, 156] | ✓ | 0.979 | 0.981 | 1.016 | 0.834 |
| `slope_updrs3_off_x_PRS` | OLS | [62, 846] | [20, 155] | ✓ | 1.430 | 1.454 | 1.120 | 0.753 |
| `slope_updrs3_off_x_PRS_LEDD_adj` | OLS | [62, 846] | [47, 155] | ✓ | 1.405 | 1.447 | 1.099 | 0.742 |
| `trajectory_CI_MOCA` | LMM_RS_lbfgs | [404, 2,670] | [114, 1,053] | ✓ | within=0.875 / between=1.288 | within=0.946 / between=1.000 | within=0.999 / between=1.313 | within=0.815 / between=1.053 |
| `trajectory_CI_PI` | LMM_RS_lbfgs | [376, 2,339] | [113, 1,028] | ✓ | within=0.916 / between=1.060 | within=0.940 / between=1.178 | within=0.950 / between=1.122 | within=0.758 / between=0.962 |
| `trajectory_HCNSDneg_vs_PDNSDpos` | LMM_RS_lbfgs (+1 variants) | [233, 1,529] | [124, 269] | ✓ | within=1.032 / between=0.981 | within=0.758 / between=0.959 | within=1.025 / between=1.020 | within=0.798 / between=0.850 |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | LMM_RS_lbfgs | [85, 776] | [81, 113] | ✓ | within=1.090 / between=1.405 | within=1.040 / between=1.054 | within=0.981 / between=1.142 | within=0.821 / between=0.746 |
| `trajectory_HC_vs_PD` | LMM_RS_lbfgs | [271, 1,759] | [158, 374] | ✓ | within=1.013 / between=1.019 | within=0.744 / between=0.911 | within=1.029 / between=1.012 | within=0.912 / between=0.842 |
| `trajectory_HC_vs_Prodromal` | LMM_RS_lbfgs | [170, 1,258] | [84, 718] | ✓ | within=1.199 / between=1.650 | within=0.682 / between=0.925 | within=1.141 / between=1.505 | within=0.874 / between=1.023 |
| `trajectory_NSD_vs_HC` | LMM_RS_lbfgs | [289, 2,009] | [176, 346] | ✓ | within=0.989 / between=1.042 | within=0.639 / between=0.897 | within=1.003 / between=1.060 | within=0.857 / between=0.853 |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | LMM_RS_lbfgs (+1 variants) | [254, 1,705] | [144, 314] | ✓ | within=1.133 / between=1.542 | within=1.267 / between=1.098 | within=0.944 / between=1.328 | within=0.684 / between=0.968 |
| `trajectory_Prodromal_vs_PD` | LMM_RS_lbfgs | [369, 2,353] | [111, 1,024] | ✓ | within=1.239 / between=1.958 | within=1.059 / between=1.014 | within=1.132 / between=1.898 | within=0.719 / between=1.206 |
| `trajectory_stage_2A_vs_2B` | LMM_RS_lbfgs | [108, 696] | [60, 105] | ✓ | within=0.825 / between=1.234 | within=1.436 / between=1.126 | within=0.700 / between=1.083 | within=0.563 / between=0.732 |

## λ by assay panel

Genomic inflation per **run × term × assay panel**, computed from the meta fixed-effect P where a meta exists and from EUR otherwise. A run-level λ averages panels with very different repeat-measure depth and platform behaviour, so the panel actually driving an inflated run is only visible here. Cells are blank where a panel has <20 tests in that run. **Bold** marks λ > 1.5.

Read small panels with the right error bar: λ estimated from *n* tests has a sampling SE of roughly `2.33/√n` under the null, so the ~130–250-test panels (p282, p288, p312) carry SE ≈ 0.15–0.21 against ≈ 0.03 for the ~5,400-test panels (p277, p293, p314). That is wide, but not wide enough to explain a λ of 2+ — those sit 5–8 SE above 1. Note also that λ > 1 is *not* by itself evidence of miscalibration: widespread true signal inflates λ too, and it grows with power. Use the negative controls and the permutation null in `calibration_check.py` to tell the two apart.

| Run | Term | Source | harmonized nulisa cns csf | harmonized nulisa cns plasma | harmonized nulisa inf csf | harmonized nulisa inf plasma | harmonized olink csf | harmonized olink plasma |
|---|---|---|---|---|---|---|---|---|
| `CI_MOCA_baseline` | effect | META_FE | 1.20 | 0.87 | 0.95 | 1.22 | 1.01 | 1.06 |
| `CI_PI_baseline` | effect | META_FE | **1.71** | 0.72 | 1.18 | 1.40 | 1.21 | 1.14 |
| `NSD_stage_2A_vs_2B` | effect | META_FE | **1.93** | 1.13 | **1.51** | 1.11 | 1.01 | 1.00 |
| `NSD_stage_2A_vs_2B_PDonly` | effect | META_FE | **2.21** | 1.16 | 1.30 | 1.36 | 1.35 | 1.15 |
| `NSD_stage_2A_vs_3` | effect | META_FE | **1.82** | 1.16 | **1.69** | 1.48 | 1.25 | 1.23 |
| `NSD_stage_2A_vs_3_PDonly` | effect | META_FE | 1.48 | 0.64 | 0.88 | 0.78 | 1.17 | 1.05 |
| `NSD_stage_2B_vs_3` | effect | META_FE | 1.10 | 1.34 | 1.13 | 1.17 | 1.21 | 1.39 |
| `NSD_stage_2B_vs_3_PDonly` | effect | META_FE | **1.65** | 1.40 | 1.11 | 0.90 | 1.14 | 1.28 |
| `NSD_stage_early_vs_late` | effect | META_FE | 1.45 | **1.92** | 1.20 | 1.31 | 1.32 | 1.43 |
| `NSD_vs_HC` | effect | META_FE | **1.56** | **2.45** | **1.65** | **2.38** | 1.11 | 1.41 |
| `NSD_vs_notNSD_GBA` | effect | META_FE | **1.57** | **1.91** | 1.28 | 1.37 | 1.37 | 1.08 |
| `NSD_vs_notNSD_LRRK2` | effect | META_FE | **1.81** | 1.18 | 1.22 | 1.27 | 1.09 | 1.29 |
| `NSD_vs_notNSD_prodromal` | effect | META_FE | **1.77** | **1.79** | 1.32 | **2.09** | 1.37 | 1.34 |
| `NSD_vs_notNSD_sPD` | effect | META_FE | **2.12** | 1.28 | **1.63** | 1.11 | 1.09 | 1.17 |
| `cox_cogstate_worsen` | effect | META_FE | 1.37 | 1.03 | 0.93 | 1.02 | 1.13 | 1.24 |
| `cox_cogstate_worsen_x_PRS` | effect | META_FE | 1.20 | 0.86 | 1.27 | 0.78 | 1.43 | 1.06 |
| `cox_moca_lt26` | effect | META_FE | **1.84** | 1.01 | 1.33 | **1.64** | 1.14 | 1.11 |
| `cox_moca_lt26_x_PRS` | effect | META_FE | 1.04 | 1.15 | 0.98 | 0.99 | 1.43 | 0.65 |
| `cox_nsd_2a_to_later` | effect | META_FE | 1.44 | 1.11 | 1.09 | 0.95 | 1.12 | 1.07 |
| `cox_nsd_2a_to_later_x_PRS` | effect | META_FE | **2.84** | **2.34** | **2.56** | **1.71** | 0.95 | **1.96** |
| `cox_nsd_2b_to_later` | effect | META_FE | 1.39 | 1.23 | 1.18 | 1.34 | 1.00 | 1.08 |
| `cox_nsd_2b_to_later_x_PRS` | effect | META_FE | 1.07 | 0.58 | 0.96 | 0.51 | 1.09 | 1.14 |
| `cox_pm_any` | effect | META_FE | **1.78** | 0.87 | 1.19 | 1.17 | 1.15 | 1.18 |
| `cox_pm_any_x_PRS` | effect | META_FE | 1.01 | 1.18 | 0.91 | 1.26 | 1.05 | 0.83 |
| `cox_pm_cog_any` | effect | META_FE | **1.96** | 1.19 | 1.37 | 1.26 | 1.17 | 1.13 |
| `cox_pm_cog_any_x_PRS` | effect | META_FE | **2.10** | 1.16 | 1.02 | 1.02 | 1.38 | 1.29 |
| `cox_pm_mc_any` | effect | META_FE | **1.55** | 0.94 | **1.60** | 1.15 | 1.16 | 1.16 |
| `cox_pm_mc_any_x_PRS` | effect | META_FE | 1.48 | 1.20 | **1.61** | 1.17 | 0.90 | 1.36 |
| `cox_stage_d` | effect | META_FE | 0.94 | 1.17 | 1.06 | 1.17 | 1.14 | 1.15 |
| `cox_stage_d_x_PRS` | effect | META_FE | 0.81 | 1.23 | 0.95 | 0.70 | 0.69 | 1.25 |
| `sPD_vs_GBA` | effect | META_FE | 1.07 | **1.89** | 1.01 | 1.03 | 1.04 | 0.98 |
| `sPD_vs_GBA_SAAadj` | effect | META_FE | 1.15 | 1.18 | 1.23 | 1.00 | 1.04 | 0.93 |
| `sPD_vs_LRRK2` | effect | META_FE | 1.38 | **1.58** | 1.28 | **1.61** | 0.99 | 1.09 |
| `sPD_vs_LRRK2_SAAadj` | effect | META_FE | **1.57** | 1.24 | 1.34 | **1.69** | 1.07 | 1.09 |
| `slope_lowput_ratio` | effect | META_FE | 1.17 | 0.94 | 0.82 | 1.02 | 0.95 | 0.97 |
| `slope_lowput_ratio_x_PRS` | effect | META_FE | 0.46 | 0.72 | 0.75 | 0.69 | 0.77 | 0.84 |
| `slope_moca` | effect | META_FE | 1.45 | 1.42 | 1.47 | **1.54** | 1.04 | 1.03 |
| `slope_moca_x_PRS` | effect | META_FE | 0.34 | 0.39 | 0.45 | 0.62 | 0.46 | 0.53 |
| `slope_updrs3_off` | effect | META_FE | 1.03 | 0.71 | 0.97 | 0.83 | 1.03 | 1.01 |
| `slope_updrs3_off_LEDD_adj` | effect | META_FE | 1.03 | 0.64 | 0.93 | 0.77 | 1.03 | 1.02 |
| `slope_updrs3_off_x_PRS` | effect | META_FE | 1.34 | 0.85 | 0.81 | 1.07 | **1.63** | 0.77 |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | META_FE | 1.14 | 0.83 | 0.81 | 0.99 | **1.59** | 0.75 |
| `trajectory_CI_MOCA` | within | META_FE | **1.55** | 1.25 | 1.36 | **2.29** | 0.96 | 0.98 |
| `trajectory_CI_MOCA` | between | META_FE | 1.14 | 1.30 | 1.29 | 1.07 | **1.78** | 1.00 |
| `trajectory_CI_PI` | within | META_FE | **1.65** | 1.24 | 1.43 | 1.08 | 0.94 | 0.93 |
| `trajectory_CI_PI` | between | META_FE | 0.97 | 0.80 | 0.99 | 1.23 | 1.40 | 0.93 |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | META_FE | **4.21** | **1.81** | **2.40** | 1.40 | 0.98 | 0.97 |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | META_FE | **1.82** | 1.23 | 1.08 | 1.06 | 1.18 | 0.85 |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | META_FE | **2.11** | **1.87** | 1.40 | **1.83** | 0.75 | 1.22 |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | META_FE | **2.53** | 1.22 | **1.53** | 1.18 | 1.11 | 1.12 |
| `trajectory_HC_vs_PD` | within | META_FE | **5.21** | **2.03** | **2.65** | **1.53** | 0.92 | 1.05 |
| `trajectory_HC_vs_PD` | between | META_FE | **2.12** | 1.28 | 1.14 | 1.13 | 1.16 | 0.85 |
| `trajectory_HC_vs_Prodromal` | within | META_FE | 1.27 | 1.10 | 0.93 | 1.17 | 0.99 | 1.33 |
| `trajectory_HC_vs_Prodromal` | between | META_FE | 0.93 | 1.27 | 1.30 | 1.30 | **1.71** | 1.35 |
| `trajectory_NSD_vs_HC` | within | META_FE | **4.75** | **1.66** | **2.53** | 1.19 | 0.92 | 1.02 |
| `trajectory_NSD_vs_HC` | between | META_FE | **2.57** | 0.99 | 1.26 | 0.95 | 1.15 | 0.91 |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | META_FE | 0.62 | **2.22** | 1.00 | **3.05** | 0.86 | 0.98 |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | META_FE | 1.04 | **2.30** | 1.17 | **2.19** | 1.24 | 1.40 |
| `trajectory_Prodromal_vs_PD` | within | META_FE | **1.55** | **1.91** | 0.78 | **1.51** | 0.99 | 1.27 |
| `trajectory_Prodromal_vs_PD` | between | META_FE | 1.23 | **3.06** | 1.18 | **1.86** | **1.88** | **1.93** |
| `trajectory_stage_2A_vs_2B` | within | META_FE | 0.79 | 0.82 | 1.08 | 0.97 | 0.71 | 0.67 |
| `trajectory_stage_2A_vs_2B` | between | META_FE | **1.83** | 1.27 | 1.09 | 1.40 | 1.06 | 1.08 |

**76 of 372** run × term × panel cells exceed λ = 1.5.


## Calibration — negative controls

The panels ship spike-in / assay control analytes (`*_CTRL_*`, `*_mCherry_*`) that cannot carry biological signal. They are **deliberately retained** in all outputs and flagged via `is_control`, because their rejection rate at α=0.05 is the empirical calibration of each run. A well-calibrated run sits near **5%**; anything above ~10% means the reported p-values are not trustworthy for that run, regardless of how small they are. Numerically-degenerate fits (`degenerate`, |β| < 1e-10) are excluded from this rate since they are a separate, separately-flagged failure mode.

| Run | control type-I EUR | control type-I AJ | status |
|---|---|---|---|
| `CI_MOCA_baseline` | — | — | ok |
| `CI_PI_baseline` | — | — | ok |
| `NSD_stage_2A_vs_2B` | — | — | ok |
| `NSD_stage_2A_vs_2B_PDonly` | — | — | ok |
| `NSD_stage_2A_vs_3` | — | — | ok |
| `NSD_stage_2A_vs_3_PDonly` | — | — | ok |
| `NSD_stage_2B_vs_3` | — | — | ok |
| `NSD_stage_2B_vs_3_PDonly` | — | — | ok |
| `NSD_stage_early_vs_late` | — | — | ok |
| `NSD_vs_HC` | — | — | ok |
| `NSD_vs_notNSD_GBA` | — | — | ok |
| `NSD_vs_notNSD_LRRK2` | — | — | ok |
| `NSD_vs_notNSD_prodromal` | — | — | ok |
| `NSD_vs_notNSD_sPD` | — | — | ok |
| `cox_cogstate_worsen` | — | — | ok |
| `cox_cogstate_worsen_x_PRS` | — | — | ok |
| `cox_moca_lt26` | — | — | ok |
| `cox_moca_lt26_x_PRS` | — | — | ok |
| `cox_nsd_2a_to_later` | — | — | ok |
| `cox_nsd_2a_to_later_x_PRS` | — | — | ok |
| `cox_nsd_2b_to_later` | — | — | ok |
| `cox_nsd_2b_to_later_x_PRS` | — | — | ok |
| `cox_pm_any` | — | — | ok |
| `cox_pm_any_x_PRS` | — | — | ok |
| `cox_pm_cog_any` | — | — | ok |
| `cox_pm_cog_any_x_PRS` | — | — | ok |
| `cox_pm_mc_any` | — | — | ok |
| `cox_pm_mc_any_x_PRS` | — | — | ok |
| `cox_stage_d` | — | — | ok |
| `cox_stage_d_x_PRS` | — | — | ok |
| `sPD_vs_GBA` | — | — | ok |
| `sPD_vs_GBA_SAAadj` | — | — | ok |
| `sPD_vs_LRRK2` | — | — | ok |
| `sPD_vs_LRRK2_SAAadj` | — | — | ok |
| `slope_lowput_ratio` | — | — | ok |
| `slope_lowput_ratio_x_PRS` | — | — | ok |
| `slope_moca` | — | — | ok |
| `slope_moca_x_PRS` | — | — | ok |
| `slope_updrs3_off` | — | — | ok |
| `slope_updrs3_off_LEDD_adj` | — | — | ok |
| `slope_updrs3_off_x_PRS` | — | — | ok |
| `slope_updrs3_off_x_PRS_LEDD_adj` | — | — | ok |
| `trajectory_CI_MOCA` | within=0% (n=2) / between=0% (n=2) | within=0% (n=2) / between=0% (n=2) | ok |
| `trajectory_CI_PI` | within=0% (n=2) / between=0% (n=2) | within=0% (n=2) / between=0% (n=2) | ok |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within=0% (n=2) / between=0% (n=2) | within=50% (n=2) / between=50% (n=2) | ⚠️ >2× nominal |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within=0% (n=2) / between=0% (n=2) | within=0% (n=2) / between=0% (n=2) | ok |
| `trajectory_HC_vs_PD` | within=0% (n=2) / between=0% (n=2) | within=50% (n=2) / between=0% (n=2) | ⚠️ >2× nominal |
| `trajectory_HC_vs_Prodromal` | within=0% (n=2) / between=0% (n=2) | within=0% (n=2) / between=0% (n=2) | ok |
| `trajectory_NSD_vs_HC` | within=0% (n=2) / between=0% (n=2) | within=50% (n=2) / between=50% (n=2) | ⚠️ >2× nominal |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within=0% (n=2) / between=0% (n=2) | within=0% (n=2) / between=50% (n=2) | ⚠️ >2× nominal |
| `trajectory_Prodromal_vs_PD` | within=0% (n=2) / between=0% (n=2) | within=0% (n=2) / between=0% (n=2) | ok |
| `trajectory_stage_2A_vs_2B` | within=0% (n=2) / between=0% (n=2) | — | ok |

**4 of 52 runs** exceed 2× the nominal rate on at least one stratum.


## PRS involvement

- **PRS-interaction runs** (name contains `_x_PRS`): **54** unique hits (reported β = `predictor:p9005_Genetic_PRS_PRS157` interaction term).
- **PRS as the looped predictor** (`p9005_Genetic_PRS_PRS157` main effect): **10** unique hits.
- **Union** (any PRS involvement): **64** unique hits.

**Interaction-run breakdown:**

| Run | # sig hits |
|---|---|
| `slope_updrs3_off_x_PRS` | 14 |
| `slope_updrs3_off_x_PRS_LEDD_adj` | 13 |
| `cox_stage_d_x_PRS` | 6 |
| `cox_pm_mc_any_x_PRS` | 4 |
| `cox_pm_cog_any_x_PRS` | 4 |
| `cox_cogstate_worsen_x_PRS` | 4 |
| `cox_moca_lt26_x_PRS` | 4 |
| `cox_pm_any_x_PRS` | 2 |
| `cox_nsd_2a_to_later_x_PRS` | 2 |
| `slope_lowput_ratio_x_PRS` | 1 |

**PRS-as-predictor breakdown:**

| Run | # sig hits |
|---|---|
| `NSD_vs_HC` | 1 |
| `NSD_stage_early_vs_late` | 1 |
| `NSD_stage_2A_vs_3` | 1 |
| `sPD_vs_GBA` | 1 |
| `sPD_vs_GBA_SAAadj` | 1 |
| `sPD_vs_LRRK2` | 1 |
| `sPD_vs_LRRK2_SAAadj` | 1 |
| `NSD_stage_2B_vs_3` | 1 |
| `NSD_vs_notNSD_prodromal` | 1 |
| `cox_stage_d` | 1 |


## Bonferroni-significant hits

Each row is one unique `(run, predictor, outcome)` tuple, deduped across EUR / AJ / META FE / META RE — the `Sources` column lists every place it was flagged significant. β/SE/P shown are from the **most-significant source** for that tuple. `α (Bonf)` is the run-wide threshold and `α (family)` the within-assay-panel one; `Lens` says which the hit clears (`both`, `run-wide`, or `family`). Sorted by P globally so the strongest signals are at the top.

For LMM runs the `Term` column says whether the hit is a **within**-subject trajectory difference or a **between**-subject cross-sectional one. `⚠` marks a numerically-degenerate fit (|β| < 1e-10), retained for transparency but never a finding.

| Run | Term | Predictor | Family | Outcome | β | SE | P | α (Bonf) | α (family) | Lens | Sources | Model |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_GPR158_NPX` | +0.182 | 0.00456 | 0.00e+00 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_KEL_NPX` | -0.239 | 0.00418 | 0.00e+00 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_MEGF9_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +4.61 | 0.211 | 1.62e-105 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `cox_pm_any` | effect | `harmonized_olink_csf_PTGDS_NPX` | `harmonized_olink_csf` | `tte_pm_any_years` | +6.65 | 0.342 | 2.34e-84 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | COX |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_NTHL1_NPX` | +0.256 | 0.0141 | 1.48e-73 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_SPAG7_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +2.76 | 0.169 | 6.81e-60 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `NSD_vs_notNSD_prodromal` | effect | `upsit` | `non_highthroughput_proteomics` | `grp_NSD_vs_notNSD_prodromal` | -1.33 | 0.0851 | 2.38e-55 | 4.32e-06 | 4.17e-03 | both | AJ, EUR, META_FE, META_RE | LOGIT |
| `NSD_vs_HC` | effect | `upsit` | `non_highthroughput_proteomics` | `grp_NSD_vs_HC` | -2.91 | 0.188 | 8.37e-54 | 4.32e-06 | 3.57e-03 | both | AJ, EUR, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_KIFC2_NPX` | -0.0972 | 0.0063 | 8.72e-54 | 4.32e-06 | 9.26e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_RBFA_NPX` | +0.159 | 0.0114 | 2.67e-44 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.18 | 0.0131 | 4.04e-43 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.179 | 0.0131 | 1.24e-42 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PGLYRP1_NPX` | -0.456 | 0.0342 | 1.44e-40 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_HC` | +2.75 | 0.207 | 2.23e-40 | 4.32e-06 | 3.91e-04 | both | AJ, EUR, META_FE, META_RE | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_MRPL24_NPX` | -0.258 | 0.0198 | 1.09e-38 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.17 | 0.0133 | 3.91e-37 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.152 | 0.012 | 1.70e-36 | 4.32e-06 | 9.26e-06 | both | AJ, EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_EPHB4_NPX` | -0.206 | 0.0163 | 2.00e-36 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_ZRANB3_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -0.989 | 0.0792 | 9.97e-36 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `NSD_vs_HC` | effect | `harmonized_olink_csf_DDC_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_HC` | +2.22 | 0.182 | 3.82e-34 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_prodromal` | +1.55 | 0.13 | 6.89e-33 | 4.32e-06 | 3.91e-04 | both | AJ, EUR, META_FE, META_RE | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PSMB11_NPX` | +0.329 | 0.0279 | 5.41e-32 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_HSPG2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +1.83 | 0.157 | 1.66e-31 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_TSSC4_NPX` | +0.0561 | 0.00507 | 2.21e-28 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_FLAD1_NPX` | +0.22 | 0.0203 | 2.08e-27 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_EIF4ENIF1_NPX` | -0.0302 | 0.0029 | 1.72e-25 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_MSR1_NPX` | -0.0684 | 0.00673 | 3.10e-24 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_TUBG2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -5.85 | 0.596 | 9.69e-23 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.142 | 0.0149 | 1.44e-21 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_C8orf33_NPX` | +0.0701 | 0.00736 | 1.80e-21 | 4.32e-06 | 9.26e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CD4_NPQ` | +0.271 | 0.0287 | 3.82e-21 | 4.32e-06 | 2.01e-04 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_sPD` | effect | `upsit` | `non_highthroughput_proteomics` | `grp_NSD_vs_notNSD_sPD` | -2.97 | 0.325 | 6.57e-20 | 4.32e-06 | 3.85e-03 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TFAP2A_NPX` | -0.0809 | 0.00891 | 1.09e-19 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.147 | 0.0159 | 1.38e-19 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_ATE1_NPX` | +0.161 | 0.0183 | 1.25e-18 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RGMB_NPX` | -0.353 | 0.034 | 1.39e-18 | 4.33e-06 | 9.25e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PAGE3_NPX` | -0.131 | 0.015 | 3.36e-18 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_LIMD1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +2.72 | 0.313 | 3.57e-18 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_KCTD19_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -0.401 | 0.0463 | 4.53e-18 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PTGDS_NPX` | +0.539 | 0.0582 | 7.82e-18 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.145 | 0.0169 | 3.46e-17 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GCFC2_NPX` | -0.138 | 0.0158 | 7.54e-17 | 4.33e-06 | 9.24e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_LRRK2` | effect | `upsit` | `non_highthroughput_proteomics` | `grp_NSD_vs_notNSD_LRRK2` | -2.03 | 0.246 | 1.99e-16 | 4.32e-06 | 3.85e-03 | both | AJ, EUR, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.147 | 0.0177 | 3.94e-16 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.173 | 0.0213 | 5.26e-16 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DBH_NPX` | +0.0977 | 0.0123 | 1.86e-15 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | within | `YEAR:grp_HC_vs_Prodromal[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ABHD14B_NPX` | -0.144 | 0.0181 | 2.14e-15 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD55_NPX` | +0.497 | 0.0591 | 2.47e-15 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.155 | 0.0192 | 2.82e-15 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_IL1R1_NPQ` | +0.138 | 0.0175 | 3.47e-15 | 4.32e-06 | 2.02e-04 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.149 | 0.0188 | 4.34e-15 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SOD1_NPX` | +0.553 | 0.0664 | 4.53e-15 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CCN3_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +10.9 | 1.41 | 1.26e-14 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NPTX2_NPQ` | +0.0434 | 0.00567 | 1.90e-14 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TNFRSF11A_NPX` | -0.335 | 0.0396 | 6.42e-14 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNASE1_NPX` | +0.54 | 0.0689 | 1.22e-13 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CST3_NPX` | +0.418 | 0.0498 | 1.38e-13 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_CCN2_NPX` | +0.126 | 0.0172 | 2.53e-13 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_HOMER2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -4.85 | 0.664 | 2.62e-13 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_SF3B4_NPX` | +0.112 | 0.0154 | 3.52e-13 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD46_NPX` | +0.478 | 0.0629 | 5.69e-13 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_NSD_vs_HC` | +0.67 | 0.0931 | 5.91e-13 | 4.32e-06 | 3.57e-03 | both | AJ, EUR, META_FE | LOGIT |
| `trajectory_CI_PI` | between | `YEAR:CI_PI[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CDV3_NPX` | +0.195 | 0.0262 | 6.10e-13 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_PLXNB2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +8.59 | 1.19 | 6.39e-13 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CD93_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +1.82 | 0.255 | 8.97e-13 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_BCAT1_NPX` | +0.107 | 0.0149 | 9.53e-13 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.158 | 0.0222 | 1.09e-12 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_FAAP20_NPX` | +0.118 | 0.0166 | 1.25e-12 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_olink_csf_DDC_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_prodromal` | +1.19 | 0.167 | 1.32e-12 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.0988 | 0.014 | 1.89e-12 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_csf_DDC_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_early_vs_late` | +0.752 | 0.107 | 1.93e-12 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LOGIT |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_TMPRSS15_NPX` | +0.0688 | 0.00978 | 2.05e-12 | 4.32e-06 | 9.26e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ANKRD52_NPX` | +0.1 | 0.0143 | 2.48e-12 | 4.32e-06 | 9.26e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.0951 | 0.0139 | 7.23e-12 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_IGFBP4_NPX` | +0.173 | 0.0201 | 7.88e-12 | 4.33e-06 | 9.26e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_TCTN3_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +4.35 | 0.639 | 9.23e-12 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `NSD_vs_notNSD_sPD` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_sPD` | +1.65 | 0.242 | 1.01e-11 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CHL1_NPX` | +0.432 | 0.0609 | 1.19e-11 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CNDP1_NPX` | +0.51 | 0.0723 | 1.65e-11 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TGOLN2_NPX` | +0.47 | 0.0669 | 1.81e-11 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.192 | 0.0287 | 2.14e-11 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_MANSC4_NPX` | +0.178 | 0.0265 | 2.21e-11 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_NSD_stage_early_vs_late` | +0.367 | 0.0549 | 2.43e-11 | 4.32e-06 | 3.57e-03 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.0884 | 0.0134 | 3.82e-11 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_3` | effect | `harmonized_olink_csf_DDC_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_2A_vs_3` | +1.09 | 0.165 | 3.83e-11 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PTGDS_NPX` | +0.45 | 0.0627 | 4.14e-11 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.109 | 0.0166 | 4.39e-11 | 4.32e-06 | 3.91e-04 | both | AJ, EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NCAM1_NPX` | +0.464 | 0.0677 | 4.83e-11 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_3` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_NSD_stage_2A_vs_3` | +0.533 | 0.0811 | 4.91e-11 | 4.32e-06 | 3.85e-03 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_RPLP1_NPX` | +0.198 | 0.0302 | 5.26e-11 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_AXL_NPX` | +0.417 | 0.0609 | 5.49e-11 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TIMP1_NPX` | +0.481 | 0.0707 | 6.68e-11 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_TSPAN15_NPX` | +0.0874 | 0.0135 | 8.37e-11 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_stage_early_vs_late` | +0.522 | 0.0806 | 9.25e-11 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SOD1_NPX` | +0.493 | 0.0703 | 9.80e-11 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.207 | 0.0321 | 1.06e-10 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_STC2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +1.66 | 0.258 | 1.06e-10 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_HSPG2_NPX` | +0.464 | 0.0693 | 1.32e-10 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SPARCL1_NPX` | +0.467 | 0.0703 | 1.82e-10 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_C1QTNF1_NPX` | +0.456 | 0.069 | 2.15e-10 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_MAEA_NPX` | +0.281 | 0.0347 | 2.46e-10 | 4.33e-06 | 9.27e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ERG_NPX` | -0.0767 | 0.0122 | 2.91e-10 | 4.32e-06 | 9.26e-06 | both | AJ, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DNAAF2_NPX` | +0.186 | 0.0298 | 4.01e-10 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNASE1_NPX` | +0.483 | 0.0717 | 4.09e-10 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GSR_NPX` | +0.404 | 0.0624 | 4.74e-10 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | within | `YEAR:grp_HC_vs_Prodromal[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_FIG4_NPX` | -0.155 | 0.024 | 5.39e-10 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `sPD_vs_GBA` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_sPD_vs_GBA` | +1.24 | 0.2 | 5.43e-10 | 4.32e-06 | 3.85e-03 | both | AJ, EUR, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.197 | 0.0318 | 6.06e-10 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_MRC1_NPX` | +0.48 | 0.0748 | 6.54e-10 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PTGDS_NPX` | +0.397 | 0.0594 | 8.19e-10 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_GABRP_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +4.09 | 0.667 | 8.57e-10 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_UQCC2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -0.339 | 0.0554 | 9.34e-10 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_TSPAN8_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -3.61 | 0.592 | 1.10e-09 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNASE6_NPX` | +0.43 | 0.068 | 1.16e-09 | 4.32e-06 | 9.24e-06 | both | AJ, EUR | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CST3_NPX` | +0.374 | 0.0547 | 1.20e-09 | 4.33e-06 | 9.27e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_PAIP1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +1.97 | 0.323 | 1.20e-09 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_IL6R_NPX` | +0.452 | 0.0727 | 2.11e-09 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD163_NPX` | +0.461 | 0.0747 | 2.53e-09 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_3` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_stage_2A_vs_3` | +0.688 | 0.116 | 2.62e-09 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DEFB4A_DEFB4B_NPX` | -0.208 | 0.0315 | 2.75e-09 | 4.33e-06 | 9.27e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_olink_plasma_ITGAV_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_HC` | -0.641 | 0.108 | 3.19e-09 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_plasma_ITGAV_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_early_vs_late` | -0.462 | 0.0782 | 3.56e-09 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ATXN7L2_NPX` | +0.141 | 0.0241 | 4.86e-09 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_LRRK2` | +1.7 | 0.291 | 4.96e-09 | 4.32e-06 | 3.91e-04 | both | AJ, META_FE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_AGRN_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +4.58 | 0.785 | 5.48e-09 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD46_NPX` | +0.481 | 0.0774 | 5.64e-09 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_ATXN2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -3.22 | 0.553 | 5.65e-09 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `upsit` | `non_highthroughput_proteomics` | `grp_NSD_vs_notNSD_GBA` | -3.51 | 0.604 | 6.40e-09 | 4.35e-06 | 5.00e-03 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CST3_NPX` | +0.413 | 0.0687 | 6.68e-09 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CTBS_NPX` | +0.449 | 0.0752 | 7.43e-09 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_MZT1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +2.33 | 0.405 | 8.67e-09 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NEFL_NPQ` | +0.0572 | 0.00995 | 8.75e-09 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SSC5D_NPX` | +0.417 | 0.0702 | 9.00e-09 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_C1S_NPX` | -0.297 | 0.0482 | 1.01e-08 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_GHRL_NPX` | -0.187 | 0.0308 | 1.42e-08 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DDC_NPX` | +0.141 | 0.0247 | 1.51e-08 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_UBXN2B_NPX` | -0.205 | 0.033 | 1.57e-08 | 4.33e-06 | 9.27e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_KIFC2_NPX` | +0.0955 | 0.0158 | 1.65e-08 | 4.33e-06 | 9.26e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_BRDT_NPX` | +0.215 | 0.0382 | 1.76e-08 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CYP4B1_NPX` | +0.219 | 0.039 | 1.85e-08 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PLB1_NPX` | +0.226 | 0.0379 | 1.87e-08 | 4.33e-06 | 9.24e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_KIAA1328_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.52 | 0.448 | 1.92e-08 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_ADA2_NPX` | +0.468 | 0.0786 | 1.98e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_ODAD4_NPX` | -0.00974 | 0.00174 | 2.02e-08 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `cox_pm_cog_any` | effect | `upsit` | `non_highthroughput_proteomics` | `tte_pm_cog_any_years` | -0.368 | 0.0657 | 2.10e-08 | 4.32e-06 | 3.57e-03 | both | AJ, EUR, META_FE, META_RE | COX |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_HAVCR2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -5.51 | 0.984 | 2.21e-08 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PODXL_NPX` | -0.463 | 0.0775 | 2.30e-08 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PPP1R35_NPX` | -0.0404 | 0.00723 | 2.38e-08 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NCAM1_NPX` | +0.436 | 0.0739 | 2.71e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DAP_NPX` | -0.056 | 0.0101 | 2.74e-08 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_PI` | within | `YEAR:CI_PI[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_VIPR1_NPX` | -0.0982 | 0.0173 | 2.75e-08 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_olink_plasma_DDC_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_notNSD_LRRK2` | +1.08 | 0.194 | 2.77e-08 | 4.32e-06 | 9.26e-06 | both | META_FE | LOGIT |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NEFL_NPQ` | +0.055 | 0.0099 | 2.78e-08 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_NPHS1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.15 | 0.387 | 2.83e-08 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNASE6_NPX` | +0.42 | 0.0716 | 3.13e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_PARK7_NPQ` | +0.119 | 0.0214 | 3.14e-08 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_olink_plasma_PEPD_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_HC` | -0.518 | 0.0936 | 3.18e-08 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_IL17A_NPX` | +0.322 | 0.0548 | 3.50e-08 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TGOLN2_NPX` | +0.419 | 0.0717 | 3.55e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_2B_PDonly` | effect | `harmonized_nulisa_inf_csf_IL17C_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_stage_2A_vs_2B` | +1.16 | 0.21 | 3.71e-08 | 4.33e-06 | 2.20e-04 | both | EUR, META_FE, META_RE | LOGIT |
| `NSD_vs_HC` | effect | `harmonized_nulisa_inf_plasma_TNFSF4_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_vs_HC` | -0.706 | 0.128 | 3.78e-08 | 4.32e-06 | 2.01e-04 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DDC_NPX` | +0.105 | 0.019 | 4.01e-08 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SMOC1_NPX` | -0.595 | 0.102 | 4.21e-08 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_FGD2_NPX` | +0.295 | 0.0505 | 4.37e-08 | 4.33e-06 | 9.25e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ABHD14B_NPX` | +0.0864 | 0.0158 | 4.66e-08 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CPM_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +4.51 | 0.826 | 4.71e-08 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNF4_NPX` | -0.354 | 0.0638 | 4.85e-08 | 4.32e-06 | 9.24e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `sPD_vs_GBA_SAAadj` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_sPD_vs_GBA` | +1.14 | 0.209 | 5.17e-08 | 4.32e-06 | 3.85e-03 | both | AJ, EUR, META_FE | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_MMP9_NPQ` | -0.19 | 0.0347 | 5.23e-08 | 4.32e-06 | 2.02e-04 | both | EUR | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CEP170B_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.94 | 0.54 | 5.31e-08 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_AXL_NPX` | +0.402 | 0.0699 | 5.39e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD55_NPX` | +0.434 | 0.0756 | 5.63e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_olink_csf_DDC_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_LRRK2` | +1.77 | 0.326 | 5.83e-08 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PCOLCE_NPX` | +0.376 | 0.0676 | 6.42e-08 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_3` | effect | `harmonized_olink_plasma_CA3_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_2A_vs_3` | +0.668 | 0.125 | 8.74e-08 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SOD1_NPX` | +0.413 | 0.0725 | 8.92e-08 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.169 | 0.0314 | 9.47e-08 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_ART3_NPX` | +0.369 | 0.0672 | 9.66e-08 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_plasma_CD276_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_early_vs_late` | +0.429 | 0.0805 | 9.78e-08 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SPARCL1_NPX` | +0.416 | 0.0739 | 9.84e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_C1QTNF1_NPX` | +0.423 | 0.0752 | 9.97e-08 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_nulisa_inf_plasma_CD276_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_stage_early_vs_late` | +0.417 | 0.0784 | 1.01e-07 | 4.32e-06 | 2.01e-04 | both | EUR, META_FE | LOGIT |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_olink_csf_ATOX1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_prodromal` | +1.52 | 0.286 | 1.05e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_KIT_NPX` | +0.371 | 0.0679 | 1.09e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_SOD1_NPQ` | +0.0872 | 0.0164 | 1.20e-07 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CHI3L1_NPX` | +0.42 | 0.0772 | 1.23e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_plasma_DDC_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_notNSD_GBA` | +1.27 | 0.241 | 1.24e-07 | 4.32e-06 | 9.26e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_MYSM1_NPX` | +0.146 | 0.026 | 1.24e-07 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CTBS_NPX` | +0.427 | 0.0766 | 1.25e-07 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CUX2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -0.146 | 0.0276 | 1.25e-07 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_S100A12_NPX` | +0.142 | 0.0268 | 1.28e-07 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_nulisa_cns_plasma_CHIT1_NPQ` | `harmonized_nulisa_cns_plasma` | `slope_updrs3_off` | -0.858 | 0.161 | 1.34e-07 | 4.32e-06 | 3.94e-04 | both | EUR | OLS |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_KIF2B_NPX` | -0.245 | 0.0443 | 1.50e-07 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNASE1_NPX` | +0.412 | 0.0737 | 1.52e-07 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DDC_NPX` | +0.101 | 0.0191 | 1.56e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CNDP1_NPX` | +0.422 | 0.0764 | 1.56e-07 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SNX33_NPX` | -0.194 | 0.0366 | 1.64e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CHL1_NPX` | +0.408 | 0.0741 | 1.75e-07 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_F11_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +1.09 | 0.21 | 1.78e-07 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_TCHHL1_NPX` | -0.388 | 0.0691 | 1.81e-07 | 4.33e-06 | 9.26e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CTNNA1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +2.55 | 0.489 | 1.85e-07 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_nulisa_cns_plasma_CHIT1_NPQ` | `harmonized_nulisa_cns_plasma` | `slope_updrs3_off` | -0.847 | 0.161 | 1.86e-07 | 4.32e-06 | 3.94e-04 | both | EUR | OLS |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_DDC_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +4.25 | 0.819 | 2.11e-07 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TIMP1_NPX` | +0.403 | 0.0739 | 2.20e-07 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_SLMAP_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.84 | 0.549 | 2.27e-07 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_VIM_NPX` | -0.251 | 0.0486 | 2.32e-07 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRL_NPX` | -0.188 | 0.0365 | 2.49e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_OPHN1_NPX` | +0.172 | 0.0334 | 2.52e-07 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRSS2_NPX` | +0.0428 | 0.00827 | 2.59e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_VASN_NPX` | +0.337 | 0.0637 | 2.59e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DDC_NPX` | +0.0999 | 0.0193 | 2.85e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PCDH9_NPX` | +0.233 | 0.04 | 2.90e-07 | 4.33e-06 | 9.26e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_HSPG2_NPX` | +0.394 | 0.0732 | 3.05e-07 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `cox_moca_lt26` | effect | `harmonized_olink_plasma_C2_NPX` | `harmonized_olink_plasma` | `tte_moca_lt26_years` | +3.33 | 0.652 | 3.11e-07 | 4.34e-06 | 9.26e-06 | both | AJ | COX |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.188 | 0.0366 | 3.50e-07 | 4.32e-06 | 9.26e-06 | both | EUR | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_nulisa_inf_csf_IL1B_NPQ` | `harmonized_nulisa_inf_csf` | `slope_updrs3_off` | +0.932 | 0.181 | 3.67e-07 | 4.32e-06 | 2.02e-04 | both | EUR, META_FE | OLS |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NEO1_NPX` | -0.495 | 0.0922 | 3.70e-07 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_USP28_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.23 | 0.438 | 3.71e-07 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_nulisa_inf_csf_IL1B_NPQ` | `harmonized_nulisa_inf_csf` | `slope_updrs3_off` | +0.93 | 0.181 | 3.84e-07 | 4.32e-06 | 2.02e-04 | both | EUR, META_FE | OLS |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_SPTBN5_NPX` | -0.188 | 0.0353 | 3.91e-07 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.184 | 0.0361 | 3.94e-07 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_VWC2L_NPX` | -0.366 | 0.0685 | 4.17e-07 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_MYO5B_NPX` | +0.0796 | 0.0158 | 4.31e-07 | 4.32e-06 | 9.23e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.109 | 0.0216 | 4.36e-07 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_FCRL3_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.68 | 0.532 | 4.61e-07 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_FOLR3_NPX` | -0.169 | 0.0334 | 4.90e-07 | 4.32e-06 | 9.26e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.0669 | 0.0133 | 4.92e-07 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_EFEMP1_NPX` | +0.324 | 0.0629 | 4.95e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRL_NPX` | -0.19 | 0.0378 | 5.01e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_plasma_CA3_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_early_vs_late` | +0.4 | 0.0797 | 5.04e-07 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_CEP164_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.768 | 0.153 | 5.17e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | OLS |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CFD_NPX` | -0.489 | 0.0824 | 5.32e-07 | 4.34e-06 | 9.28e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_ADA2_NPX` | +0.425 | 0.0829 | 5.62e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.151 | 0.03 | 5.63e-07 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_NYAP1_NPX` | -0.138 | 0.0272 | 5.74e-07 | 4.32e-06 | 9.24e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_PRDX6_NPQ` | +0.0787 | 0.0158 | 6.05e-07 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `cox_pm_mc_any_x_PRS` | effect | `harmonized_olink_plasma_PNMA2_NPX` | `harmonized_olink_plasma` | `tte_pm_mc_any_years` | -0.441 | 0.0884 | 6.24e-07 | 4.32e-06 | 9.26e-06 | both | AJ | COX |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRL_NPX` | -0.175 | 0.035 | 6.33e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `cox_stage_d` | effect | `upsit` | `non_highthroughput_proteomics` | `tte_stage_d_years` | -0.342 | 0.0686 | 6.34e-07 | 4.32e-06 | 4.17e-03 | both | EUR, META_FE, META_RE | COX |
| `NSD_stage_2A_vs_2B_PDonly` | effect | `harmonized_nulisa_inf_csf_IL13RA2_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_stage_2A_vs_2B` | +2.38 | 0.479 | 6.46e-07 | 4.33e-06 | 2.20e-04 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_EPHB6_NPX` | -0.467 | 0.089 | 6.47e-07 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_csf_TNFRSF11B_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_early_vs_late` | +0.742 | 0.149 | 6.50e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD93_NPX` | +0.349 | 0.0684 | 6.50e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_IL18BP_NPX` | +0.334 | 0.0655 | 6.67e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_SYNGAP1_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -1.08 | 0.214 | 6.88e-07 | 4.32e-06 | 9.23e-06 | both | EUR | OLS |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_csf_CEP164_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.762 | 0.154 | 7.06e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | OLS |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_olink_csf_STMN2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_prodromal` | -0.568 | 0.115 | 7.50e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PCDH7_NPX` | -0.469 | 0.09 | 7.51e-07 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_csf_SYNGAP1_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -1.07 | 0.214 | 7.57e-07 | 4.32e-06 | 9.23e-06 | both | EUR | OLS |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_LATS1_NPX` | +0.23 | 0.0442 | 7.74e-07 | 4.33e-06 | 9.25e-06 | both | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_STARD13_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -0.904 | 0.183 | 7.77e-07 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRL_NPX` | -0.155 | 0.0313 | 7.83e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ABHD14B_NPX` | -0.123 | 0.025 | 7.91e-07 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CCN1_NPX` | -0.152 | 0.0305 | 7.96e-07 | 4.32e-06 | 9.23e-06 | both | EUR | LMM_RS_lbfgs |
| `cox_stage_d` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `tte_stage_d_years` | +0.472 | 0.0957 | 8.32e-07 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE, META_RE | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_CCL8_NPX` | -0.607 | 0.117 | 8.43e-07 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_SCRG1_NPX` | +0.161 | 0.0325 | 8.51e-07 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NPTX2_NPQ` | +0.0423 | 0.0086 | 8.53e-07 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_csf_PRELP_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_early_vs_late` | +0.744 | 0.151 | 8.83e-07 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_GBA` | +4.2 | 0.856 | 9.14e-07 | 4.35e-06 | 3.91e-04 | both | AJ, META_FE, META_RE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_SOX14_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +0.472 | 0.0961 | 9.18e-07 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PROC_NPX` | +0.36 | 0.0717 | 9.37e-07 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.165 | 0.0335 | 9.44e-07 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.177 | 0.0359 | 9.48e-07 | 4.32e-06 | 9.26e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_NME3_NPX` | +0.167 | 0.0341 | 9.82e-07 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_TRIAP1_NPX` | -0.526 | 0.102 | 1.00e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ASS1_NPX` | -0.0805 | 0.0165 | 1.00e-06 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ASS1_NPX` | -0.0811 | 0.0166 | 1.00e-06 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_TSGA10_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.773 | 0.158 | 1.01e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | OLS |
| `NSD_stage_2A_vs_3` | effect | `harmonized_olink_csf_TMED1_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_2A_vs_3` | -0.52 | 0.106 | 1.04e-06 | 4.32e-06 | 9.23e-06 | both | META_FE | LOGIT |
| `sPD_vs_LRRK2` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_sPD_vs_LRRK2` | +2.29 | 0.47 | 1.06e-06 | 4.32e-06 | 3.85e-03 | both | EUR, META_FE, META_RE | LOGIT |
| `CI_PI_baseline` | effect | `harmonized_olink_plasma_VGF_NPX` | `harmonized_olink_plasma` | `CI_PI` | -0.608 | 0.125 | 1.07e-06 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_S100A12_NPQ` | +0.101 | 0.0208 | 1.12e-06 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `ptau` | `non_highthroughput_proteomics` | `grp_NSD_vs_HC` | -0.587 | 0.121 | 1.14e-06 | 4.32e-06 | 3.57e-03 | both | EUR, META_FE, META_RE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SPINK2_NPX` | -0.373 | 0.0731 | 1.19e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_TNFSF14_NPQ` | -0.0999 | 0.0206 | 1.20e-06 | 4.32e-06 | 2.01e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_HECW2_NPX` | -0.252 | 0.052 | 1.22e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRR20B_NPX` | +0.185 | 0.0382 | 1.23e-06 | 4.32e-06 | 9.23e-06 | both | META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_BAG3_NPX` | +0.0771 | 0.0158 | 1.26e-06 | 4.32e-06 | 9.26e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_REG3A_NPX` | +0.0636 | 0.0131 | 1.26e-06 | 4.32e-06 | 9.26e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.175 | 0.036 | 1.27e-06 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_THPO_NPX` | -0.0965 | 0.0199 | 1.30e-06 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_MCAM_NPX` | +0.332 | 0.0669 | 1.30e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_DDC_NPX` | +0.157 | 0.0323 | 1.30e-06 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SV2A_NPX` | +0.472 | 0.087 | 1.31e-06 | 4.33e-06 | 9.26e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_DDC_NPX` | +0.157 | 0.0323 | 1.32e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `slope_moca` | effect | `abeta` | `non_highthroughput_proteomics` | `slope_moca` | +0.204 | 0.0422 | 1.38e-06 | 4.32e-06 | 3.57e-03 | both | AJ, EUR, META_FE, META_RE | OLS |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_IL1RL1_NPQ` | -0.19 | 0.0345 | 1.42e-06 | 4.33e-06 | 2.03e-04 | both | AJ, META_FE | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_nulisa_cns_csf_CHIT1_NPQ` | `harmonized_nulisa_cns_csf` | `slope_updrs3_off` | -0.832 | 0.171 | 1.42e-06 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE | OLS |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_DDC_NPX` | +0.202 | 0.0418 | 1.42e-06 | 4.32e-06 | 9.23e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_ST6GAL2_NPX` | +0.136 | 0.0283 | 1.45e-06 | 4.32e-06 | 9.23e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TGFA_NPX` | -0.141 | 0.0291 | 1.45e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GFRA2_NPX` | -0.49 | 0.0969 | 1.46e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD46_NPX` | +0.389 | 0.0767 | 1.46e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_LMNB1_NPX` | +0.175 | 0.0361 | 1.48e-06 | 4.32e-06 | 9.23e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CA14_NPX` | +0.568 | 0.105 | 1.50e-06 | 4.33e-06 | 9.26e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_VNN2_NPX` | -0.304 | 0.0603 | 1.52e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_nulisa_inf_csf_IL11_NPQ` | `harmonized_nulisa_inf_csf` | `slope_updrs3_off` | +1.03 | 0.211 | 1.53e-06 | 4.32e-06 | 2.02e-04 | both | EUR | OLS |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_NME3_NPX` | +0.106 | 0.0217 | 1.55e-06 | 4.32e-06 | 9.26e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_3` | effect | `harmonized_nulisa_cns_plasma_pTau-181_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_stage_2A_vs_3` | +0.879 | 0.183 | 1.57e-06 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NCAM1_NPX` | +0.38 | 0.0752 | 1.60e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.0581 | 0.0121 | 1.66e-06 | 4.32e-06 | 9.26e-06 | both | META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNPC3_NPX` | -0.282 | 0.056 | 1.67e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_CXCL3_NPX` | -0.373 | 0.0741 | 1.68e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_plasma_GOLGA8B_NPX` | `harmonized_olink_plasma` | `slope_updrs3_off` | -0.875 | 0.181 | 1.71e-06 | 4.32e-06 | 9.26e-06 | both | EUR | OLS |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CEP170_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.66 | 0.556 | 1.72e-06 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_csf_TSGA10_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.768 | 0.161 | 1.74e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | OLS |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_TSLP_NPQ` | +0.181 | 0.0375 | 1.81e-06 | 4.33e-06 | 2.01e-04 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_BEND2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +0.972 | 0.204 | 1.86e-06 | 4.35e-06 | 9.37e-06 | both | AJ, META_FE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SOD3_NPX` | +0.423 | 0.0869 | 1.90e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_nulisa_cns_plasma_pTau-231_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_stage_early_vs_late` | +0.631 | 0.132 | 1.91e-06 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LOGIT |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_TFF3_NPX` | +0.0643 | 0.0135 | 1.91e-06 | 4.32e-06 | 9.26e-06 | both | META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRL_NPX` | -0.159 | 0.0334 | 1.93e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `sPD_vs_LRRK2` | effect | `harmonized_nulisa_inf_csf_TLR3_NPQ` | `harmonized_nulisa_inf_csf` | `grp_sPD_vs_LRRK2` | +1.7 | 0.358 | 1.95e-06 | 4.32e-06 | 2.02e-04 | both | EUR, META_FE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_ADD1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -3.21 | 0.676 | 2.01e-06 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TGOLN2_NPX` | +0.369 | 0.0738 | 2.01e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ADAMTS13_NPX` | -0.0651 | 0.0137 | 2.02e-06 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NEFL_NPQ` | +0.0461 | 0.00972 | 2.06e-06 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNASE4_NPX` | +0.322 | 0.0663 | 2.10e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD163_NPX` | +0.394 | 0.0795 | 2.11e-06 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_nulisa_inf_csf_IL11_NPQ` | `harmonized_nulisa_inf_csf` | `slope_updrs3_off` | +1.01 | 0.211 | 2.13e-06 | 4.32e-06 | 2.02e-04 | both | EUR | OLS |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_plasma_GOLGA8B_NPX` | `harmonized_olink_plasma` | `slope_updrs3_off` | -0.865 | 0.181 | 2.18e-06 | 4.32e-06 | 9.26e-06 | both | EUR | OLS |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TGFBI_NPX` | +0.309 | 0.0638 | 2.18e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CD164_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +6.05 | 1.28 | 2.18e-06 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_nulisa_cns_csf_CHIT1_NPQ` | `harmonized_nulisa_cns_csf` | `slope_updrs3_off` | -0.814 | 0.17 | 2.21e-06 | 4.32e-06 | 3.91e-04 | both | EUR, META_FE | OLS |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CTSD_NPX` | +0.387 | 0.0799 | 2.24e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_NCK2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -1.71 | 0.361 | 2.29e-06 | 4.35e-06 | 9.37e-06 | both | AJ | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CCL15_NPX` | -0.119 | 0.0251 | 2.40e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_ALPG_NPX` | -0.0832 | 0.0176 | 2.43e-06 | 4.32e-06 | 9.23e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_PI` | between | `YEAR:CI_PI[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_MRPL58_NPX` | -0.196 | 0.0414 | 2.44e-06 | 4.32e-06 | 9.23e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PROC_NPX` | -0.407 | 0.0822 | 2.46e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CENPJ_NPX` | +0.0673 | 0.0142 | 2.48e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_GDI1_NPQ` | +0.105 | 0.0222 | 2.48e-06 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NEFL_NPQ` | +0.0468 | 0.00995 | 2.48e-06 | 4.32e-06 | 3.94e-04 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_DNAAF8_NPX` | +0.254 | 0.0537 | 2.56e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_MKNK1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -2.03 | 0.431 | 2.59e-06 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SMOC1_NPX` | +0.53 | 0.101 | 2.63e-06 | 4.33e-06 | 9.26e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SOD1_NPX` | +0.358 | 0.0712 | 2.67e-06 | 4.33e-06 | 9.27e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_CELA2A_NPX` | -0.0521 | 0.0111 | 2.70e-06 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `sPD_vs_LRRK2_SAAadj` | effect | `harmonized_nulisa_inf_csf_TLR3_NPQ` | `harmonized_nulisa_inf_csf` | `grp_sPD_vs_LRRK2` | +1.88 | 0.4 | 2.73e-06 | 4.32e-06 | 2.02e-04 | both | EUR, META_FE | LOGIT |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_PAMR1_NPX` | -0.0815 | 0.0174 | 2.75e-06 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_NME3_NPX` | +0.109 | 0.0233 | 2.81e-06 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_MFAP3_NPX` | +0.253 | 0.054 | 2.94e-06 | 4.32e-06 | 9.23e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ASS1_NPX` | -0.076 | 0.0163 | 2.94e-06 | 4.32e-06 | 9.26e-06 | both | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_DNAJC6_NPX` | -0.514 | 0.105 | 3.04e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.155 | 0.033 | 3.16e-06 | 4.32e-06 | 9.26e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_CTRB1_NPX` | -0.0465 | 0.00997 | 3.17e-06 | 4.32e-06 | 9.26e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_MDM1_NPX` | +0.227 | 0.0487 | 3.21e-06 | 4.32e-06 | 9.23e-06 | both | META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GFRA2_NPX` | +0.449 | 0.0868 | 3.26e-06 | 4.33e-06 | 9.26e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_EPHA10_NPX` | -0.419 | 0.086 | 3.33e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_HLA-DRA_NPQ` | +0.0785 | 0.0169 | 3.34e-06 | 4.32e-06 | 2.02e-04 | both | AJ, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_CXCL8_NPX` | -0.542 | 0.111 | 3.37e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GPI_NPX` | +0.392 | 0.0827 | 3.44e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_PLEKHG3_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.8 | 0.172 | 3.50e-06 | 4.32e-06 | 9.23e-06 | both | META_FE, META_RE | OLS |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_MMP3_NPX` | -0.142 | 0.0307 | 3.53e-06 | 4.32e-06 | 9.23e-06 | both | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_THPO_NPQ` | +0.119 | 0.0256 | 3.57e-06 | 4.32e-06 | 2.01e-04 | both | EUR, META_FE | LMM_RS_lbfgs |
| `slope_moca` | effect | `harmonized_olink_csf_PTPRN2_NPX` | `harmonized_olink_csf` | `slope_moca` | +0.623 | 0.135 | 3.59e-06 | 4.32e-06 | 9.23e-06 | both | META_FE | OLS |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_IL6R_NPX` | +0.384 | 0.0796 | 3.64e-06 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `cox_pm_any` | effect | `harmonized_olink_csf_JAM2_NPX` | `harmonized_olink_csf` | `tte_pm_any_years` | -1.71 | 0.369 | 3.64e-06 | 4.32e-06 | 9.24e-06 | both | AJ | COX |
| `trajectory_HC_vs_Prodromal` | within | `YEAR:grp_HC_vs_Prodromal[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_PTGDS_NPX` | -0.29 | 0.0613 | 3.68e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CSF1R_NPX` | +0.355 | 0.0751 | 3.76e-06 | 4.32e-06 | 9.24e-06 | both | AJ, EUR | LMM_RS_lbfgs |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_olink_csf_TDGF1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_prodromal` | -0.769 | 0.166 | 3.83e-06 | 4.32e-06 | 9.23e-06 | both | EUR | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_CRB2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +1.41 | 0.306 | 3.83e-06 | 4.32e-06 | 9.24e-06 | both | AJ, META_FE, META_RE | LOGIT |
| `cox_pm_cog_any_x_PRS` | effect | `harmonized_olink_csf_GCNT4_NPX` | `harmonized_olink_csf` | `tte_pm_cog_any_years` | +0.571 | 0.124 | 3.83e-06 | 4.33e-06 | 9.26e-06 | both | AJ | COX |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CNTN4_NPX` | +0.0867 | 0.0187 | 3.85e-06 | 4.32e-06 | 9.23e-06 | both | EUR | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_IGFBP6_NPX` | +0.358 | 0.0745 | 3.91e-06 | 4.33e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD55_NPX` | +0.381 | 0.0786 | 4.00e-06 | 4.33e-06 | 9.25e-06 | both | AJ | LMM_RS_lbfgs |
| `cox_stage_d_x_PRS` | effect | `harmonized_olink_plasma_HNRNPUL1_NPX` | `harmonized_olink_plasma` | `tte_stage_d_years` | +0.569 | 0.123 | 4.01e-06 | 4.32e-06 | 9.26e-06 | both | EUR | COX |
| `NSD_vs_HC` | effect | `harmonized_olink_csf_WIF1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_HC` | -0.905 | 0.196 | 4.04e-06 | 4.32e-06 | 9.23e-06 | both | META_FE, META_RE | LOGIT |
| `NSD_stage_2A_vs_3` | effect | `harmonized_olink_csf_PRELP_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_2A_vs_3` | +1.18 | 0.255 | 4.05e-06 | 4.32e-06 | 9.23e-06 | both | EUR | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TGFA_NPX` | +0.111 | 0.024 | 4.06e-06 | 4.32e-06 | 9.23e-06 | both | EUR | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_3` | effect | `harmonized_olink_csf_MEF2C_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_2A_vs_3` | +0.498 | 0.108 | 4.10e-06 | 4.32e-06 | 9.23e-06 | both | EUR | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CTSZ_NPX` | +0.365 | 0.0777 | 4.14e-06 | 4.32e-06 | 9.24e-06 | both | AJ | LMM_RS_lbfgs |
| `cox_cogstate_worsen_x_PRS` | effect | `harmonized_olink_csf_CYP26C1_NPX` | `harmonized_olink_csf` | `tte_cogstate_worsen_years` | -0.346 | 0.0751 | 4.14e-06 | 4.32e-06 | 9.24e-06 | both | META_FE, META_RE | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_POMC_NPX` | -0.473 | 0.0983 | 4.35e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PRL_NPX` | -0.144 | 0.0313 | 4.38e-06 | 4.32e-06 | 9.23e-06 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `cox_stage_d_x_PRS` | effect | `harmonized_olink_plasma_GLYR1_NPX` | `harmonized_olink_plasma` | `tte_stage_d_years` | +0.523 | 0.114 | 4.48e-06 | 4.32e-06 | 9.26e-06 | family | EUR, META_FE | COX |
| `cox_pm_any` | effect | `harmonized_olink_plasma_SAMD15_NPX` | `harmonized_olink_plasma` | `tte_pm_any_years` | +0.686 | 0.15 | 4.55e-06 | 4.32e-06 | 9.26e-06 | family | AJ | COX |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_PSMB11_NPX` | +0.0615 | 0.0134 | 4.56e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_THPO_NPQ` | +0.152 | 0.0332 | 4.59e-06 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_GEMIN2_NPX` | +0.0786 | 0.0171 | 4.59e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_FLT1_NPQ` | +0.0946 | 0.0205 | 4.61e-06 | 4.32e-06 | 2.03e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_KNG1_NPQ` | +0.148 | 0.0322 | 4.64e-06 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_csf_PLEKHG3_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.8 | 0.175 | 4.69e-06 | 4.32e-06 | 9.23e-06 | family | META_FE, META_RE | OLS |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_DNAJC6_NPX` | +0.538 | 0.106 | 4.72e-06 | 4.33e-06 | 9.26e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_olink_plasma_GOLGA4_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_HC` | +0.628 | 0.137 | 4.72e-06 | 4.32e-06 | 9.26e-06 | family | META_FE | LOGIT |
| `cox_pm_cog_any` | effect | `harmonized_olink_csf_DYNC1H1_NPX` | `harmonized_olink_csf` | `tte_pm_cog_any_years` | -0.836 | 0.183 | 4.74e-06 | 4.32e-06 | 9.24e-06 | family | AJ | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_BAG3_NPX` | +0.109 | 0.0238 | 4.86e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CD248_NPX` | -0.112 | 0.0243 | 4.88e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NID1_NPX` | +0.323 | 0.0691 | 4.90e-06 | 4.32e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `cox_moca_lt26_x_PRS` | effect | `harmonized_olink_csf_CEP350_NPX` | `harmonized_olink_csf` | `tte_moca_lt26_years` | +0.338 | 0.0739 | 4.94e-06 | 4.32e-06 | 9.23e-06 | family | META_FE, META_RE | COX |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_LTBP3_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | -4.51 | 0.987 | 4.99e-06 | 4.32e-06 | 9.24e-06 | family | META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CTSZ_NPX` | +0.414 | 0.0872 | 5.00e-06 | 4.33e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_plasma_ANKRD2_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_early_vs_late` | +0.331 | 0.0726 | 5.05e-06 | 4.32e-06 | 9.26e-06 | family | META_FE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SUSD4_NPX` | -0.0928 | 0.0203 | 5.06e-06 | 4.32e-06 | 9.23e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_HMOX2_NPX` | +0.109 | 0.0239 | 5.21e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_GIGYF2_NPX` | +0.0871 | 0.019 | 5.21e-06 | 4.32e-06 | 9.26e-06 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_ISLR2_NPX` | -0.435 | 0.0913 | 5.25e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_CPVL_NPX` | +0.0568 | 0.0125 | 5.37e-06 | 4.32e-06 | 9.26e-06 | family | META_FE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `tau` | `non_highthroughput_proteomics` | `grp_NSD_vs_HC` | -0.5 | 0.11 | 5.41e-06 | 4.32e-06 | 3.57e-03 | family | EUR, META_FE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_LRTM2_NPX` | -0.402 | 0.0846 | 5.44e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `cox_moca_lt26` | effect | `harmonized_nulisa_cns_csf_PDGFRB_NPQ` | `harmonized_nulisa_cns_csf` | `tte_moca_lt26_years` | -0.937 | 0.206 | 5.50e-06 | 4.34e-06 | 3.91e-04 | family | AJ | COX |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_IL18BP_NPX` | +0.338 | 0.0714 | 5.53e-06 | 4.33e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NCS1_NPX` | -0.401 | 0.0844 | 5.56e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_nulisa_cns_plasma_pTau-181_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_stage_early_vs_late` | +0.525 | 0.116 | 5.56e-06 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE | LOGIT |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_CBX4_NPX` | +0.106 | 0.0224 | 5.58e-06 | 4.33e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PALM_NPX` | -0.405 | 0.0853 | 5.62e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_olink_plasma_DDC_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_HC` | +0.42 | 0.0925 | 5.65e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_COL9A1_NPX` | +0.0545 | 0.012 | 5.71e-06 | 4.32e-06 | 9.26e-06 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ACTA2_NPX` | +0.0696 | 0.0153 | 5.71e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_HDGFL2_NPX` | -0.235 | 0.0495 | 5.72e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_ODAD4_NPX` | +0.0741 | 0.0163 | 5.82e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_TSHZ2_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.739 | 0.163 | 5.83e-06 | 4.32e-06 | 9.23e-06 | family | META_FE | OLS |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SCARB2_NPX` | -0.447 | 0.0945 | 5.88e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GPX1_NPX` | +0.122 | 0.0271 | 6.11e-06 | 4.32e-06 | 9.23e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_LMOD1_NPX` | +0.055 | 0.0122 | 6.14e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TFF3_NPX` | +0.333 | 0.0723 | 6.22e-06 | 4.32e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_olink_csf_GAS2_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_notNSD_GBA` | +2.01 | 0.446 | 6.25e-06 | 4.35e-06 | 9.37e-06 | family | AJ | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_ICAM1_NPX` | +0.335 | 0.0727 | 6.25e-06 | 4.32e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GCC1_NPX` | +0.22 | 0.0484 | 6.34e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_KDR_NPX` | -0.398 | 0.0844 | 6.38e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ZER1_NPX` | +0.164 | 0.0361 | 6.41e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_UBXN2B_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.816 | 0.179 | 6.44e-06 | 4.32e-06 | 9.23e-06 | family | EUR | OLS |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_VEGFA_NPX` | -0.104 | 0.0232 | 6.51e-06 | 4.32e-06 | 9.23e-06 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CHL1_NPX` | +0.359 | 0.0761 | 6.61e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ADAMTS13_NPX` | -0.0684 | 0.0152 | 6.62e-06 | 4.32e-06 | 9.26e-06 | family | META_FE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_olink_csf_SDC1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_HC` | -0.845 | 0.188 | 6.65e-06 | 4.32e-06 | 9.23e-06 | family | META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CHI3L1_NPX` | +0.364 | 0.0777 | 6.66e-06 | 4.33e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_2B_PDonly` | effect | `harmonized_nulisa_cns_csf_IL6_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_stage_2A_vs_2B` | -1.78 | 0.394 | 6.66e-06 | 4.33e-06 | 4.63e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CHRDL1_NPX` | -0.34 | 0.0723 | 6.71e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NRP1_NPX` | +0.301 | 0.0655 | 6.76e-06 | 4.32e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_stage_2B_vs_3_PDonly` | effect | `harmonized_olink_plasma_TPX2_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_2B_vs_3` | +0.601 | 0.134 | 6.77e-06 | 4.32e-06 | 9.26e-06 | family | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_HPD_NPX` | -0.11 | 0.0241 | 6.87e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_csf_UBXN2B_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.814 | 0.179 | 6.88e-06 | 4.32e-06 | 9.23e-06 | family | EUR | OLS |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SNX33_NPX` | -0.219 | 0.0483 | 7.10e-06 | 4.32e-06 | 9.24e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RSPO3_NPX` | -0.451 | 0.0962 | 7.11e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_plasma_IGFBP2_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_early_vs_late` | +0.348 | 0.0774 | 7.13e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LOGIT |
| `NSD_vs_HC` | effect | `harmonized_olink_plasma_PPP1R18_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_HC` | -0.715 | 0.159 | 7.15e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_CALB2_NPX` | +0.162 | 0.0308 | 7.18e-06 | 4.34e-06 | 9.28e-06 | family | AJ, META_FE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_olink_plasma_CLEC10A_NPX` | `harmonized_olink_plasma` | `grp_NSD_vs_HC` | -0.45 | 0.1 | 7.19e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LOGIT |
| `trajectory_CI_PI` | within | `YEAR:CI_PI[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_PLAAT3_NPX` | -0.159 | 0.0354 | 7.23e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_LMOD1_NPX` | +0.0524 | 0.0117 | 7.26e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_COL9A1_NPX` | +0.0535 | 0.0119 | 7.38e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_ACP1_NPX` | -0.199 | 0.0441 | 7.43e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `slope_lowput_ratio` | effect | `harmonized_olink_plasma_KIAA1614_NPX` | `harmonized_olink_plasma` | `slope_lowput_ratio` | +0.0218 | 0.00461 | 7.48e-06 | 4.32e-06 | 9.26e-06 | family | AJ | OLS |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_AOC3_NPX` | -0.122 | 0.0271 | 7.52e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_PTGDS_NPX` | +0.32 | 0.0672 | 7.56e-06 | 4.33e-06 | 9.27e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | within | `YEAR:grp_HC_vs_Prodromal[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_FN1_NPX` | -0.121 | 0.0269 | 7.58e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_SOD1_NPQ` | -0.0873 | 0.0195 | 7.58e-06 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_SLITRK1_NPX` | -0.402 | 0.086 | 7.58e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_IRAK4_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | +0.756 | 0.169 | 7.63e-06 | 4.32e-06 | 9.23e-06 | family | META_FE | OLS |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_RNASE6_NPX` | +0.362 | 0.0773 | 7.63e-06 | 4.33e-06 | 9.25e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | `YEAR:grp_HCNSDneg_vs_PDNSDpos[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.126 | 0.0283 | 7.66e-06 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_csf_EDARADD_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.936 | 0.207 | 7.70e-06 | 4.32e-06 | 9.23e-06 | family | EUR | OLS |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_FLT1_NPQ` | +0.0953 | 0.0212 | 7.70e-06 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_FLT1_NPQ` | +0.0832 | 0.0185 | 7.73e-06 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CHL1_NPX` | +0.481 | 0.0975 | 7.75e-06 | 4.33e-06 | 9.26e-06 | family | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ACP3_NPX` | +0.0738 | 0.0163 | 7.85e-06 | 4.32e-06 | 9.26e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_olink_plasma_ANPEP_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_early_vs_late` | +0.35 | 0.0782 | 7.87e-06 | 4.32e-06 | 9.26e-06 | family | EUR | LOGIT |
| `NSD_stage_2B_vs_3` | effect | `harmonized_nulisa_inf_plasma_CD276_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_stage_2B_vs_3` | +0.43 | 0.0962 | 7.90e-06 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE | LOGIT |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_olink_csf_EDARADD_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.934 | 0.207 | 7.97e-06 | 4.32e-06 | 9.23e-06 | family | EUR | OLS |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_olink_csf` | `harmonized_olink_csf_DNAJC30_NPX` | -0.189 | 0.0423 | 7.97e-06 | 4.32e-06 | 9.23e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `cox_pm_any_x_PRS` | effect | `harmonized_olink_plasma_GCK_NPX` | `harmonized_olink_plasma` | `tte_pm_any_years` | +0.279 | 0.0625 | 7.97e-06 | 4.32e-06 | 9.26e-06 | family | AJ | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_LPL_NPX` | +0.148 | 0.0331 | 8.04e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_2B_vs_3` | effect | `harmonized_olink_plasma_CD276_NPX` | `harmonized_olink_plasma` | `grp_NSD_stage_2B_vs_3` | +0.442 | 0.0992 | 8.28e-06 | 4.32e-06 | 9.26e-06 | family | META_FE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_IGFBP6_NPX` | +0.331 | 0.0727 | 8.31e-06 | 4.32e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_olink_csf_DLK1_NPX` | `harmonized_olink_csf` | `grp_NSD_vs_HC` | -0.898 | 0.201 | 8.31e-06 | 4.32e-06 | 9.23e-06 | family | META_FE | LOGIT |
| `NSD_stage_2A_vs_3` | effect | `harmonized_nulisa_cns_plasma_pTau-231_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_stage_2A_vs_3` | +0.938 | 0.211 | 8.45e-06 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_MSMP_NPX` | -0.111 | 0.0249 | 8.48e-06 | 4.32e-06 | 9.26e-06 | family | META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_TNFSF12_NPX` | -0.111 | 0.0248 | 8.58e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CTHRC1_NPX` | +0.109 | 0.0244 | 8.63e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_ADAMTS13_NPX` | -0.063 | 0.0142 | 8.71e-06 | 4.32e-06 | 9.26e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_PI` | between | `YEAR:CI_PI[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_HSD17B3_NPX` | +0.335 | 0.0742 | 8.74e-06 | 4.32e-06 | 9.24e-06 | family | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_nulisa_inf_csf_IL5_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_vs_notNSD_prodromal` | +2.33 | 0.524 | 8.77e-06 | 4.32e-06 | 2.02e-04 | family | AJ | LOGIT |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_olink_csf_TSHZ2_NPX` | `harmonized_olink_csf` | `slope_updrs3_off` | -0.732 | 0.165 | 8.85e-06 | 4.32e-06 | 9.23e-06 | family | META_FE | OLS |
| `NSD_vs_notNSD_sPD` | effect | `harmonized_nulisa_cns_csf_NEFL_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_sPD` | -1.09 | 0.246 | 8.86e-06 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_olink_plasma` | `harmonized_olink_plasma_INCENP_NPX` | +0.109 | 0.0244 | 8.89e-06 | 4.32e-06 | 9.26e-06 | family | META_FE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_HSD17B3_NPX` | +0.141 | 0.0317 | 8.89e-06 | 4.32e-06 | 9.23e-06 | family | META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_CTSV_NPX` | -0.153 | 0.0344 | 8.92e-06 | 4.32e-06 | 9.23e-06 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_GPC2_NPX` | +0.0948 | 0.0213 | 8.98e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `cox_moca_lt26_x_PRS` | effect | `harmonized_olink_csf_FDXR_NPX` | `harmonized_olink_csf` | `tte_moca_lt26_years` | +1.04 | 0.234 | 9.02e-06 | 4.74e-06 | 1.14e-05 | family | AJ | COX |
| `cox_nsd_2a_to_later` | effect | `harmonized_olink_csf_L1CAM_NPX` | `harmonized_olink_csf` | `tte_nsd_2a_to_later_years` | +1.59 | 0.358 | 9.03e-06 | 4.32e-06 | 9.24e-06 | family | EUR, META_FE, META_RE | COX |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_NODAL_NPX` | +0.304 | 0.0677 | 9.11e-06 | 4.32e-06 | 9.24e-06 | family | EUR | LMM_RS_lbfgs |
| `cox_stage_d` | effect | `harmonized_olink_csf_PPME1_NPX` | `harmonized_olink_csf` | `tte_stage_d_years` | -0.345 | 0.0777 | 9.12e-06 | 4.32e-06 | 9.24e-06 | family | META_FE, META_RE | COX |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_olink_csf` | `harmonized_olink_csf_MEGF9_NPX` | +0.0753 | 0.0169 | 9.20e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LMM_RS_lbfgs |
| `NSD_stage_2B_vs_3` | effect | `harmonized_olink_csf_FST_NPX` | `harmonized_olink_csf` | `grp_NSD_stage_2B_vs_3` | +0.822 | 0.185 | 9.20e-06 | 4.32e-06 | 9.23e-06 | family | EUR | LOGIT |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_LGALS9_NPQ` | +0.136 | 0.0306 | 9.42e-06 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `cox_pm_mc_any` | effect | `harmonized_nulisa_inf_plasma_TNFSF8_NPQ` | `harmonized_nulisa_inf_plasma` | `tte_pm_mc_any_years` | +0.373 | 0.0843 | 9.74e-06 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE | COX |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_cns_csf_NPTX2_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_cog_any_years` | -0.609 | 0.139 | 1.15e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | COX |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_BD-pTau-231_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | -0.346 | 0.0789 | 1.17e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | OLS |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_FLT1_NPQ` | +0.0809 | 0.0184 | 1.19e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_inf_csf_IL23_NPQ` | `harmonized_nulisa_inf_csf` | `tte_pm_cog_any_years` | +0.382 | 0.0874 | 1.21e-05 | 4.32e-06 | 2.02e-04 | family | META_FE, META_RE | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.147 | 0.0336 | 1.24e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_2B_vs_3_PDonly` | effect | `harmonized_nulisa_cns_csf_CCL17_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_stage_2B_vs_3` | +0.634 | 0.145 | 1.25e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `sPD_vs_LRRK2_SAAadj` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_sPD_vs_LRRK2` | +2.52 | 0.576 | 1.27e-05 | 4.32e-06 | 4.17e-03 | family | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_FLT1_NPQ` | +0.0838 | 0.0192 | 1.31e-05 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_SFRP1_NPQ` | -0.243 | 0.0494 | 1.43e-05 | 4.34e-06 | 3.94e-04 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.206 | 0.047 | 1.43e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `slope_lowput_ratio` | effect | `harmonized_nulisa_cns_csf_BDNF_NPQ` | `harmonized_nulisa_cns_csf` | `slope_lowput_ratio` | -0.019 | 0.00416 | 1.46e-05 | 4.32e-06 | 3.91e-04 | family | AJ | OLS |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_nulisa_cns_plasma_DDC_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_vs_notNSD_LRRK2` | +0.754 | 0.174 | 1.48e-05 | 4.32e-06 | 3.94e-04 | family | META_FE | LOGIT |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_NPTX2_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | +0.364 | 0.0843 | 1.54e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | OLS |
| `cox_stage_d` | effect | `harmonized_nulisa_inf_csf_IFNA1; IFNA13_NPQ` | `harmonized_nulisa_inf_csf` | `tte_stage_d_years` | -0.316 | 0.0732 | 1.54e-05 | 4.32e-06 | 2.02e-04 | family | EUR | COX |
| `NSD_vs_HC` | effect | `harmonized_nulisa_cns_csf_UCHL1_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_HC` | +0.514 | 0.12 | 1.87e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LOGIT |
| `NSD_vs_HC` | effect | `asyn` | `non_highthroughput_proteomics` | `grp_NSD_vs_HC` | -0.443 | 0.104 | 1.94e-05 | 4.32e-06 | 3.57e-03 | family | EUR, META_FE | LOGIT |
| `cox_stage_d` | effect | `harmonized_nulisa_inf_csf_EPO_NPQ` | `harmonized_nulisa_inf_csf` | `tte_stage_d_years` | +0.412 | 0.0965 | 1.96e-05 | 4.32e-06 | 2.02e-04 | family | AJ, META_FE | COX |
| `cox_pm_mc_any_x_PRS` | effect | `harmonized_nulisa_cns_csf_IFNG_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_mc_any_years` | -0.419 | 0.0985 | 2.10e-05 | 4.32e-06 | 3.91e-04 | family | AJ | COX |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_FLT1_NPQ` | +0.0761 | 0.0178 | 2.12e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CD3E_NPQ` | -0.0908 | 0.0212 | 2.15e-05 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `sPD_vs_LRRK2` | effect | `total_di_22_6_BMP` | `non_highthroughput_proteomics` | `grp_sPD_vs_LRRK2` | +3.44 | 0.81 | 2.16e-05 | 4.32e-06 | 3.85e-03 | family | EUR, META_FE | LOGIT |
| `cox_stage_d_x_PRS` | effect | `harmonized_nulisa_cns_csf_PDGFRB_NPQ` | `harmonized_nulisa_cns_csf` | `tte_stage_d_years` | +0.436 | 0.103 | 2.17e-05 | 4.32e-06 | 3.91e-04 | family | EUR | COX |
| `slope_moca` | effect | `harmonized_nulisa_cns_plasma_ACHE_NPQ` | `harmonized_nulisa_cns_plasma` | `slope_moca` | -0.525 | 0.121 | 2.55e-05 | 4.32e-06 | 3.94e-04 | family | AJ | OLS |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_GDF2_NPQ` | +0.0728 | 0.0173 | 2.58e-05 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS` | effect | `harmonized_nulisa_cns_plasma_SNAP25_NPQ` | `harmonized_nulisa_cns_plasma` | `slope_updrs3_off` | -0.594 | 0.141 | 2.66e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE | OLS |
| `cox_nsd_2b_to_later` | effect | `harmonized_nulisa_inf_csf_MIF_NPQ` | `harmonized_nulisa_inf_csf` | `tte_nsd_2b_to_later_years` | -0.667 | 0.159 | 2.70e-05 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | COX |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_nulisa_cns_csf_BD-pTau-181_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_prodromal` | +1.06 | 0.255 | 3.17e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LOGIT |
| `cox_pm_any` | effect | `harmonized_nulisa_cns_csf_GDNF_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_any_years` | +0.186 | 0.0449 | 3.33e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | COX |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CXCL6_NPQ` | -0.381 | 0.0907 | 3.34e-05 | 4.32e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.108 | 0.0261 | 3.37e-05 | 4.32e-06 | 3.91e-04 | family | AJ, META_FE | LMM_RS_lbfgs |
| `sPD_vs_LRRK2` | effect | `harmonized_nulisa_cns_csf_ENO2_NPQ` | `harmonized_nulisa_cns_csf` | `grp_sPD_vs_LRRK2` | -2.09 | 0.504 | 3.44e-05 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | LOGIT |
| `NSD_stage_2B_vs_3` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_NSD_stage_2B_vs_3` | +0.299 | 0.0722 | 3.46e-05 | 4.32e-06 | 3.57e-03 | family | EUR, META_FE, META_RE | LOGIT |
| `cox_pm_mc_any` | effect | `harmonized_nulisa_inf_csf_IL11_NPQ` | `harmonized_nulisa_inf_csf` | `tte_pm_mc_any_years` | -0.336 | 0.0811 | 3.52e-05 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | COX |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_MDH1_NPQ` | +0.0666 | 0.016 | 3.54e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_MMP9_NPQ` | +0.159 | 0.0385 | 3.55e-05 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_SFTPD_NPQ` | -0.122 | 0.0294 | 3.55e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `slope_moca` | effect | `upsit` | `non_highthroughput_proteomics` | `slope_moca` | +0.161 | 0.0389 | 3.59e-05 | 4.32e-06 | 3.57e-03 | family | AJ, EUR, META_FE | OLS |
| `slope_moca` | effect | `harmonized_nulisa_inf_csf_IL23_NPQ` | `harmonized_nulisa_inf_csf` | `slope_moca` | -0.179 | 0.0434 | 3.62e-05 | 4.32e-06 | 2.02e-04 | family | META_FE | OLS |
| `cox_pm_any` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_any_years` | +0.25 | 0.0605 | 3.67e-05 | 4.32e-06 | 3.91e-04 | family | META_FE | COX |
| `cox_moca_lt26` | effect | `harmonized_nulisa_cns_csf_NPTX2_NPQ` | `harmonized_nulisa_cns_csf` | `tte_moca_lt26_years` | -0.596 | 0.144 | 3.67e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | COX |
| `sPD_vs_LRRK2_SAAadj` | effect | `harmonized_nulisa_cns_csf_ENO2_NPQ` | `harmonized_nulisa_cns_csf` | `grp_sPD_vs_LRRK2` | -2.42 | 0.59 | 3.93e-05 | 4.32e-06 | 3.91e-04 | family | META_FE | LOGIT |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_NPTXR_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | +0.442 | 0.108 | 3.97e-05 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | OLS |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_KNG1_NPQ` | +0.187 | 0.0455 | 4.16e-05 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_VEGFD_NPQ` | +0.0513 | 0.0125 | 4.17e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CD46_NPQ` | +0.112 | 0.0273 | 4.19e-05 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_TNFRSF17_NPQ` | -0.289 | 0.0683 | 4.22e-05 | 4.33e-06 | 2.02e-04 | family | AJ | LMM_RS_lbfgs |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_inf_plasma_GFAP_NPQ` | `harmonized_nulisa_inf_plasma` | `tte_pm_cog_any_years` | +0.414 | 0.101 | 4.23e-05 | 4.32e-06 | 2.01e-04 | family | EUR | COX |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_PARK7_NPQ` | +0.0579 | 0.0141 | 4.24e-05 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_ICOSLG_NPQ` | -0.109 | 0.0266 | 4.37e-05 | 4.32e-06 | 2.02e-04 | family | EUR | LMM_RS_lbfgs |
| `slope_lowput_ratio` | effect | `harmonized_nulisa_inf_plasma_IL10_NPQ` | `harmonized_nulisa_inf_plasma` | `slope_lowput_ratio` | -0.0173 | 0.00407 | 4.50e-05 | 4.32e-06 | 2.01e-04 | family | AJ | OLS |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_IL5RA_NPQ` | -0.144 | 0.0352 | 4.68e-05 | 4.32e-06 | 2.02e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_MICA_NPQ` | +0.141 | 0.0344 | 4.70e-05 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `NSD_stage_early_vs_late` | effect | `harmonized_nulisa_inf_plasma_AGRP_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_stage_early_vs_late` | +0.298 | 0.0732 | 4.74e-05 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE | LOGIT |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_nulisa_inf_csf_TGFB1_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_vs_notNSD_GBA` | +1.66 | 0.408 | 4.79e-05 | 4.32e-06 | 2.02e-04 | family | META_FE, META_RE | LOGIT |
| `NSD_stage_2A_vs_3` | effect | `harmonized_nulisa_inf_csf_TNFRSF11B_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_stage_2A_vs_3` | +0.655 | 0.161 | 4.83e-05 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_CCL19_NPQ` | -0.135 | 0.0332 | 4.85e-05 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_TAFA5_NPQ` | -0.147 | 0.0352 | 4.98e-05 | 4.33e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_IGF1R_NPQ` | +0.0835 | 0.0206 | 5.18e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_NPY_NPQ` | -0.498 | 0.11 | 5.28e-05 | 4.34e-06 | 3.94e-04 | family | AJ | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_LIF_NPQ` | -0.344 | 0.0842 | 5.31e-05 | 4.32e-06 | 2.02e-04 | family | EUR | LMM_RS_lbfgs |
| `NSD_vs_notNSD_prodromal` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `grp_NSD_vs_notNSD_prodromal` | -0.29 | 0.0719 | 5.40e-05 | 4.32e-06 | 4.17e-03 | family | EUR, META_FE | LOGIT |
| `trajectory_HC_vs_PD` | between | `YEAR:grp_HC_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CXCL6_NPQ` | -0.386 | 0.0934 | 5.50e-05 | 4.33e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `NSD_stage_2B_vs_3` | effect | `harmonized_nulisa_inf_plasma_BMP7_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_stage_2B_vs_3` | -0.377 | 0.0936 | 5.57e-05 | 4.32e-06 | 2.01e-04 | family | EUR | LOGIT |
| `slope_moca` | effect | `harmonized_nulisa_inf_plasma_GFAP_NPQ` | `harmonized_nulisa_inf_plasma` | `slope_moca` | -0.207 | 0.0513 | 5.60e-05 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE, META_RE | OLS |
| `trajectory_HC_vs_Prodromal` | within | `YEAR:grp_HC_vs_Prodromal[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_SQSTM1_NPQ` | -0.0641 | 0.0159 | 5.62e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `slope_updrs3_off_x_PRS_LEDD_adj` | effect | `harmonized_nulisa_cns_plasma_SNAP25_NPQ` | `harmonized_nulisa_cns_plasma` | `slope_updrs3_off` | -0.569 | 0.141 | 5.67e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE | OLS |
| `cox_moca_lt26` | effect | `harmonized_nulisa_cns_plasma_PDGFRB_NPQ` | `harmonized_nulisa_cns_plasma` | `tte_moca_lt26_years` | -0.626 | 0.155 | 5.69e-05 | 4.34e-06 | 3.94e-04 | family | AJ | COX |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_IL7R_NPQ` | -0.0598 | 0.0149 | 5.77e-05 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_POSTN_NPQ` | +0.053 | 0.0132 | 5.99e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_nulisa_cns_csf_pTau-181_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_prodromal` | +0.858 | 0.214 | 6.10e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_POSTN_NPQ` | +0.0532 | 0.0133 | 6.10e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_PI` | within | `YEAR:CI_PI[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_IL34_NPQ` | +0.098 | 0.0242 | 6.13e-05 | 4.32e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_pTau-181_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | -0.385 | 0.0961 | 6.16e-05 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | OLS |
| `cox_pm_any_x_PRS` | effect | `harmonized_nulisa_cns_plasma_IL15_NPQ` | `harmonized_nulisa_cns_plasma` | `tte_pm_any_years` | +0.194 | 0.0485 | 6.28e-05 | 4.32e-06 | 3.94e-04 | family | AJ, META_FE | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_SELP_NPQ` | +0.122 | 0.0305 | 6.30e-05 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_nulisa_inf_plasma_TREM1_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_vs_HC` | -0.547 | 0.137 | 6.50e-05 | 4.32e-06 | 2.01e-04 | family | EUR | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_IL1B_NPQ` | +0.0761 | 0.019 | 6.53e-05 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_pm_cog_any_x_PRS` | effect | `harmonized_nulisa_inf_csf_SELP_NPQ` | `harmonized_nulisa_inf_csf` | `tte_pm_cog_any_years` | +0.42 | 0.105 | 6.59e-05 | 4.33e-06 | 2.02e-04 | family | AJ | COX |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NPTX2_NPQ` | +0.0394 | 0.00988 | 6.74e-05 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_AGRN_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | +0.38 | 0.0954 | 6.82e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | OLS |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_cog_any_years` | +0.344 | 0.0864 | 6.91e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | COX |
| `NSD_stage_2B_vs_3_PDonly` | effect | `harmonized_nulisa_inf_csf_CHI3L1_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_stage_2B_vs_3` | -0.772 | 0.194 | 7.15e-05 | 4.32e-06 | 2.02e-04 | family | EUR | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_IL33_NPQ` | +0.103 | 0.0258 | 7.33e-05 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_stage_2A_vs_2B` | between | `YEAR:grp_NSD_stage_2A_vs_2B[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_TNFSF15_NPQ` | +0.491 | 0.112 | 7.38e-05 | 4.34e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_nulisa_inf_csf_TLR3_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_vs_HC` | +0.467 | 0.118 | 7.40e-05 | 4.32e-06 | 2.02e-04 | family | META_FE | LOGIT |
| `NSD_stage_2A_vs_2B_PDonly` | effect | `harmonized_nulisa_inf_csf_TNFRSF9_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_stage_2A_vs_2B` | -1.61 | 0.406 | 7.45e-05 | 4.33e-06 | 2.20e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_cns_csf_BD-pTau-231_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_cog_any_years` | +0.631 | 0.159 | 7.47e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | COX |
| `NSD_vs_HC` | effect | `harmonized_nulisa_inf_plasma_IL23_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_vs_HC` | +0.411 | 0.104 | 7.53e-05 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | LOGIT |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_nulisa_inf_plasma_GFAP_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_vs_notNSD_prodromal` | +0.399 | 0.101 | 8.00e-05 | 4.32e-06 | 2.01e-04 | family | EUR | LOGIT |
| `NSD_stage_2A_vs_3` | effect | `harmonized_nulisa_cns_csf_AÎ²38_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_stage_2A_vs_3` | -0.817 | 0.207 | 8.13e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `cox_stage_d_x_PRS` | effect | `harmonized_nulisa_inf_csf_CSF1_NPQ` | `harmonized_nulisa_inf_csf` | `tte_stage_d_years` | +0.596 | 0.151 | 8.22e-05 | 4.32e-06 | 2.02e-04 | family | EUR | COX |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_TNFSF9_NPQ` | -0.347 | 0.0872 | 8.32e-05 | 4.32e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `slope_lowput_ratio_x_PRS` | effect | `harmonized_nulisa_cns_plasma_UBB_NPQ` | `harmonized_nulisa_cns_plasma` | `slope_lowput_ratio` | +0.0102 | 0.00249 | 8.66e-05 | 4.32e-06 | 3.94e-04 | family | AJ | OLS |
| `cox_pm_any` | effect | `harmonized_nulisa_inf_plasma_CCL1_NPQ` | `harmonized_nulisa_inf_plasma` | `tte_pm_any_years` | -0.269 | 0.0685 | 8.71e-05 | 4.32e-06 | 2.01e-04 | family | EUR | COX |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_UCHL1_NPQ` | -0.114 | 0.0291 | 8.87e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `NSD_stage_2B_vs_3_PDonly` | effect | `harmonized_nulisa_inf_plasma_LGALS9_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_stage_2B_vs_3` | -0.658 | 0.168 | 8.94e-05 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | LOGIT |
| `cox_pm_mc_any` | effect | `harmonized_nulisa_cns_csf_NRGN_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_mc_any_years` | +1.02 | 0.261 | 9.00e-05 | 4.32e-06 | 3.91e-04 | family | AJ | COX |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_BD-pTau-217_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | -0.392 | 0.1 | 9.16e-05 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | OLS |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_VSNL1_NPQ` | +0.0764 | 0.0195 | 9.26e-05 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `NSD_stage_2B_vs_3_PDonly` | effect | `harmonized_nulisa_inf_plasma_CX3CL1_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_stage_2B_vs_3` | +0.612 | 0.157 | 9.36e-05 | 4.32e-06 | 2.01e-04 | family | EUR | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_pTDP43-409_NPQ` | +0.0763 | 0.0195 | 9.59e-05 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_VEGFD_NPQ` | +0.147 | 0.0377 | 9.61e-05 | 4.32e-06 | 2.03e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `cox_stage_d` | effect | `harmonized_nulisa_inf_plasma_CSF3_NPQ` | `harmonized_nulisa_inf_plasma` | `tte_stage_d_years` | +0.309 | 0.0798 | 1.06e-04 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | COX |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_POSTN_NPQ` | +0.0496 | 0.0128 | 1.08e-04 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_VSNL1_NPQ` | +0.0701 | 0.0181 | 1.09e-04 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `sPD_vs_LRRK2` | effect | `total_di_18_1_BMP` | `non_highthroughput_proteomics` | `grp_sPD_vs_LRRK2` | +2.19 | 0.565 | 1.09e-04 | 4.32e-06 | 3.85e-03 | family | EUR, META_FE | LOGIT |
| `cox_nsd_2a_to_later_x_PRS` | effect | `harmonized_nulisa_inf_csf_TNFRSF9_NPQ` | `harmonized_nulisa_inf_csf` | `tte_nsd_2a_to_later_years` | +0.783 | 0.202 | 1.09e-04 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | COX |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NPY_NPQ` | -0.0544 | 0.014 | 1.12e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_IL17RA_NPQ` | -0.0721 | 0.0185 | 1.14e-04 | 4.32e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NPY_NPQ` | -0.0536 | 0.0138 | 1.14e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_cogstate_worsen_x_PRS` | effect | `harmonized_nulisa_inf_csf_CCL8_NPQ` | `harmonized_nulisa_inf_csf` | `tte_cogstate_worsen_years` | +0.463 | 0.12 | 1.14e-04 | 4.32e-06 | 2.02e-04 | family | AJ | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_SOD1_NPQ` | +0.084 | 0.0218 | 1.15e-04 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_FURIN_NPQ` | -0.164 | 0.0427 | 1.15e-04 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `slope_moca` | effect | `harmonized_nulisa_inf_plasma_CD80_NPQ` | `harmonized_nulisa_inf_plasma` | `slope_moca` | -0.561 | 0.142 | 1.15e-04 | 4.32e-06 | 2.01e-04 | family | AJ | OLS |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_PGK1_NPQ` | +0.0594 | 0.0154 | 1.16e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_pm_cog_any_x_PRS` | effect | `harmonized_nulisa_inf_csf_BST2_NPQ` | `harmonized_nulisa_inf_csf` | `tte_pm_cog_any_years` | +0.468 | 0.121 | 1.16e-04 | 4.33e-06 | 2.02e-04 | family | AJ | COX |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_CD200_NPQ` | +0.154 | 0.0398 | 1.16e-04 | 4.32e-06 | 2.03e-04 | family | EUR | LMM_RS_lbfgs |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_nulisa_inf_csf_IL19_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_vs_notNSD_LRRK2` | -0.938 | 0.244 | 1.18e-04 | 4.32e-06 | 2.02e-04 | family | AJ, META_FE, META_RE | LOGIT |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_nulisa_inf_csf_TAFA5_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_vs_notNSD_LRRK2` | -2.62 | 0.681 | 1.22e-04 | 4.32e-06 | 2.02e-04 | family | AJ | LOGIT |
| `trajectory_stage_2A_vs_2B` | within | `YEAR:grp_NSD_stage_2A_vs_2B[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_IL6_NPQ` | +0.309 | 0.0805 | 1.22e-04 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_GBA` | effect | `harmonized_nulisa_cns_plasma_DDC_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_vs_notNSD_GBA` | +0.883 | 0.23 | 1.27e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LOGIT |
| `slope_lowput_ratio` | effect | `harmonized_nulisa_inf_csf_IL4_NPQ` | `harmonized_nulisa_inf_csf` | `slope_lowput_ratio` | +0.00789 | 0.00205 | 1.29e-04 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | OLS |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_VEGFD_NPQ` | +0.0525 | 0.0137 | 1.29e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[within]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_TNFSF12_NPQ` | -0.181 | 0.046 | 1.33e-04 | 4.33e-06 | 2.02e-04 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_DDC_NPQ` | +0.109 | 0.0285 | 1.39e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_stage_2A_vs_3_PDonly` | effect | `harmonized_nulisa_inf_csf_TNFRSF11B_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_stage_2A_vs_3` | +1.33 | 0.35 | 1.42e-04 | 4.33e-06 | 2.11e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_CD63_NPQ` | -0.0992 | 0.026 | 1.43e-04 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE | LMM_RS_lbfgs |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_IL15_NPQ` | +0.0855 | 0.0225 | 1.43e-04 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_SELP_NPQ` | +0.0805 | 0.0211 | 1.44e-04 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_nulisa_cns_plasma_IL17A_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_vs_notNSD_LRRK2` | +0.639 | 0.168 | 1.44e-04 | 4.32e-06 | 3.94e-04 | family | AJ, META_FE, META_RE | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_IKBKG_NPQ` | +0.0822 | 0.0216 | 1.46e-04 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_pm_mc_any_x_PRS` | effect | `harmonized_nulisa_inf_csf_CEACAM5_NPQ` | `harmonized_nulisa_inf_csf` | `tte_pm_mc_any_years` | -0.497 | 0.131 | 1.47e-04 | 4.32e-06 | 2.02e-04 | family | AJ | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_PARK7_NPQ` | +0.119 | 0.0314 | 1.49e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_nulisa_cns_plasma_TREM1_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_vs_HC` | -0.481 | 0.127 | 1.50e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_TARDBP_NPQ` | +0.0605 | 0.0159 | 1.52e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | within | `YEAR:grp_HC_vs_Prodromal[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_VEGFD_NPQ` | +0.0481 | 0.0127 | 1.52e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_LRRK2` | effect | `harmonized_nulisa_cns_csf_IL6_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_LRRK2` | -1.04 | 0.273 | 1.53e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | LOGIT |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_S100A12_NPQ` | +0.0826 | 0.0218 | 1.53e-04 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `cox_cogstate_worsen_x_PRS` | effect | `harmonized_nulisa_cns_plasma_BD-pTau-231_NPQ` | `harmonized_nulisa_cns_plasma` | `tte_cogstate_worsen_years` | +0.318 | 0.084 | 1.56e-04 | 4.32e-06 | 3.94e-04 | family | AJ | COX |
| `cox_nsd_2b_to_later` | effect | `harmonized_nulisa_inf_csf_MMP8_NPQ` | `harmonized_nulisa_inf_csf` | `tte_nsd_2b_to_later_years` | +0.377 | 0.0997 | 1.58e-04 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | COX |
| `NSD_vs_HC` | effect | `harmonized_nulisa_inf_csf_VEGFC_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_vs_HC` | +0.549 | 0.146 | 1.66e-04 | 4.32e-06 | 2.02e-04 | family | EUR | LOGIT |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_IGF1R_NPQ` | +0.118 | 0.0311 | 1.68e-04 | 4.32e-06 | 3.91e-04 | family | EUR | LMM_RS_lbfgs |
| `slope_moca` | effect | `harmonized_nulisa_inf_csf_TNFRSF18_NPQ` | `harmonized_nulisa_inf_csf` | `slope_moca` | +0.256 | 0.0681 | 1.69e-04 | 4.32e-06 | 2.02e-04 | family | META_FE, META_RE | OLS |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_IL1RN_NPQ` | +0.112 | 0.0299 | 1.73e-04 | 4.32e-06 | 2.01e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_PTX3_NPQ` | +0.0479 | 0.0128 | 1.77e-04 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_VEGFD_NPQ` | +0.0495 | 0.0132 | 1.78e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_sPD` | effect | `harmonized_nulisa_inf_csf_TNFSF13_NPQ` | `harmonized_nulisa_inf_csf` | `grp_NSD_vs_notNSD_sPD` | +0.949 | 0.253 | 1.79e-04 | 4.32e-06 | 2.02e-04 | family | EUR, META_FE, META_RE | LOGIT |
| `NSD_stage_2A_vs_3` | effect | `harmonized_nulisa_inf_plasma_CD276_NPQ` | `harmonized_nulisa_inf_plasma` | `grp_NSD_stage_2A_vs_3` | +0.42 | 0.112 | 1.80e-04 | 4.32e-06 | 2.01e-04 | family | META_FE | LOGIT |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_MMP3_NPQ` | +0.0407 | 0.0109 | 1.80e-04 | 4.32e-06 | 2.01e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | between | `YEAR:CI_MOCA[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_TNFSF12_NPQ` | -0.162 | 0.0429 | 1.82e-04 | 4.33e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_nulisa_cns_csf_AGRN_NPQ` | `harmonized_nulisa_cns_csf` | `grp_NSD_vs_notNSD_prodromal` | -0.675 | 0.18 | 1.82e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | LOGIT |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CD3E_NPQ` | -0.0505 | 0.0135 | 1.84e-04 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CD40_NPQ` | +0.0557 | 0.0149 | 1.87e-04 | 4.32e-06 | 2.01e-04 | family | EUR | LMM_RS_lbfgs |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_pTau-231_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | -0.317 | 0.085 | 1.92e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | OLS |
| `trajectory_NSD_vs_HC` | between | `YEAR:grp_NSD_vs_HC[between]` | `harmonized_nulisa_inf_plasma` | `harmonized_nulisa_inf_plasma_CXCL6_NPQ` | -0.362 | 0.095 | 1.93e-04 | 4.33e-06 | 2.01e-04 | family | AJ | LMM_RS_lbfgs |
| `cox_moca_lt26` | effect | `harmonized_nulisa_cns_csf_pTau-217_NPQ` | `harmonized_nulisa_cns_csf` | `tte_moca_lt26_years` | +0.604 | 0.162 | 1.94e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | COX |
| `cox_pm_any` | effect | `harmonized_nulisa_inf_plasma_CCL20_NPQ` | `harmonized_nulisa_inf_plasma` | `tte_pm_any_years` | +0.195 | 0.0524 | 2.01e-04 | 4.32e-06 | 2.01e-04 | family | EUR | COX |
| `sPD_vs_GBA_SAAadj` | effect | `harmonized_nulisa_cns_csf_DDC_NPQ` | `harmonized_nulisa_cns_csf` | `grp_sPD_vs_GBA` | +1.02 | 0.274 | 2.02e-04 | 4.32e-06 | 3.91e-04 | family | META_FE | LOGIT |
| `slope_moca` | effect | `harmonized_nulisa_cns_plasma_BD-pTau-217_NPQ` | `harmonized_nulisa_cns_plasma` | `slope_moca` | -0.569 | 0.149 | 2.02e-04 | 4.32e-06 | 3.94e-04 | family | AJ | OLS |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_nulisa_inf_csf` | `harmonized_nulisa_inf_csf_VEGFD_NPQ` | +0.188 | 0.0506 | 2.02e-04 | 4.32e-06 | 2.02e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_BD-pTau-181_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | -0.408 | 0.11 | 2.06e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | OLS |
| `slope_moca` | effect | `harmonized_nulisa_cns_csf_pTau-217_NPQ` | `harmonized_nulisa_cns_csf` | `slope_moca` | -0.344 | 0.0926 | 2.08e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | OLS |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_pSNCA-129_NPQ` | -0.0689 | 0.0186 | 2.10e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `NSD_vs_notNSD_prodromal` | effect | `harmonized_nulisa_cns_plasma_TARDBP_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_vs_notNSD_prodromal` | -0.898 | 0.243 | 2.18e-04 | 4.32e-06 | 3.94e-04 | family | META_FE | LOGIT |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_PGK1_NPQ` | -0.0656 | 0.0177 | 2.20e-04 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_NPY_NPQ` | -0.0531 | 0.0143 | 2.24e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_PARK7_NPQ` | +0.0997 | 0.0268 | 2.27e-04 | 4.33e-06 | 3.94e-04 | family | AJ | LMM_RS_lbfgs |
| `sPD_vs_LRRK2_SAAadj` | effect | `total_di_22_6_BMP` | `non_highthroughput_proteomics` | `grp_sPD_vs_LRRK2` | +4.44 | 1.2 | 2.29e-04 | 4.32e-06 | 4.17e-03 | family | EUR, META_FE | LOGIT |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_DDC_NPQ` | +0.0978 | 0.0264 | 2.29e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_moca_lt26` | effect | `harmonized_nulisa_cns_csf_pTau-181_NPQ` | `harmonized_nulisa_cns_csf` | `tte_moca_lt26_years` | +0.619 | 0.168 | 2.35e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | COX |
| `trajectory_CI_MOCA` | within | `YEAR:CI_MOCA[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_SNAP25_NPQ` | -0.0666 | 0.0181 | 2.35e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_SOD1_NPQ` | -0.0933 | 0.0254 | 2.42e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LMM_RS_lbfgs |
| `cox_stage_d_x_PRS` | effect | `harmonized_nulisa_cns_plasma_IL5_NPQ` | `harmonized_nulisa_cns_plasma` | `tte_stage_d_years` | -0.347 | 0.0953 | 2.65e-04 | 4.33e-06 | 3.94e-04 | family | AJ | COX |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_cns_csf_pTau-181_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_cog_any_years` | +0.731 | 0.2 | 2.66e-04 | 4.32e-06 | 3.91e-04 | family | EUR | COX |
| `cox_pm_mc_any_x_PRS` | effect | `harmonized_nulisa_cns_csf_PDGFRB_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_mc_any_years` | +0.275 | 0.0755 | 2.70e-04 | 4.32e-06 | 3.91e-04 | family | META_FE, META_RE | COX |
| `trajectory_NSD_vs_HC` | within | `YEAR:grp_NSD_vs_HC[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_VSNL1_NPQ` | +0.0712 | 0.0195 | 2.76e-04 | 4.32e-06 | 3.91e-04 | family | EUR, META_FE, META_RE | LMM_RS_lbfgs |
| `trajectory_HC_vs_PD` | within | `YEAR:grp_HC_vs_PD[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_PGK1_NPQ` | +0.085 | 0.0234 | 2.76e-04 | 4.32e-06 | 3.91e-04 | family | META_FE | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | `YEAR:grp_HCNSDneg_vs_PDNSDpos[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_FABP3_NPQ` | +0.0761 | 0.0209 | 2.92e-04 | 4.32e-06 | 3.91e-04 | family | EUR | LMM_RS_lbfgs |
| `NSD_vs_HC` | effect | `harmonized_nulisa_cns_plasma_GDI1_NPQ` | `harmonized_nulisa_cns_plasma` | `grp_NSD_vs_HC` | +0.677 | 0.187 | 2.92e-04 | 4.32e-06 | 3.94e-04 | family | META_FE, META_RE | LOGIT |
| `cox_nsd_2a_to_later_x_PRS` | effect | `harmonized_nulisa_cns_plasma_APOE_NPQ` | `harmonized_nulisa_cns_plasma` | `tte_nsd_2a_to_later_years` | +0.606 | 0.167 | 2.93e-04 | 4.32e-06 | 3.94e-04 | family | EUR, META_FE, META_RE | COX |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_IL1B_NPQ` | +0.0897 | 0.0248 | 2.96e-04 | 4.32e-06 | 3.94e-04 | family | META_FE | LMM_RS_lbfgs |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | `YEAR:grp_ProdNSDpos_vs_PDNSDpos[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_BDNF_NPQ` | -0.171 | 0.0475 | 3.26e-04 | 4.32e-06 | 3.91e-04 | family | META_FE | LMM_RS_lbfgs |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_cns_csf_pTau-231_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_cog_any_years` | +0.632 | 0.176 | 3.32e-04 | 4.32e-06 | 3.91e-04 | family | EUR | COX |
| `cox_moca_lt26_x_PRS` | effect | `harmonized_nulisa_cns_plasma_VGF_NPQ` | `harmonized_nulisa_cns_plasma` | `tte_moca_lt26_years` | +0.462 | 0.129 | 3.34e-04 | 4.74e-06 | 3.94e-04 | family | AJ | COX |
| `cox_stage_d` | effect | `harmonized_nulisa_cns_plasma_CD40LG_NPQ` | `harmonized_nulisa_cns_plasma` | `tte_stage_d_years` | -0.619 | 0.173 | 3.37e-04 | 4.32e-06 | 3.94e-04 | family | EUR | COX |
| `trajectory_Prodromal_vs_PD` | within | `YEAR:grp_Prodromal_vs_PD[within]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_VGF_NPQ` | -0.0654 | 0.0181 | 3.44e-04 | 4.32e-06 | 3.94e-04 | family | AJ | LMM_RS_lbfgs |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_S100B_NPQ` | +0.184 | 0.0511 | 3.45e-04 | 4.32e-06 | 3.91e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_stage_d_x_PRS` | effect | `harmonized_nulisa_cns_csf_IL7_NPQ` | `harmonized_nulisa_cns_csf` | `tte_stage_d_years` | +0.412 | 0.115 | 3.49e-04 | 4.32e-06 | 3.91e-04 | family | EUR | COX |
| `cox_pm_cog_any_x_PRS` | effect | `harmonized_nulisa_cns_csf_APOE4_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_cog_any_years` | +0.333 | 0.0933 | 3.50e-04 | 4.33e-06 | 3.91e-04 | family | AJ | COX |
| `trajectory_Prodromal_vs_PD` | between | `YEAR:grp_Prodromal_vs_PD[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_BDNF_NPQ` | -0.129 | 0.036 | 3.58e-04 | 4.32e-06 | 3.91e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_pm_cog_any` | effect | `harmonized_nulisa_cns_csf_ARSA_NPQ` | `harmonized_nulisa_cns_csf` | `tte_pm_cog_any_years` | +0.312 | 0.0878 | 3.77e-04 | 4.32e-06 | 3.91e-04 | family | EUR | COX |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | `YEAR:grp_HCNSDneg_vs_ProdNSDpos[between]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_IL5_NPQ` | +0.116 | 0.0325 | 3.78e-04 | 4.32e-06 | 3.91e-04 | family | META_FE | LMM_RS_lbfgs |
| `trajectory_HC_vs_Prodromal` | between | `YEAR:grp_HC_vs_Prodromal[between]` | `harmonized_nulisa_cns_plasma` | `harmonized_nulisa_cns_plasma_S100A12_NPQ` | -0.126 | 0.0353 | 3.82e-04 | 4.32e-06 | 3.94e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_cogstate_worsen_x_PRS` | effect | `harmonized_nulisa_cns_csf_pSNCA-129_NPQ` | `harmonized_nulisa_cns_csf` | `tte_cogstate_worsen_years` | +0.326 | 0.0917 | 3.84e-04 | 4.32e-06 | 3.91e-04 | family | AJ | COX |
| `trajectory_CI_PI` | within | `YEAR:CI_PI[within]` | `harmonized_nulisa_cns_csf` | `harmonized_nulisa_cns_csf_CCL2_NPQ` | -0.157 | 0.0442 | 3.90e-04 | 4.32e-06 | 3.91e-04 | family | EUR | LMM_RS_lbfgs |
| `cox_stage_d` | effect | `p9005_Genetic_PRS_PRS157` | `non_highthroughput_proteomics` | `tte_stage_d_years` | +0.232 | 0.0667 | 4.93e-04 | 4.32e-06 | 4.17e-03 | family | AJ, META_FE | COX |
| `sPD_vs_LRRK2_SAAadj` | effect | `total_di_18_1_BMP` | `non_highthroughput_proteomics` | `grp_sPD_vs_LRRK2` | +1.99 | 0.576 | 5.54e-04 | 4.32e-06 | 4.17e-03 | family | EUR, META_FE | LOGIT |
| `cox_pm_cog_any` | effect | `abeta` | `non_highthroughput_proteomics` | `tte_pm_cog_any_years` | -0.255 | 0.0746 | 6.39e-04 | 4.32e-06 | 3.57e-03 | family | META_FE, META_RE | COX |
| `sPD_vs_LRRK2` | effect | `upsit` | `non_highthroughput_proteomics` | `grp_sPD_vs_LRRK2` | +1.15 | 0.341 | 7.59e-04 | 4.32e-06 | 3.85e-03 | family | EUR, META_FE | LOGIT |
| `NSD_vs_notNSD_LRRK2` | effect | `abeta` | `non_highthroughput_proteomics` | `grp_NSD_vs_notNSD_LRRK2` | -0.983 | 0.302 | 1.15e-03 | 4.32e-06 | 3.85e-03 | family | EUR, META_FE, META_RE | LOGIT |
| `cox_stage_d` | effect | `asyn` | `non_highthroughput_proteomics` | `tte_stage_d_years` | -1.1 | 0.339 | 1.21e-03 | 4.32e-06 | 4.17e-03 | family | EUR, META_FE, META_RE | COX |
| `NSD_vs_notNSD_LRRK2` | effect | `tau` | `non_highthroughput_proteomics` | `grp_NSD_vs_notNSD_LRRK2` | -1.24 | 0.394 | 1.61e-03 | 4.32e-06 | 3.85e-03 | family | EUR, META_FE, META_RE | LOGIT |
| `cox_stage_d` | effect | `abeta` | `non_highthroughput_proteomics` | `tte_stage_d_years` | -0.702 | 0.223 | 1.63e-03 | 4.32e-06 | 4.17e-03 | family | EUR | COX |
| `NSD_stage_early_vs_late` | effect | `upsit` | `non_highthroughput_proteomics` | `grp_NSD_stage_early_vs_late` | -0.543 | 0.174 | 1.83e-03 | 4.32e-06 | 4.17e-03 | family | AJ | LOGIT |
| `NSD_stage_early_vs_late` | effect | `urate` | `non_highthroughput_proteomics` | `grp_NSD_stage_early_vs_late` | -0.652 | 0.211 | 1.99e-03 | 4.32e-06 | 4.17e-03 | family | AJ | LOGIT |
| `cox_stage_d` | effect | `tau` | `non_highthroughput_proteomics` | `tte_stage_d_years` | -0.85 | 0.275 | 2.01e-03 | 4.32e-06 | 4.17e-03 | family | EUR | COX |
| `cox_moca_lt26` | effect | `abeta` | `non_highthroughput_proteomics` | `tte_moca_lt26_years` | -0.22 | 0.073 | 2.62e-03 | 4.32e-06 | 3.57e-03 | family | META_FE, META_RE | COX |
| `cox_pm_any` | effect | `asyn` | `non_highthroughput_proteomics` | `tte_pm_any_years` | -0.213 | 0.0709 | 2.65e-03 | 4.32e-06 | 3.57e-03 | family | EUR | COX |
| `cox_moca_lt26_x_PRS` | effect | `hemohi` | `non_highthroughput_proteomics` | `tte_moca_lt26_years` | +0.965 | 0.325 | 2.99e-03 | 4.74e-06 | 4.17e-03 | family | AJ | COX |
| `cox_stage_d` | effect | `ptau` | `non_highthroughput_proteomics` | `tte_stage_d_years` | -1.01 | 0.348 | 3.83e-03 | 4.32e-06 | 4.17e-03 | family | EUR | COX |
