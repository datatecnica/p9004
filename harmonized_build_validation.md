# Harmonized blocks — validation of the built dataset

Reads the shipped `harmonized_*` columns in `Project_9004_Unified_Emerging_Biomarkers.tab`, so this tests what a downstream user consumes rather than a reconstruction.

## 1. Fill-rule integrity

Where the reference project has a value, the harmonized column must equal it exactly. The corrected mapped value may appear **only** where the reference is null. A violation would mean real reference data was silently overwritten.

| block | analytes checked | cells with reference value | **mismatches** | cells filled from mapped |
|---|---|---|---|---|
| `nulisa_cns_csf` | 60 | 136,920 | **0** | 71,820 |
| `nulisa_cns_plasma` | 60 | 159,900 | **0** | 97,080 |
| `nulisa_inf_csf` | 60 | 113,940 | **0** | 74,040 |
| `nulisa_inf_plasma` | 60 | 160,140 | **0** | 96,960 |
| `olink_csf` | 60 | 123,909 | **0** | 68,172 |
| `olink_plasma` | 60 | 174,098 | **0** | 43,673 |

**PASS** — every reference value is preserved exactly.

### 1b. Filled values match the recorded correction

Check 1 proves reference values survived. This proves the *filled* cells are exactly `applied_slope * mapped + applied_intercept`, i.e. the export matches the slope and intercept recorded per analyte in the data dictionary. Any drift here would mean the dictionary documents a correction the data does not carry.

| block | analytes checked | filled cells | **mismatches** | max abs deviation |
|---|---|---|---|---|
| `nulisa_cns_csf` | 40 | 47,880 | **0** | 1.07e-14 |
| `nulisa_cns_plasma` | 40 | 64,720 | **0** | 7.11e-15 |
| `nulisa_inf_csf` | 40 | 49,360 | **0** | 7.11e-15 |
| `nulisa_inf_plasma` | 40 | 64,640 | **0** | 1.42e-14 |
| `olink_csf` | 40 | 45,456 | **0** | 3.55e-15 |
| `olink_plasma` | 40 | 29,102 | **0** | 1.78e-15 |

**PASS** — every filled cell reproduces the recorded slope and intercept, so the export matches the dictionary.

## 2. Coverage and correction dependence

| block | reference | mapped | rows, reference alone | rows, harmonized | gain | % of rows corrected |
|---|---|---|---|---|---|---|
| `nulisa_cns_csf` | p282_CNS_CSF | p312_Neuro_CSF | 2,282 | **3,479** | +1,197 | 34.4% |
| `nulisa_cns_plasma` | p288_CNS_plasma | p312_Neuro_Plasma | 2,665 | **4,283** | +1,618 | 37.8% |
| `nulisa_inf_csf` | p282_Inflammation_CSF | p312_Inflammation_CSF | 1,899 | **3,133** | +1,234 | 39.4% |
| `nulisa_inf_plasma` | p288_Inflammation_plasma | p312_Inflammation_Plasma | 2,669 | **4,285** | +1,616 | 37.7% |
| `olink_csf` | p314_CSF | p277_CSF | 2,098 | **3,235** | +1,137 | 35.1% |
| `olink_plasma` | p314_Plasma | p293_olink_plasma | 2,942 | **3,670** | +728 | 19.8% |

## 3. Positive and negative controls on the built columns

Contrast `grp_NSD_vs_HC`, OLS with PATNO-clustered SEs over all visits, outcome z-scored. A random-intercept LMM is not identifiable — these panels average ~1.2 observations per participant.

### Positive controls

| block | analyte | ref n | ref beta | ref P | harm n | harm beta | harm P | P better? |
|---|---|---|---|---|---|---|---|---|
| `nulisa_cns_csf` | CHI3L1 | 1774 | -0.179 | 2.5e-02 | 2382 | -0.104 | 7.6e-02 | no |
| `nulisa_cns_csf` | CST3 | 1774 | -0.048 | 5.7e-01 | 2382 | -0.016 | 7.9e-01 | no |
| `nulisa_cns_csf` | DDC | 1774 | +1.307 | 4.7e-74 | 2382 | +1.235 | 1.8e-146 | **yes** |
| `nulisa_cns_csf` | GFAP | 1774 | +0.048 | 5.6e-01 | 2382 | +0.031 | 6.3e-01 | no |
| `nulisa_cns_csf` | MAPT | 1774 | -0.306 | 3.0e-04 | 2382 | -0.193 | 6.1e-03 | no |
| `nulisa_cns_csf` | NEFL | 1774 | +0.130 | 6.8e-02 | 2382 | +0.171 | 1.3e-03 | **yes** |
| `nulisa_cns_csf` | SNCA | 1774 | -0.304 | 5.1e-04 | 2382 | -0.123 | 5.5e-02 | no |
| `nulisa_cns_plasma` | CHI3L1 | 2120 | -0.175 | 4.6e-02 | 2744 | -0.225 | 1.2e-03 | **yes** |
| `nulisa_cns_plasma` | CST3 | 2120 | +0.052 | 5.0e-01 | 2744 | +0.056 | 3.5e-01 | **yes** |
| `nulisa_cns_plasma` | DDC | 2120 | +0.373 | 7.9e-06 | 2744 | +0.351 | 1.1e-08 | **yes** |
| `nulisa_cns_plasma` | GFAP | 2120 | +0.063 | 4.4e-01 | 2744 | +0.031 | 6.0e-01 | no |
| `nulisa_cns_plasma` | MAPT | 2120 | +0.014 | 9.1e-01 | 2744 | -0.069 | 3.8e-01 | **yes** |
| `nulisa_cns_plasma` | NEFL | 2120 | +0.220 | 1.2e-03 | 2744 | +0.205 | 2.4e-05 | **yes** |
| `nulisa_cns_plasma` | SNCA | 2120 | +0.031 | 7.3e-01 | 2744 | +0.030 | 6.6e-01 | **yes** |
| `nulisa_inf_csf` | CHI3L1 | 1462 | -0.201 | 2.4e-02 | 2090 | -0.082 | 2.1e-01 | no |
| `nulisa_inf_csf` | GFAP | 1462 | +0.052 | 5.7e-01 | 2090 | +0.037 | 5.8e-01 | no |
| `nulisa_inf_plasma` | CHI3L1 | 2119 | -0.179 | 5.0e-02 | 2746 | -0.223 | 1.4e-03 | **yes** |
| `nulisa_inf_plasma` | GFAP | 2119 | +0.029 | 7.1e-01 | 2746 | +0.017 | 7.7e-01 | no |
| `olink_csf` | CHI3L1 | 1156 | +0.126 | 1.3e-01 | 2271 | -0.073 | 2.6e-01 | no |
| `olink_csf` | CST3 | 1100 | +0.206 | 3.4e-02 | 2218 | +0.170 | 1.9e-02 | **yes** |
| `olink_csf` | DDC | 1199 | +1.175 | 2.2e-54 | 2312 | +1.004 | 9.1e-77 | **yes** |
| `olink_csf` | GFAP | 1200 | +0.028 | 7.1e-01 | 2314 | -0.120 | 2.3e-02 | **yes** |
| `olink_csf` | MAPT | 1200 | +0.023 | 8.1e-01 | 2314 | -0.266 | 2.7e-04 | **yes** |
| `olink_csf` | NEFL | 1196 | +0.201 | 8.9e-03 | 2310 | +0.016 | 7.7e-01 | no |
| `olink_csf` | SNCA | 1199 | +0.075 | 3.7e-01 | 2312 | -0.080 | 1.2e-01 | **yes** |
| `olink_plasma` | CHI3L1 | 1626 | -0.138 | 1.3e-01 | 2330 | -0.175 | 9.6e-03 | **yes** |
| `olink_plasma` | CST3 | 1626 | +0.102 | 1.8e-01 | 2330 | -0.014 | 8.3e-01 | no |
| `olink_plasma` | DDC | 1624 | +0.545 | 1.4e-15 | 2326 | +0.575 | 4.8e-27 | **yes** |
| `olink_plasma` | GFAP | 1624 | +0.083 | 2.6e-01 | 2326 | +0.034 | 5.9e-01 | no |
| `olink_plasma` | MAPT | 1579 | -0.201 | 4.4e-03 | 2283 | -0.120 | 2.6e-02 | no |
| `olink_plasma` | NEFL | 1623 | +0.267 | 3.4e-05 | 2325 | +0.189 | 4.3e-04 | no |
| `olink_plasma` | SNCA | 1624 | -0.039 | 6.2e-01 | 2326 | -0.084 | 2.4e-01 | **yes** |

### Negative controls

| block | analyte | ref n | ref beta | ref P | harm n | harm beta | harm P | P better? |
|---|---|---|---|---|---|---|---|---|
| `olink_csf` | CTRL | 1196 | -0.010 | 9.0e-01 | 2309 | -0.180 | 6.5e-04 | **yes** |
| `olink_plasma` | CTRL | 1619 | +0.053 | 5.3e-01 | 2321 | +0.057 | 4.1e-01 | **yes** |

### Summary

- positive-control fits: **32**
- observations: median 1700 -> **2356**
- median SE: 0.0831 -> **0.0620** (25% smaller)
- p-value improved in **17 of 32**
- beta shift beyond sampling noise (|z| > 1.96): **1 of 32**

DDC, the anchor signal, across every block that carries it:

- `olink_plasma`: +0.545 (P=1.4e-15, n=1624) -> **+0.575 (P=4.8e-27, n=2326)**
- `olink_csf`: +1.175 (P=2.2e-54, n=1199) -> **+1.004 (P=9.1e-77, n=2312)**
- `nulisa_cns_plasma`: +0.373 (P=7.9e-06, n=2120) -> **+0.351 (P=1.1e-08, n=2744)**
- `nulisa_cns_csf`: +1.307 (P=4.7e-74, n=1774) -> **+1.235 (P=1.8e-146, n=2382)**

