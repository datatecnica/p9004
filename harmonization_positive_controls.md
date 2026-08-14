# Does harmonization work? Longitudinal positive controls

Contrast `grp_NSD_vs_HC`. OLS `analyte ~ group + age_at_visit + SEX` with PATNO-clustered SEs, fitted over **all visits** carrying a value. Outcome z-scored so betas are per-SD. A random-intercept LMM is not identifiable here — these panels average ~1.2 observations per participant, so the variance component is unestimable and statsmodels returns a singular matrix.

**The reference is the project with more rows**, so the harmonized column is mostly uncorrected native values and correction error is minimised.

| block | **reference** | rows | mapped | rows |
|---|---|---|---|---|
| `olink_plasma` | **p314_Plasma** | 2,942 | p293_olink_plasma | 784 |
| `olink_csf` | **p314_CSF** | 2,098 | p277_CSF | 1,167 |
| `nulisa_cns_plasma` | **p288_CNS_plasma** | 2,665 | p312_Neuro_Plasma | 1,618 |
| `nulisa_cns_csf` | **p282_CNS_CSF** | 2,282 | p312_Neuro_CSF | 1,197 |
| `nulisa_inf_plasma` | **p288_Inflammation_plasma** | 2,669 | p312_Inflammation_Plasma | 1,616 |
| `nulisa_inf_csf` | **p282_Inflammation_CSF** | 1,899 | p312_Inflammation_CSF | 1,234 |

Harmonization should add both participants and timepoints, so a working correction gives **smaller p-values than the reference alone**. If p-values do not improve, the added samples are not carrying the effect.

## Positive controls

| block | analyte | slope | **ref** obs / ppl | ref beta | ref P | **harm** obs / ppl | harm beta | harm P | P better? | mapped beta | mapped P |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `nulisa_cns_csf` | CHI3L1 | 0.68 | 1774 / 1448 | -0.179 | 2.5e-02 | 2382 / 1730 | -0.104 | 7.6e-02 | no | +0.024 | 7.5e-01 |
| `nulisa_cns_csf` | CST3 | 0.63 | 1774 / 1448 | -0.048 | 5.7e-01 | 2382 / 1730 | -0.016 | 7.9e-01 | no | +0.000 | 1.0e+00 |
| `nulisa_cns_csf` | DDC | 1.10 | 1774 / 1448 | +1.307 | 4.7e-74 | 2382 / 1730 | +1.249 | 4.5e-146 | **yes** | +1.344 | 7.0e-109 |
| `nulisa_cns_csf` | GFAP | 0.93 | 1774 / 1448 | +0.048 | 5.6e-01 | 2382 / 1730 | +0.031 | 6.2e-01 | no | +0.087 | 2.9e-01 |
| `nulisa_cns_csf` | MAPT | 0.88 | 1774 / 1448 | -0.306 | 3.0e-04 | 2382 / 1730 | -0.191 | 4.7e-03 | no | -0.093 | 2.8e-01 |
| `nulisa_cns_csf` | NEFL | 0.93 | 1774 / 1448 | +0.130 | 6.8e-02 | 2382 / 1730 | +0.169 | 1.2e-03 | **yes** | +0.282 | 5.2e-05 |
| `nulisa_cns_csf` | SNCA | 0.82 | 1774 / 1448 | -0.304 | 5.1e-04 | 2382 / 1730 | -0.137 | 2.4e-02 | no | -0.090 | 2.7e-01 |
| `nulisa_cns_plasma` | CHI3L1 | 0.82 | 2120 / 1398 | -0.175 | 4.6e-02 | 2744 / 1785 | -0.225 | 1.2e-03 | **yes** | -0.228 | 3.6e-02 |
| `nulisa_cns_plasma` | CST3 | 0.93 | 2120 / 1398 | +0.052 | 5.0e-01 | 2744 / 1785 | +0.055 | 3.5e-01 | **yes** | +0.134 | 1.4e-01 |
| `nulisa_cns_plasma` | DDC | 1.01 | 2120 / 1398 | +0.373 | 7.9e-06 | 2744 / 1785 | +0.351 | 1.1e-08 | **yes** | +0.334 | 4.0e-04 |
| `nulisa_cns_plasma` | GFAP | 0.90 | 2120 / 1398 | +0.063 | 4.4e-01 | 2744 / 1785 | +0.035 | 5.5e-01 | no | +0.055 | 5.1e-01 |
| `nulisa_cns_plasma` | MAPT | 0.98 | 2120 / 1398 | +0.014 | 9.1e-01 | 2744 / 1785 | -0.068 | 3.9e-01 | **yes** | -0.107 | 2.8e-01 |
| `nulisa_cns_plasma` | NEFL | 0.92 | 2120 / 1398 | +0.220 | 1.2e-03 | 2744 / 1785 | +0.206 | 2.1e-05 | **yes** | +0.261 | 2.3e-04 |
| `nulisa_cns_plasma` | SNCA | 1.17 | 2120 / 1398 | +0.031 | 7.3e-01 | 2744 / 1785 | +0.035 | 6.1e-01 | **yes** | -0.020 | 8.4e-01 |
| `nulisa_inf_csf` | CHI3L1 | 0.86 | 1462 / 1249 | -0.201 | 2.4e-02 | 2090 / 1582 | -0.085 | 1.8e-01 | no | +0.058 | 4.6e-01 |
| `nulisa_inf_csf` | GFAP | 0.88 | 1462 / 1249 | +0.052 | 5.7e-01 | 2090 / 1582 | +0.038 | 5.6e-01 | **yes** | +0.097 | 2.3e-01 |
| `nulisa_inf_plasma` | CHI3L1 | 0.91 | 2119 / 1398 | -0.179 | 5.0e-02 | 2746 / 1788 | -0.216 | 1.8e-03 | **yes** | -0.222 | 3.7e-02 |
| `nulisa_inf_plasma` | GFAP | 0.90 | 2119 / 1398 | +0.029 | 7.1e-01 | 2746 / 1788 | +0.019 | 7.4e-01 | no | +0.085 | 3.2e-01 |
| `olink_csf` | CHI3L1 | 1.38 | 1156 / 976 | +0.126 | 1.3e-01 | 2271 / 1569 | -0.074 | 2.6e-01 | no | -0.158 | 3.4e-02 |
| `olink_csf` | CST3 | 0.83 | 1100 / 937 | +0.206 | 3.4e-02 | 2218 / 1535 | +0.172 | 1.7e-02 | **yes** | +0.101 | 1.7e-01 |
| `olink_csf` | DDC | 1.18 | 1199 / 1014 | +1.175 | 2.2e-54 | 2312 / 1596 | +1.009 | 1.3e-74 | **yes** | +0.860 | 6.9e-34 |
| `olink_csf` | GFAP | 1.56 | 1200 / 1015 | +0.028 | 7.1e-01 | 2314 / 1598 | -0.142 | 1.1e-02 | **yes** | -0.237 | 1.4e-03 |
| `olink_csf` | MAPT | 1.33 | 1200 / 1015 | +0.023 | 8.1e-01 | 2314 / 1598 | -0.266 | 2.7e-04 | **yes** | -0.385 | 2.1e-05 |
| `olink_csf` | NEFL | 1.48 | 1196 / 1011 | +0.201 | 8.9e-03 | 2310 / 1594 | +0.016 | 7.7e-01 | no | -0.094 | 1.9e-01 |
| `olink_csf` | SNCA | 1.11 | 1199 / 1014 | +0.075 | 3.7e-01 | 2312 / 1596 | -0.083 | 1.1e-01 | **yes** | -0.130 | 4.9e-02 |
| `olink_plasma` | CHI3L1 | 1.09 | 1626 / 1301 | -0.138 | 1.3e-01 | 2330 / 1531 | -0.181 | 8.1e-03 | **yes** | -0.176 | 7.9e-02 |
| `olink_plasma` | CST3 | 0.95 | 1626 / 1301 | +0.102 | 1.8e-01 | 2330 / 1531 | -0.010 | 8.8e-01 | no | -0.115 | 2.6e-01 |
| `olink_plasma` | DDC | 0.85 | 1624 / 1302 | +0.545 | 1.4e-15 | 2326 / 1531 | +0.554 | 3.3e-26 | **yes** | +0.790 | 2.3e-23 |
| `olink_plasma` | GFAP | 1.00 | 1624 / 1302 | +0.083 | 2.6e-01 | 2326 / 1531 | +0.034 | 5.9e-01 | no | +0.045 | 6.4e-01 |
| `olink_plasma` | MAPT | 0.48 | 1579 / 1279 | -0.201 | 4.4e-03 | 2283 / 1510 | -0.120 | 2.6e-02 | no | +0.045 | 6.4e-01 |
| `olink_plasma` | NEFL | 0.76 | 1623 / 1301 | +0.267 | 3.4e-05 | 2325 / 1530 | +0.189 | 4.3e-04 | no | +0.189 | 2.6e-02 |
| `olink_plasma` | SNCA | 0.82 | 1624 / 1302 | -0.039 | 6.2e-01 | 2326 / 1531 | -0.082 | 2.3e-01 | **yes** | -0.109 | 3.4e-01 |

## Negative controls

| block | analyte | slope | **ref** obs / ppl | ref beta | ref P | **harm** obs / ppl | harm beta | harm P | P better? | mapped beta | mapped P |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `olink_csf` | CTRL | 1.48 | 1196 / 1012 | -0.010 | 9.0e-01 | 2309 / 1595 | -0.205 | 2.2e-04 | **yes** | -0.262 | 3.9e-04 |
| `olink_plasma` | CTRL | 0.93 | 1619 / 1299 | +0.053 | 5.3e-01 | 2321 / 1528 | +0.056 | 4.0e-01 | **yes** | +0.055 | 6.2e-01 |

## Summary

- positive-control fits: **32**
- observations: reference median 1700 -> harmonized 2356 (**+656**)
- participants: reference median 1350 -> harmonized 1664
- median SE: 0.0831 -> **0.0610** (27% smaller)
- **p-value improved in 18 of 32 fits overall**
- among the 15 where the reference was already significant, improved in **8**
- beta shift beyond sampling noise (|z| > 1.96): **1 of 32**
- constituent projects differ (P_het < 0.05): **8 of 32** (1.6 expected by chance)

