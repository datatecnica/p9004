# Calibration report

## Negative controls

_No control analytes found in the results directory._


## Blocked permutation null

Group labels shuffled between subjects within blocks of (n_obs, mean-time tertile), preserving the real follow-up imbalance. A calibrated run returns λ ≈ 1 and ~5% at P<0.05.

Each run is permuted **10×** over the same fixed sample of 100 analytes, so the spread below is permutation variability alone, not analyte sampling. A single replicate reports a point with no error bar; the **range** column is what says whether that point could be trusted.

| Run | Term | perms | n fits | λ mean ± sd | λ range | % P<0.05 mean ± sd | % range | warn |
|---|---|---|---|---|---|---|---|---|
| `trajectory_NSD_vs_HC` | within | 10 | 999 | 1.028 ± 0.415 | 0.623–2.038 | 5.5% ± 3.5 | 1.0–12.0% | ⚠️ unstable |
| `trajectory_NSD_vs_HC` | between | 10 | 1,000 | 1.128 ± 0.228 | 0.684–1.416 | 7.5% ± 1.3 | 5.0–9.0% | ok |
| `trajectory_CI_MOCA` | within | 10 | 999 | 0.764 ± 0.162 | 0.599–0.985 | 3.6% ± 1.9 | 2.0–7.0% | ok |
| `trajectory_CI_MOCA` | between | 10 | 1,000 | 1.099 ± 0.289 | 0.522–1.568 | 5.0% ± 2.1 | 2.0–9.0% | ok |
| `trajectory_stage_2A_vs_2B` | within | 10 | 995 | 0.802 ± 0.214 | 0.488–1.190 | 4.7% ± 2.8 | 1.0–9.2% | ok |
| `trajectory_stage_2A_vs_2B` | between | 10 | 1,000 | 1.156 ± 0.388 | 0.618–1.783 | 6.6% ± 3.2 | 2.0–11.0% | ⚠️ unstable |
| `trajectory_CI_PI` | within | 10 | 998 | 0.863 ± 0.224 | 0.511–1.262 | 5.0% ± 3.3 | 2.0–13.0% | ⚠️ unstable |
| `trajectory_CI_PI` | between | 10 | 1,000 | 1.124 ± 0.430 | 0.540–1.778 | 7.9% ± 4.3 | 3.0–17.0% | ⚠️ unstable |
| `trajectory_HC_vs_Prodromal` | within | 10 | 1,000 | 0.991 ± 0.219 | 0.699–1.292 | 3.1% ± 1.3 | 1.0–5.0% | ok |
| `trajectory_HC_vs_Prodromal` | between | 10 | 1,000 | 1.049 ± 0.259 | 0.780–1.555 | 4.9% ± 2.8 | 2.0–11.0% | ⚠️ unstable |
| `trajectory_HC_vs_PD` | within | 10 | 1,000 | 1.053 ± 0.393 | 0.646–2.030 | 5.5% ± 6.3 | 0.0–21.0% | ⚠️ unstable |
| `trajectory_HC_vs_PD` | between | 10 | 1,000 | 0.917 ± 0.195 | 0.661–1.219 | 4.4% ± 2.2 | 1.0–8.0% | ok |
| `trajectory_Prodromal_vs_PD` | within | 10 | 999 | 1.082 ± 0.409 | 0.630–1.965 | 4.1% ± 2.2 | 1.0–8.0% | ok |
| `trajectory_Prodromal_vs_PD` | between | 10 | 1,000 | 1.227 ± 0.443 | 0.954–2.449 | 7.3% ± 4.4 | 3.0–18.0% | ⚠️ unstable |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | within | 10 | 999 | 0.798 ± 0.235 | 0.407–1.120 | 3.4% ± 2.7 | 0.0–10.0% | ok |
| `trajectory_HCNSDneg_vs_ProdNSDpos` | between | 10 | 1,000 | 1.090 ± 0.281 | 0.777–1.583 | 5.2% ± 3.4 | 1.0–11.0% | ⚠️ unstable |
| `trajectory_HCNSDneg_vs_PDNSDpos` | within | 10 | 998 | 0.810 ± 0.317 | 0.505–1.502 | 4.0% ± 2.4 | 1.0–8.0% | ok |
| `trajectory_HCNSDneg_vs_PDNSDpos` | between | 10 | 1,000 | 0.967 ± 0.302 | 0.649–1.724 | 5.0% ± 2.5 | 2.0–10.0% | ok |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | within | 10 | 999 | 0.886 ± 0.363 | 0.399–1.424 | 5.0% ± 4.0 | 1.0–15.0% | ⚠️ unstable |
| `trajectory_ProdNSDpos_vs_PDNSDpos` | between | 10 | 1,000 | 1.192 ± 0.269 | 0.824–1.612 | 7.8% ± 2.5 | 3.0–12.0% | ⚠️ unstable |

**13** permutation replicate(s) exceeded the 10% type-I limit across all runs tested. A run flagged `unstable` had replicates on both sides of it.

