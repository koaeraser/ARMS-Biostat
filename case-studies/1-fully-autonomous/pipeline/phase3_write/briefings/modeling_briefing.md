# Modeling Briefing

## Evaluations Run

| Evaluation Type | Scenarios | B | Wall Time | Output CSV |
|----------------|-----------|---|-----------|------------|
| LOO-CV (all 5 methods) | 35 arms x 5 methods | N/A (deterministic) | ~2s | data/loocv_results.csv, data/loocv_summary.csv |
| Sensitivity to beta | 9 values [1, 3, 5, 8, 10, 15, 20, 30, 50] | N/A | ~10s | data/sensitivity_beta.csv |
| Sensitivity to w0 | 6 values [0, 0.05, 0.1, 0.2, 0.3, 0.5] | N/A | ~6s | data/sensitivity_w0.csv |
| Contamination (external) | 4 scenarios (Clean, +Adj, +For, +Both) | N/A | ~15s | data/contamination_results.csv |
| Contamination (adversarial) | 9 scenarios (varying n_fake, sim, rate) | N/A | ~20s | data/contamination_stress_results.csv |
| Ablation | 7 variants | N/A | ~8s | data/ablation_results.csv |
| Simulation study | 3 scenarios (Favorable, Adverse, Mixed) | 500 | ~TBD | data/simulation_summary.csv, data/simulation_arm_level.csv |

**Note**: KG-DAP is fully deterministic (no MCMC sampling). The LOO-CV, sensitivity, contamination, and ablation evaluations produce identical results across runs. Monte Carlo replications (B) are only needed for the simulation study where data is randomly generated.

## Formula-Code Consistency Check

| Equation | Spec Location | Code Location | Status |
|----------|--------------|---------------|--------|
| Composite similarity S_ij | Step 1 | kg_dap.py:compute_similarity_matrix() -> build_kg.py | MATCH |
| Power weights omega_k = S^beta/sum(S^beta) | Step 2 | kg_dap.py:compute_power_weights() lines 175-220 | MATCH |
| Historical arm posteriors pi_k = Beta(a0+y_k, b0+n_k-y_k) | Step 4 | kg_dap.py:construct_prior() lines 267-270 | MATCH |
| Sample-size cap scaling | Step 4 | kg_dap.py:construct_prior() line 268 | MATCH |
| Vague component Beta(1,1) | Step 5 | kg_dap.py:construct_prior() lines 263-264 | MATCH |
| KG-DAP prior mixture weights | Step 6 | kg_dap.py:construct_prior() lines 254-256 | MATCH |
| Posterior update via log-marginal | Step 7 | kg_dap.py:posterior_update() lines 310-321 | MATCH |
| ESS via moment-matching | Step 8 | kg_dap.py:compute_ess() lines 388-398 | MATCH (uses kappa = mu*(1-mu)/var - 1) |
| Beta-Binomial predictive PMF | Step 9 | kg_dap.py:log_beta_binomial_pmf() lines 404-411 | MATCH |
| Mixture predictive (log-sum-exp) | Step 9 | kg_dap.py:log_predictive_pmf() lines 414-427 | MATCH |
| Mixture mean | Misc | kg_dap.py:mixture_mean() lines 329-335 | MATCH |
| Mixture variance (law of total variance) | Misc | kg_dap.py:mixture_variance() lines 338-349 | MATCH |
| rMAP: DerSimonian-Laird tau^2 | Comparator 1 | comparators.py:rmap_prior() lines 141-170 | MATCH |
| rMAP: Logit-normal moment matching | Comparator 1 | comparators.py:rmap_prior() lines 172-200 | MATCH |
| Pooled prior | Comparator 3 | comparators.py:pooled_prior() lines 45-70 | MATCH |
| Equal-weight mixture | Comparator 4 | comparators.py:equal_weight_prior() lines 76-102 | MATCH |

### Mismatches Found
**NONE**. All 16 equations in the methodology specification match the validated code exactly.

### ESS Formula Detail
The spec says: ESS = mu*(1-mu)/var - 1. The code computes kappa = mu*(1-mu)/var - 1 and returns max(kappa, 0). This is consistent — ESS = kappa where kappa is the precision parameter of the moment-matched Beta.

## Comparator Coverage

| Evaluation | Comparators Included | All Present? |
|-----------|---------------------|-------------|
| LOO-CV | KG-DAP, rMAP, Uniform, Pooled, EqualWeight | YES |
| Sensitivity (beta) | KG-DAP only (parameter varies) | N/A (single-method sweep) |
| Sensitivity (w0) | KG-DAP only (parameter varies) | N/A (single-method sweep) |
| Contamination (external) | KG-DAP only | N/A (robustness test) |
| Contamination (adversarial) | KG-DAP only | N/A (robustness test) |
| Ablation | 7 KG-DAP variants + EqualWeight | YES (EqualWeight = no-graph baseline) |
| Simulation study | KG-DAP, rMAP, Uniform, Pooled, EqualWeight | YES |

## Consistency Check vs Phase 2

| Metric | Phase 2 (B=N/A, deterministic) | Production | Within 2x MC SE? |
|--------|-------------------------------|------------|-------------------|
| KG-DAP MAE | 0.1257 | 0.1257 | YES (identical, deterministic) |
| rMAP MAE | 0.1697 | 0.1697 | YES (identical) |
| KG-DAP Coverage | 0.943 | 0.943 | YES (identical) |
| KG-DAP RMSE | 0.1559 | 0.1559 | YES (identical) |

**Note**: Since KG-DAP and all comparators (for LOO-CV) are fully deterministic, Phase 2 and production results are identical. No MC SE applies.

## Monte Carlo Standard Errors

MC-SEs are relevant ONLY for the simulation study (B=500).
Will be reported after simulation completes.

All MC-SEs reported: PENDING (simulation running)
Any SE > 50% of effect size: PENDING

## Figures Produced

- figures/fig1_loocv_mae.pdf — LOO-CV MAE comparison bar chart (5 methods)
- figures/fig2_arm_level.pdf — Arm-level scatter: prior mean vs observed rate
- figures/fig3_sensitivity_beta.pdf — MAE and coverage vs beta (dual axis, with safe range)
- figures/fig4_ablation.pdf — Component contribution analysis horizontal bar chart
- figures/fig5_contamination.pdf — Contamination robustness (sim=0.8 and sim=0.9 panels)
- figures/fig6_sensitivity_w0.pdf — MAE and coverage vs w0 (supplementary)

## Anomalies or Warnings

1. **Pooled prior has 2.9% coverage** — Expected: the pooled prior is far too concentrated (ESS=9075) and misses the true rate for nearly all arms. This is not a bug; it demonstrates the danger of naive pooling.

2. **KG-DAP does not beat EqualWeight on log-predictive score** — EqualWeight: -5.070 vs KG-DAP: -5.028. KG-DAP IS better (less negative = higher score). The difference is modest because log-predictive score depends on both mean accuracy and interval calibration; KG-DAP's narrower intervals give it an edge.

3. **rMAP has sub-nominal coverage (88.6%)** — This reflects rMAP's tendency to produce overly concentrated priors when the DerSimonian-Laird estimator underestimates heterogeneity. The logit-normal moment-matching approximation may also contribute.

4. **Sample-size cap has negligible effect** — In the ablation, "No Cap" produces identical MAE to "Full KG-DAP" (0.1257). This is because most NDMM trials have n < 200, so the cap rarely binds. Document as a finding, not a problem.
