## Paper Assessment

### Manuscript: Knowledge-Graph-Driven Adaptive Prior for Bayesian Borrowing in Clinical Trials
### Date: 2026-03-23
### Round: 2 (FINAL) — Previous scores: Round 0 = 40/50, Round 1 = 41/50

---

### Verification of Round 2 Fixes

Before grading, the four claimed fixes from Round 2 are verified against ground truth:

1. **Drug-Only delta corrected to +7.3%.** VERIFIED. Table 3 line 444 shows `$+$7.3` and the prose at line 459 says `Drug-only ($+$7.3\%)`. The validated data gives (0.1348 - 0.1257)/0.1257 = 7.2523%, which rounds to 7.3% at one decimal place. This is correct; Round 1 had incorrectly changed it to 7.2%.

2. **ji2022dynamic and kaizer2021bayesian citations verified as present.** VERIFIED. Line 98 cites both: `\citep[iMEM;][]{kaizer2021bayesian}` and `\citep[dMEM;][]{ji2022dynamic}`. Both entries are in references.bib (lines 31-49). Cross-reference check confirms zero orphaned bib entries and zero missing citations -- perfect alignment between manuscript and bib file.

3. **Bootstrap CIs updated to 4 decimal places.** VERIFIED. Table 4 (lines 540, 547, 554) now shows `$[-0.0440, -0.0438]$`, `$[0.0012, 0.0014]$`, and `$[-0.0302, -0.0300]$`. The prose (lines 525, 529) matches. The validated CSV data confirms: Favorable CI = [-0.0440, -0.0438], Adverse CI = [0.0012, 0.0014], Mixed CI = [-0.0302, -0.0300]. All values match exactly.

4. **Proposition 2 proof sketch expanded.** VERIFIED. Lines 301-303 now include the law of total variance decomposition, the heuristic step acknowledgment ("The argument is heuristic at the step relating mixture-mean..."), and numerical verification ("observed maximum ESS = 10.2; bound = 204"). This is a substantial improvement over the previous terse sketch.

**All four fixes are genuinely implemented and correct. No regressions detected.**

---

### 1. Correctness: 4/5

All mathematical derivations remain correct. The five propositions and their proofs are valid. The LOO-CV table (Table 2) matches the validated `loocv_summary.csv` to 4 significant figures: KG-DAP MAE=0.1257 (CSV: 0.12567), rMAP MAE=0.1697 (CSV: 0.16967), coverage 94.3% (CSV: 0.9429). The simulation table (Table 4) matches `simulation_summary_B1000.csv`: Favorable KG-DAP MAE=0.1258 (CSV: 0.12580), Mixed KG-DAP MAE=0.1560 (CSV: 0.15598). The design OC table (Table 6) matches `design_oc_representative.csv` exactly to 3 decimal places for all 30 cells. The Drug-Only delta is now correctly +7.3%. The bootstrap CIs at 4dp are now informative and match the CSV data. All bibliography entries are cited and all citations have corresponding bib entries. The improvement percentages (25.9%, 22.4%, 21.8%, 16.2%) are all verified correct. The MC SE claim ("below 0.0002") is accurate. The one remaining issue preventing a score of 5 is that Proposition 2's proof sketch explicitly acknowledges a heuristic step without providing a full proof -- but the label is honest, numerical verification is provided, and the bound holds with substantial slack (observed max ESS=10.2 vs bound=204), so the incompleteness is transparently communicated rather than hidden.

### 2. Completeness: 5/5

The manuscript covers all standard and extended components expected for a Biometrics methods paper. Primary evaluation: LOO-CV across 35 arms with all 5 comparators (Table 2). Sensitivity analyses: power parameter beta (Table S3, Figure 3) and robustness weight w0 (Table 3, Figure S1). Robustness: contamination test at 8 scenarios varying number of fake arms, similarity level, and fake rate (Table 5, Figure 5). Component contribution analysis: 6 ablation variants isolating each component (Table 3, Figure 4). Simulation study: B=1,000 replications, 3 scenarios (Favorable/Adverse/Mixed), 5 methods, with paired bootstrap CIs (Table 4). Design operating characteristics: go/no-go power curves at 10 true rates, 3 methods, 6 representative arms, B=5,000 (Table 6). Both prior and posterior ESS reported (Table 2). Supplementary material includes trial arms table (S1), borrowing weights for representative arm (S2), and full beta sensitivity (S3). Six figures referenced and present. This coverage is exhaustive.

### 3. Rigor: 4/5

Uncertainty quantification is thorough. The B=1,000 simulation with deterministic seeds (SCENARIO_SEEDS dict, formula `42 + b*997 + scenario_seed*9973`) ensures full cross-session reproducibility. Paired bootstrap 95% CIs establish statistical significance for the KG-DAP vs rMAP MAE difference in Favorable and Mixed scenarios, and correctly show non-significance in Adverse (CI includes positive values). The CIs at 4dp are now informative rather than appearing degenerate. MC standard errors are accurately reported below 0.0002. Sensitivity analyses cover both key hyperparameters with sufficient resolution. The type I error discussion in Section 5.6 is honest and notes that threshold recalibration is standard practice. Proposition 2's expanded proof sketch with law of total variance and numerical verification is much stronger than before, though still explicitly heuristic. Prior sensitivity to the composite similarity weights (alpha, beta_w, gamma) is acknowledged but not systematically explored. No formal power analysis for frequentist sample size is provided, though the design OC section partially addresses this gap. Score remains 4/5 -- the gap to 5 would require a fully rigorous Proposition 2 proof and systematic alpha/beta_w/gamma sensitivity.

### 4. Clarity: 5/5

The manuscript is well-organized with a logical flow from problem statement through related work, methods, theory, evaluation, and discussion. Notation is consistent throughout. The abstract is quantitative and informative. The contribution statement (four numbered items) is precise. Related work is comprehensive with fair characterization of all major competing approaches. Remark 1 (power weights vs diffusion kernel) is exemplary transparent scientific reporting. The terminology is now standardized: "component contribution analysis" is used consistently in all reader-visible text. All bibliography entries align perfectly with citations -- no orphans or missing entries. The limitations section is honest and specific (five numbered items). The discussion thoughtfully addresses where KG-DAP does not excel (few historical trials, patient-level data available). The bootstrap CIs at 4dp are now clearly intervals rather than appearing as point estimates. The Proposition 2 proof sketch's honest acknowledgment of its heuristic step ("a fully rigorous proof would require...") is a model of intellectual honesty. Upgrading from 4 to 5 because the Round 2 fixes resolved the three remaining clarity issues: dead bib entries (now zero), degenerate-looking CIs (now 4dp), and Drug-Only delta rounding (now correct).

### 5. Novelty: 4/5

Unchanged from prior rounds. KG-DAP represents a genuine cross-disciplinary contribution: importing knowledge graph structure into Bayesian borrowing weight construction. The combination of power-sharpened KG-derived weights, finite Beta mixture architecture, full pre-specifiability, and closed-form computation is non-obvious. Table 1 clearly demonstrates that no existing method simultaneously provides trial-specific weights, full pre-specifiability, external structure, and closed-form computation. The key insight ("what you borrow from matters more than how much") is well-supported by the modest ESS (~6) producing disproportionately large MAE improvements (26%). The cross-domain bridge from knowledge graphs to Bayesian priors requires insights from multiple fields. However, the power-sharpened weights are a relatively simple softmax function, not the more theoretically novel diffusion kernel originally proposed. Score: 4/5.

### 6. Impact: 4/5

The design operating characteristics analysis (Section 5.6) quantifies the method's practical value: KG-DAP increases the probability of a correct go decision by 7-12 percentage points at clinically relevant alternatives (theta = 0.35-0.40) compared to no borrowing. Full pre-specifiability aligns with FDA guidance on pre-specified informative priors, providing a clear regulatory path. Closed-form computation enables the extensive sensitivity analyses that responsible application demands. The 25.9% MAE improvement is practically meaningful. The method remains limited by single-domain evaluation (NDMM only), binary endpoints only, and the requirement for a curated knowledge graph. The design OC analysis uses an average over 6 representative arms rather than all 35, which somewhat limits its generalizability. Score: 4/5 -- a second disease area or endpoint type would be needed for 5.

### 7. Performance: 4/5

The empirical evidence is comprehensive and well-quantified. In LOO-CV, KG-DAP achieves 25.9% MAE improvement over rMAP with well-calibrated 94.3% coverage. The B=1,000 simulation confirms these findings under known ground truth, with paired bootstrap CIs establishing statistical significance in Favorable (CI: [-0.0440, -0.0438]) and Mixed (CI: [-0.0302, -0.0300]) scenarios. Under Adverse conditions, KG-DAP gracefully degrades (0.8% worse than rMAP), confirming the safety property. The component contribution analysis cleanly isolates contributions of graph structure (+32.6%) and power weights (+30.8%). Design OC translates these improvements into decision-relevant quantities. All four non-trivial comparators are consistently present across all analyses. Score: 4/5 because (a) the evaluation is entirely prior-predictive (LOO-CV) or simulation-based -- no posterior evaluation with actual current-trial data is performed, and (b) no formal frequentist power/sample-size analysis is provided despite the design OC section.

---

### Raw Score: 29/35
### Weighted Score: 42/50
  Execution (x1.0): Correctness 4 + Completeness 5 + Rigor 4 + Clarity 5 = 18/20
  Research  (x2.0): (Novelty 4 + Impact 4 + Performance 4) x 2.0 = 24/30
  Weighted total = 18 + 24 = 42/50
### Grade: A (>=42)

### Score Progression
| Dimension | Round 0 | Round 1 | Round 2 | Change R1->R2 |
|-----------|---------|---------|---------|---------------|
| Correctness | 4 | 4 | 4 | -- |
| Completeness | 4 | 5 | 5 | -- |
| Rigor | 4 | 4 | 4 | -- |
| Clarity | 4 | 4 | **5** | +1 |
| Novelty | 4 | 4 | 4 | -- |
| Impact | 4 | 4 | 4 | -- |
| Performance | 4 | 4 | 4 | -- |
| **Weighted** | **40** | **41** | **42** | **+1** |

### Code Audit Summary

| Category | PASS | FLAG | FAIL |
|----------|------|------|------|
| Formula Correctness (1-6) | 5 | 0 | 0 |
| Parameter Consistency (7-11) | 5 | 0 | 0 |
| Methodological Completeness (12-15) | 4 | 0 | 0 |

**Formula Correctness**: RMSE (PASS -- sqrt of mean of squared errors in validated code), Coverage (PASS -- non-strict inequality via mixture_interval), ESS (PASS -- mu*(1-mu)/var - 1 matches Morita 2008 and Equation 6), log-predictive (PASS -- sum of log Beta-Binomial PMFs), posterior weights (PASS -- log-sum-exp normalization in posterior_update).

**Parameter Consistency**: All hyperparameters match between code (DEFAULT_PARAMS in run_validation.py) and manuscript: beta=15.0, w0=0.20, alpha/beta_w/gamma=0.20/0.20/0.60, n_cap=200, alpha0=beta0=1.0 (PASS). B=1,000 in simulation matches table caption (PASS). Deterministic seeds via SCENARIO_SEEDS dict (PASS). Decision threshold theta_0=0.30, gamma_go=0.80 consistent between Section 5.6 text and table caption (PASS). All citation keys in manuscript match bib file exactly (PASS).

**Methodological Completeness**: All 5 comparators present in all evaluations (PASS). Both prior and posterior ESS reported (PASS). Power analysis via design OC section (PASS). Paired bootstrap CIs for method comparisons (PASS).

### Critical Findings
None (no FAIL items, no FLAG items).

### Top 3 Strengths

1. **Internal consistency is now impeccable.** All numbers in the manuscript match validated CSV data to 4 significant figures across all 7 tables. All citation keys align perfectly with the bib file. The Drug-Only delta, bootstrap CIs, and MC SE claims are all verified correct. This level of internal consistency is rare even in published papers.

2. **Design operating characteristics close the gap to practice.** Section 5.6 directly translates MAE improvements into go/no-go decision probabilities that regulators and trial designers can evaluate. The honest discussion of type I error inflation with the standard recalibration solution is exemplary.

3. **Exhaustive evaluation with honest scientific reporting.** LOO-CV, sensitivity analyses (2 hyperparameters), contamination robustness (8 scenarios), component contribution analysis (6 variants), simulation (3 scenarios, B=1000, bootstrap CIs), and design OC (10 true rates, B=5000) tell a coherent, quantified story. Remark 1 on the failed diffusion kernel and the five specific limitations are models of transparent reporting.

### Top 3 Weaknesses (with suggested improvement)

1. **Proposition 2 remains a proof sketch, not a proof.** The law of total variance argument and numerical verification are substantial improvements, but the heuristic step (bounding mu(1-mu) relative to component terms) is acknowledged but unresolved. -- Suggested improvement: Either provide a fully rigorous proof using Jensen's inequality on the concave function f(x) = x(1-x) (which would give mu(1-mu) >= sum_k w_k mu_k(1-mu_k)), or downgrade the proposition to a "conjecture with numerical evidence" and add the bound as an empirical observation rather than a theoretical result.

2. **Evaluation is entirely prior-predictive or simulation-based.** The LOO-CV evaluates prior predictive accuracy, not posterior accuracy after observing current-trial data. No experiment simulates the arrival of actual current-trial data and evaluates how the posterior performs. -- Suggested improvement: Add a "posterior calibration" experiment where for each held-out arm, you simulate observing partial data (e.g., first 50% of subjects) and evaluate the posterior prediction for the remaining subjects.

3. **Single therapeutic area limits generalizability claims.** All evidence comes from NDMM MRD negativity trials. The method's performance in therapeutic areas with different graph densities, different similarity distributions, or different baseline heterogeneity is unknown. -- Suggested improvement: The planned future work (Section 6.4) mentions this; a second disease area in revision would substantially strengthen the paper.

### Fixable Issues (for paper-fixer)

1. **Proposition 2 proof note (optional).** Consider adding a brief footnote noting that for the Beta family, the concavity of f(x) = x(1-x) on [0,1] gives mu(1-mu) >= sum_k w_k mu_k(1-mu_k) by Jensen's inequality, which would close the heuristic gap. This is a 1-line fix.

2. **Abstract length (minor).** The abstract is 16 lines of dense text. Consider trimming by 2-3 lines to improve readability. This does not affect the score but would improve the reading experience.

3. **No remaining correctness or consistency issues.** All previous fixable issues have been resolved.

### Verdict

Round 2 achieves a Grade A score of 42/50, meeting the target threshold. The +1 improvement over Round 1 comes from the Clarity dimension, where all three previously identified issues (dead bib entries, degenerate-looking bootstrap CIs, Drug-Only delta rounding error) have been correctly resolved. The manuscript now has zero internal inconsistencies between tables, text, and validated data -- a level of quality that is unusual even among published papers. All 17 citation keys align perfectly between manuscript and bib file. The remaining gaps to higher scores are structural (Proposition 2 proof completeness, single disease area, prior-predictive-only evaluation) and cannot be addressed by text editing alone. Current assessment: **acceptable for submission** to Biometrics, with the expectation that reviewers may request (a) a second disease area, (b) a rigorous Proposition 2 proof, and/or (c) a posterior calibration experiment. The paper's strengths -- genuine cross-disciplinary novelty, regulatory-aligned design, exhaustive evaluation, and honest reporting -- position it well for a constructive review process.
