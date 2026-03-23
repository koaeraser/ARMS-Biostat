## Paper Assessment

### Manuscript: Knowledge-Graph-Driven Adaptive Prior for Bayesian Borrowing in Clinical Trials
### Date: 2026-03-23

---

### 1. Correctness: 4/5

The mathematical derivations are correct throughout. All five propositions and their proofs are valid: conjugate closure follows from standard Beta-Binomial conjugacy, bounded ESS from mixture variance bounds, asymptotic dominance from Bernstein--von Mises, and monotonicity/limit properties from elementary analysis of the power function. The numbers in the manuscript tables match the validated CSV data to at least 3-4 significant figures for all primary metrics (MAE, RMSE, coverage, ESS). Relative improvement percentages (25.9%, 22.4%, 21.8%) are verified correct. One correctness issue: the manuscript states "All Monte Carlo standard errors are below 0.02% of the corresponding effect sizes" (Section 5.5, line 521). This is incorrect by approximately an order of magnitude -- the maximum SE is 0.000162, which is approximately 0.37% of the KG-DAP vs rMAP effect size (0.0438), not below 0.02%. The SEs are genuinely small and the simulation conclusions are reliable, but the specific claim is wrong. Code audit confirms: RMSE formula is correct (sqrt of mean of squared errors), coverage uses non-strict inequality on both sides (standard), ESS formula matches specification.

### 2. Completeness: 4/5

The manuscript covers all standard sections expected for a Biometrics methods paper: introduction, related work, methods, theory, empirical evaluation, and discussion. The empirical evaluation is comprehensive, including LOO-CV (primary), sensitivity analyses for both key hyperparameters (beta and w0), a contamination robustness test, a component contribution analysis (ablation), and a simulation study with three scenarios (favorable/adverse/mixed). All four comparators (rMAP, Uniform, Pooled, EqualWeight) appear in all result tables. The manuscript has six figures, all referenced in text. However, several completeness gaps exist: (a) no design operating characteristics (power, type I error) are reported -- these are standard for Biometrics borrowing papers and practically essential for regulatory applications; (b) no posterior ESS is reported, only prior ESS; (c) the contamination test in the manuscript (Table 5) uses fabricated arms with fixed similarity, not the original external adjacent/foreign trials from the validated results (contamination_results.csv), so the external trial contamination results are lost; (d) three bibliography entries (lin2018borrowing, tierney1986accurate, clemen1989combining) are included in references.bib but never cited in the manuscript text; (e) supplementary material is minimal (one table, one figure) given the depth of analyses available.

### 3. Rigor: 4/5

Uncertainty quantification is present throughout: 95% coverage is reported for all LOO-CV and sensitivity results, credible interval widths are tabulated, and Monte Carlo standard errors are reported for the simulation study (even if the prose claim about their magnitude is incorrect). The sensitivity analyses cover both key hyperparameters with sufficient resolution to establish the operating plateau. The contamination analysis systematically varies both the number and similarity of misleading arms. The component contribution analysis is well-designed, isolating individual contributions. However, the "proof sketch" for Proposition 2 (bounded ESS) is not a complete proof -- the claim that "mixture variance is bounded below by the variance of a single component" requires more careful argument (the variance of a mixture is not simply bounded this way; it depends on the spread of component means). The simulation study uses hash(sname) for seed construction (run_simulation_v3.py line 81), which is non-reproducible across Python sessions since Python randomizes hash seeds by default. This is a reproducibility concern, though the overall conclusions are sound given the deterministic LOO-CV.

### 4. Clarity: 4/5

The manuscript is well-organized with a logical flow from problem statement through method, theory, evaluation, and discussion. Notation is consistent throughout (the use of beta for both the power parameter and Beta distributions is mildly confusing but standard). The abstract is informative and appropriately quantitative. The contribution statement (four numbered items in the introduction) is precise. The related work section is comprehensive, covering all major competing approaches with fair characterization. Remark 1 (choice of power weights over diffusion kernel) is a model of transparent scientific reporting -- explaining why the original approach failed and what replaced it. The discussion honestly acknowledges limitations, including the pre-specification trade-off and single-domain evaluation. Minor issues: (a) the keyword "Beta mixture" in the abstract keywords could be more specific ("Beta mixture prior"); (b) the manuscript uses "component contribution analysis" and "ablation" interchangeably (Section 5.4 header says "Component Contribution Analysis" but Figure 4 caption says "ablation variant"); (c) the supplementary material feels hastily assembled rather than integral to the paper.

### 5. Novelty: 4/5

KG-DAP represents a genuine cross-disciplinary contribution: importing knowledge graph structure into Bayesian borrowing weight construction. While individual components are standard (Beta mixtures, power-sharpened softmax weights, robustness via vague component mixing), the combination is non-obvious and addresses a real gap in the literature. The manuscript's Table 1 clearly demonstrates that no existing method simultaneously provides trial-specific weights, full pre-specifiability, external structure, and closed-form computation. The key insight that "what you borrow from matters more than how much" is well-supported by the modest ESS (~6) producing disproportionately large MAE improvements (26%). The novelty is somewhat tempered by the fact that the power-sharpened weights are a simple softmax function, not the more theoretically elegant diffusion kernel originally proposed -- but the manuscript is commendably transparent about this pragmatic choice. The novelty is clearly above "meaningful extension" (score 3) because the cross-domain bridge from knowledge graphs to Bayesian priors is non-trivial and the resulting framework opens new research directions.

### 6. Impact: 4/5

KG-DAP addresses a genuine practical problem: how to incorporate structural knowledge about trial relationships into borrowing weights while maintaining full pre-specifiability for regulatory SAPs. The FDA's emphasis on pre-specified informative priors (cited as fda2019adaptive) gives KG-DAP a clear regulatory path that data-adaptive methods like LEAP and MEM lack. The method is computationally accessible -- full LOO-CV in under 2 seconds with no MCMC -- making it practical for the extensive sensitivity analyses that responsible prior elicitation requires. The 26% MAE improvement is practically meaningful for small-trial settings. The impact is limited by: (a) single-domain evaluation (NDMM only), making it hard to assess generalizability; (b) no design OC analysis, which is what regulators and trial designers actually need to evaluate the method; (c) binary endpoints only, excluding the common continuous and time-to-event settings; (d) requirement for a curated knowledge graph, which adds an upfront burden. The method opens clear future directions (other endpoints, other disease areas, graph perturbation studies) but the current paper does not demonstrate breadth beyond NDMM.

### 7. Performance: 4/5

The empirical evidence is strong for the primary comparison: KG-DAP achieves a 25.9% MAE improvement over rMAP in LOO-CV across 35 NDMM arms, with well-calibrated 94.3% coverage versus rMAP's under-covering 88.6%. All four comparators are consistently included across all analyses. The simulation study confirms these findings under known ground truth, and importantly shows graceful degradation under the Adverse scenario (only 0.7% worse than rMAP when the graph is uninformative). The component contribution analysis cleanly isolates the contributions of graph structure (+32.6%) and power weights (+30.8%). The sensitivity analyses demonstrate a broad operating plateau for beta ([8, 30]). However, the LOO-CV is a prior predictive evaluation, not a posterior evaluation -- the method is never tested with actual current-trial data arriving to update the posterior. The simulation B=200 replications is modest by Biometrics standards (B=1000+ is typical). The MC SEs, while reported, have an incorrect prose characterization. No formal statistical tests (e.g., paired comparisons, bootstrap confidence intervals for the MAE difference) are provided to establish that the improvement is statistically significant. Performance is capped at 4 because no power/type I error analysis is conducted.

---

### Raw Score: 28/35
### Weighted Score: 40/50
  Execution (x1.0): Correctness 4 + Completeness 4 + Rigor 4 + Clarity 4 = 16/20
  Research  (x2.0): (Novelty 4 + Impact 4 + Performance 4) x 2.0 = 24/30
  Weighted total = 16 + 24 = 40/50
### Grade: B (>=34)

### Code Audit Summary

| Category | PASS | FLAG | FAIL |
|----------|------|------|------|
| Formula Correctness (1-6) | 5 | 0 | 0 |
| Parameter Consistency (7-11) | 3 | 2 | 0 |
| Methodological Completeness (12-15) | 2 | 2 | 0 |

**Formula Correctness**: RMSE (PASS), Coverage (PASS), ESS (PASS), log-predictive (PASS), posterior weights (PASS). All formulas match their labels.

**Parameter Consistency**: Hyperparameters (beta=15, w0=0.20, alpha/beta_w/gamma=0.20/0.20/0.60, n_cap=200, alpha0=beta0=1) are consistent between code and manuscript (PASS). Replicate count B=200 in simulation code matches manuscript table caption (PASS). Random seeds are deterministic per-replication but use `hash(sname)` which is non-reproducible across Python sessions (FLAG). The same seed base (42) is reused between LOO-CV and simulation but these are independent evaluations (FLAG -- minor).

**Methodological Completeness**: All 5 comparators present in all evaluations (PASS). Prior ESS reported but no posterior ESS (FLAG). No power analysis / design OC (FLAG).

### Critical Findings
None (no FAIL items).

### Flags for Human Review
1. `hash(sname)` in simulation seed construction (non-reproducible across Python sessions)
2. No posterior ESS reported
3. No design operating characteristics (power, type I error)
4. Incorrect MC SE magnitude claim in Section 5.5

### Top 3 Strengths
1. **Strong empirical evidence with comprehensive evaluation.** The 25.9% MAE improvement over rMAP is supported by LOO-CV, simulation, sensitivity analyses, contamination tests, and ablation -- all telling a consistent story. All four comparators are present in every analysis, and the simulation study tests favorable, adverse, and mixed scenarios.
2. **Honest scientific reporting.** The manuscript is unusually transparent about the failed diffusion kernel approach, the pre-specification trade-off (inability to detect prior-data conflict), and specific limitations (single domain, binary endpoints, graph quality dependence). Remark 1 and the Limitations section are exemplary.
3. **Clear regulatory path.** Full pre-specifiability with closed-form computation is a genuine practical advantage over LEAP, MEM, and commensurate priors. The connection to FDA adaptive design guidance is appropriately articulated.

### Top 3 Weaknesses (with suggested improvement)
1. **No design operating characteristics.** Regulators and trial designers need power curves and type I error assessments, not just prior predictive MAE. -- Fix: Add a design OC section computing P(go | true_rate) and P(reject H0 | true_rate) under KG-DAP vs comparators for representative scenarios.
2. **Simulation study is modest and has reproducibility issues.** B=200 replications is below Biometrics standards; `hash()` usage breaks cross-session reproducibility; no formal statistical tests for the MAE differences. -- Fix: Increase to B=1000, replace `hash(sname)` with a deterministic scenario-to-seed mapping, add paired bootstrap CIs for MAE differences.
3. **Incomplete supplementary material and dead references.** Three bib entries are uncited (lin2018, tierney1986, clemen1989), the supplementary material has only one table and one figure despite the wealth of available results, and no trial-level data table or weight heatmap is provided. -- Fix: Remove dead references, add a trial arms table (S.1), a weight heatmap for a representative arm (S.2), and move the contamination detailed table to supplement.

### Fixable Issues (for paper-fixer)

1. **Incorrect MC SE claim** (Section 5.5): "All Monte Carlo standard errors are below 0.02% of the corresponding effect sizes" is wrong by ~20x. Replace with accurate statement, e.g., "All Monte Carlo standard errors are below 0.0002 (Table 4)."
2. **Dead bibliography entries**: Remove lin2018borrowing, tierney1986accurate, and clemen1989combining from references.bib (they are never cited in the manuscript).
3. **Inconsistent terminology**: The text alternates between "component contribution analysis" (Section 5.4 header) and "ablation" (Figure 4 caption, abstract). Standardize to one term throughout.
4. **Minor rounding inconsistency in ablation table**: Drug-Only delta is reported as +7.3% in the manuscript table but the validation report says +7.2%. Use consistent rounding (the computed value is 7.24%).
5. **Missing posterior ESS reporting**: Add posterior ESS to the LOO-CV table or a supplementary table (the code already computes it via posterior_update).
6. **Supplementary material**: Add trial arms summary table (trial, n, y, p_hat) and a representative weight heatmap figure to supplement.
7. **Simulation seed reproducibility**: Replace `hash(sname)` in run_simulation_v3.py line 81 with a deterministic mapping (e.g., `{'Favorable': 0, 'Adverse': 1, 'Mixed': 2}[sname]`).

### Verdict

This is a solid methods paper with a clear contribution (knowledge-graph-derived borrowing weights with full pre-specifiability), comprehensive primary evaluation, and honest scientific reporting. The main gaps preventing an A grade are the absence of design operating characteristics (which Biometrics reviewers will likely require), the modest simulation study, and the incomplete supplementary material. With design OC added, simulation strengthened (B>=1000, reproducible seeds, formal tests), and the fixable issues addressed, this manuscript has a realistic path to acceptance at Biometrics after one round of major revision. Current assessment: **revise and resubmit** at Biometrics standards.
