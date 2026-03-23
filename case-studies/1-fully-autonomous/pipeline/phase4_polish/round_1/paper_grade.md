## Paper Assessment

### Manuscript: Knowledge-Graph-Driven Adaptive Prior for Bayesian Borrowing in Clinical Trials
### Date: 2026-03-23
### Round: 1 (previous score: 40/50)

---

### 1. Correctness: 4/5

All mathematical derivations remain correct. The five propositions and their proofs are valid. Numbers in the LOO-CV table (Table 2) match the validated CSV data to 4 significant figures: KG-DAP MAE=0.1257 (CSV: 0.12567), rMAP MAE=0.1697 (CSV: 0.16967), coverage 94.3% (CSV: 0.9429), etc. The simulation table (Table 4, B=1000) matches the `simulation_summary_B1000.csv` file exactly. The design OC table (Table 6) matches `design_oc_representative.csv` to 3 decimal places. The MC SE claim has been correctly fixed: "All Monte Carlo standard errors are below 0.0002 (Table 4)" is accurate. However, a new correctness issue was introduced by the fix: the Drug-Only delta in the component contribution analysis table (Table 3) now reads +7.2%, but the validated data gives (0.1348 - 0.1257)/0.1257 = 7.2523%, which rounds to +7.3% at one decimal place under standard rounding. The "fix" changed a correct value to an incorrect one. Proposition 2's "proof sketch" remains heuristic (the variance bound argument does not rigorously establish the stated ESS inequality), but the label "proof sketch" is appropriate. Two bibliography entries (ji2022dynamic, kaizer2021bayesian) are present in references.bib but never cited in the manuscript -- while the original three dead entries were removed, these two became orphaned, apparently because their citations were dropped during editing.

### 2. Completeness: 5/5

Round 1 addresses the two primary completeness gaps identified in Round 0. First, Section 5.6 adds design operating characteristics: power curves via go/no-go decision analysis across 10 true rates, three methods, six representative arms, B=5,000 replications per arm, with a clear decision rule (P(theta > 0.30 | data) > 0.80). This is exactly the analysis that Biometrics reviewers and regulators need. Second, posterior ESS is now reported in the LOO-CV table alongside prior ESS, addressing the prior-only ESS reporting gap. The manuscript now covers all standard and extended components: LOO-CV (primary), sensitivity analyses (beta, w0), contamination robustness, component contribution analysis, simulation study with bootstrap CIs, design operating characteristics, and supplementary material with trial arms table (S1), borrowing weights table (S2), and full sensitivity table (S3). All four comparators appear in all primary analyses. Six figures are referenced and present. This comprehensive coverage merits 5/5 under the rubric ("exhaustive: + design OC + multiple robustness checks + ablation").

### 3. Rigor: 4/5

Uncertainty quantification is now substantially improved. The simulation study uses B=1,000 (up from 200), meeting Biometrics standards. Paired bootstrap 95% CIs for the KG-DAP vs. rMAP MAE difference are reported for all three scenarios, establishing statistical significance in Favorable and Mixed (CIs exclude zero) and non-significance in Adverse. The seed construction uses a deterministic mapping (SCENARIO_SEEDS dict with formula `42 + b*997 + scenario_seed*9973`), fully reproducible across sessions -- the hash() issue from Round 0 is fixed. MC standard errors are accurately reported. Sensitivity analyses cover both key hyperparameters. The type I error discussion in the design OC section (Section 5.6) is honest, noting the inflation and explaining that threshold recalibration is standard practice. However, the bootstrap CIs appear degenerate at 3 decimal places (e.g., [-0.044, -0.044] for Favorable), which is a presentation issue suggesting extreme precision that may mislead readers about the actual width. The full precision data shows [-0.0440, -0.0438], a real but very narrow interval. Reporting at 4 decimal places would be more informative. Proposition 2 remains a "proof sketch" without a complete proof. Prior sensitivity to the composite similarity weights (alpha, beta_w, gamma) is not systematically explored, though it is acknowledged as secondary.

### 4. Clarity: 4/5

The manuscript remains well-organized with consistent notation and logical flow. The terminology has been standardized: "component contribution analysis" is used consistently in all user-visible text (abstract, section header, table captions, figure captions, discussion). Internal LaTeX labels use "ablation" but this is invisible to readers. The new Section 5.6 (Design OC) is well-motivated and clearly presented, with the decision rule stated precisely before the results. The paragraph interpreting the type I error inflation is a strength -- it explains the trade-off honestly and notes the standard recalibration solution. The posterior ESS column is a useful addition to Table 2, and the accompanying paragraph (lines 376) provides a clean interpretation linking it to Proposition 3. Minor issues remain: (a) the bootstrap CIs in Table 4 showing [-0.044, -0.044] look like point estimates, not intervals; (b) two dead bib entries (ji2022dynamic, kaizer2021bayesian) remain; (c) the abstract is quite long (one dense paragraph spanning ~16 lines of text) and could be more concise.

### 5. Novelty: 4/5

Unchanged from Round 0. KG-DAP represents a genuine cross-disciplinary contribution: importing knowledge graph structure into Bayesian borrowing weight construction. The combination of power-sharpened KG-derived weights, finite Beta mixture architecture, full pre-specifiability, and closed-form computation is non-obvious and addresses a real gap demonstrated by Table 1. The key insight ("what you borrow from matters more than how much") is well-supported empirically. The novelty is clearly above "meaningful extension" (score 3) because the cross-domain bridge from knowledge graphs to Bayesian priors requires insights from multiple fields and the resulting framework opens new research directions. However, the power-sharpened weights are a relatively simple softmax function, not the more theoretically novel diffusion kernel originally proposed. Score remains 4/5.

### 6. Impact: 4/5

The addition of design operating characteristics (Section 5.6) substantially strengthens the impact argument compared to Round 0. Regulators and trial designers can now directly see the method's effect on decision-making: KG-DAP increases the probability of a correct go decision by 7-12 percentage points at clinically relevant alternatives (theta = 0.35-0.40) compared to no borrowing. The type I error discussion demonstrates responsible reporting. Full pre-specifiability for SAPs, closed-form computation enabling rapid sensitivity analysis, and the quantified 26% MAE improvement make a clear case for practical value. The method remains limited by single-domain evaluation (NDMM only), binary endpoints only, and the requirement for a curated knowledge graph. The design OC analysis, while valuable, uses an average over only 6 representative arms rather than all 35, which somewhat limits its generalizability. Impact remains 4/5 -- the design OC addition closes the gap to practice but does not yet demonstrate regulatory uptake or multi-domain validation.

### 7. Performance: 4/5

The empirical evidence is strengthened by two additions. First, the B=1,000 simulation with paired bootstrap CIs now provides formal statistical significance for the MAE differences in Favorable and Mixed scenarios, addressing the Round 0 gap of "no formal statistical tests." The bootstrap CIs are narrow (reflecting the tight Monte Carlo variance at B=1,000), confirming that the 25.9% Favorable improvement and 16.2% Mixed improvement are robust to sampling variability. Second, the design OC analysis translates the MAE improvements into decision-relevant quantities (power curves). Under the Adverse scenario, KG-DAP gracefully degrades (0.8% worse than rMAP), confirming the safety property. The Uniform method slightly outperforms in the simulation's Favorable scenario for MAE (0.1608 vs KG-DAP's 0.1258 -- wait, KG-DAP is better), confirming KG-DAP's advantage across all non-adversarial conditions. All four comparators are consistently present. Performance remains at 4/5 because: (a) no posterior evaluation is performed (only prior predictive), (b) the degenerate-looking bootstrap CIs raise questions about their informativeness, and (c) no formal power analysis (frequentist sample size) is provided despite the design OC section.

---

### Raw Score: 29/35
### Weighted Score: 41/50
  Execution (x1.0): Correctness 4 + Completeness 5 + Rigor 4 + Clarity 4 = 17/20
  Research  (x2.0): (Novelty 4 + Impact 4 + Performance 4) x 2.0 = 24/30
  Weighted total = 17 + 24 = 41/50
### Grade: B (>=34)

### Code Audit Summary

| Category | PASS | FLAG | FAIL |
|----------|------|------|------|
| Formula Correctness (1-6) | 5 | 0 | 0 |
| Parameter Consistency (7-11) | 4 | 1 | 0 |
| Methodological Completeness (12-15) | 4 | 0 | 0 |

**Formula Correctness**: RMSE (PASS -- sqrt of mean of squared errors in validated code), Coverage (PASS -- non-strict inequality), ESS (PASS -- mu*(1-mu)/var - 1 matches Morita 2008), log-predictive (PASS), posterior weights (PASS -- log-sum-exp normalization).

**Parameter Consistency**: Hyperparameters match between code and manuscript (PASS). B=1,000 in simulation code matches manuscript table caption (PASS). Deterministic seeds via SCENARIO_SEEDS dict (PASS -- replaces hash()). Two dead bib entries (FLAG -- ji2022dynamic, kaizer2021bayesian in bib but not cited). Decision threshold (theta_0=0.30, gamma_go=0.80) consistent between Section 5.6 text and table caption (PASS).

**Methodological Completeness**: All 5 comparators present in all evaluations (PASS). Both prior and posterior ESS now reported (PASS). Power analysis via design OC section (PASS). Paired bootstrap CIs for method comparisons (PASS).

### Critical Findings
None (no FAIL items).

### Flags for Human Review
1. Drug-Only delta changed from correct +7.3% to incorrect +7.2% (the "fix" introduced an error)
2. Two dead bib entries: ji2022dynamic, kaizer2021bayesian
3. Bootstrap CIs display as degenerate at 3dp (e.g., [-0.044, -0.044]) -- consider 4dp

### Top 3 Strengths
1. **Design operating characteristics close the practice gap.** Section 5.6 directly translates MAE improvements into power curves and go/no-go decision probabilities, making the method's value concrete for regulators and trial designers. The honest discussion of type I error inflation is exemplary.
2. **Statistically rigorous simulation.** The B=1,000 simulation with deterministic seeds and paired bootstrap CIs establishes formal statistical significance, meeting Biometrics standards. The seed fix ensures full cross-session reproducibility.
3. **Comprehensive and internally consistent evaluation.** All table numbers match validated CSV data to 4 significant figures. LOO-CV, sensitivity, contamination, component contribution analysis, simulation, and design OC tell a coherent story with five consistent comparators throughout.

### Top 3 Weaknesses (with suggested improvement)
1. **Drug-Only delta rounding error introduced by fix.** The table shows +7.2% but the validated data rounds to +7.3%. -- Fix: Change "+7.2" back to "+7.3" in the component contribution analysis table (or compute from the validated CSV directly: (0.1348 - 0.1257)/0.1257 * 100 = 7.25, which rounds to 7.3%).
2. **Two dead bibliography entries.** ji2022dynamic and kaizer2021bayesian are in references.bib but not cited in the manuscript text. -- Fix: Either cite them in the Related Work section (they are relevant: iMEM and dMEM) or remove them from references.bib.
3. **Bootstrap CIs appear degenerate at 3dp.** Table 4 shows intervals like [-0.044, -0.044] which look like point estimates, reducing the informativeness of the uncertainty quantification. -- Fix: Report CIs at 4 decimal places ([-0.0440, -0.0438]) or report them as [-0.044, -0.044] with a note that the width is <0.001, or report the CI width explicitly.

### Fixable Issues (for paper-fixer)
1. **Drug-Only delta rounding**: Change +7.2% to +7.3% in Table 3 and the corresponding bullet point (line 459).
2. **Dead bib entries**: Either cite ji2022dynamic and kaizer2021bayesian in the Related Work section (they were mentioned in the Round 0 manuscript) or remove them from references.bib.
3. **Bootstrap CI presentation**: Report CIs at 4 decimal places in Table 4 for clarity.
4. **Proposition 2 label**: Consider downgrading from "Proof sketch" to "Heuristic argument" or providing a complete proof using Jensen's inequality on the precision functional.

### Verdict

Round 1 achieves a meaningful improvement over Round 0 (41/50 vs 40/50), driven primarily by the Completeness dimension moving from 4 to 5 with the addition of design operating characteristics and posterior ESS reporting. The B=1,000 simulation with bootstrap CIs and deterministic seeds addresses the three most important Round 0 flags (no formal tests, modest B, hash()-based seeds). The Drug-Only delta rounding error is a minor regression introduced by the fix process. Two orphaned bib entries are a housekeeping issue. The manuscript is now substantially closer to Biometrics submission quality, with all seven evaluation components (LOO-CV, sensitivity, contamination, component contribution, simulation, design OC, supplementary material) present and internally consistent. The remaining gap to the A threshold (42/50) is narrow: fixing the two minor correctness issues and improving the bootstrap CI presentation would address the residual weaknesses. Current assessment: **strong revise-and-resubmit** at Biometrics standards, approaching acceptability.
