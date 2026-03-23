# Literature Briefing for KG-DAP Manuscript

## Executive Summary

The clinical trial borrowing landscape is dominated by methods that either treat all historical trials as exchangeable (rMAP, Schmidli et al. 2014), use data-driven adaptation that cannot be pre-specified (LEAP, Alt et al. 2024; commensurate priors, Hobbs et al. 2011; MEM, Kaizer et al. 2018), or rely on binary exchangeability decisions (iMEM/dMEM). KG-DAP fills a specific gap: it is the first method to derive trial-specific, continuous borrowing weights from an external knowledge graph while maintaining full pre-specifiability and closed-form computation. The key novelty is importing graph kernel techniques from machine learning (Kondor and Lafferty 2002; Smola and Kondor 2003) into the Bayesian historical borrowing framework, operationalized as power-sharpened composite similarity weights for a finite Beta mixture prior.

## Comparative Table

| Aspect | rMAP (Schmidli 2014) | LEAP (Alt 2024) | MEM (Kaizer 2018) | Commensurate (Hobbs 2011) | Power Prior (Ibrahim 2000) | **KG-DAP (Ours)** |
|--------|---------------------|-----------------|-------------------|--------------------------|---------------------------|-------------------|
| Prior specification | Partially pre-specifiable (tau estimated) | Not pre-specifiable | Not pre-specifiable | Not pre-specifiable | Pre-specifiable if a_0 fixed | **Fully pre-specifiable** |
| Borrowing mechanism | Global heterogeneity tau | Latent exchangeability indicators | Binary exchangeability partitions | Per-trial shrinkage tau_k | Global discount a_0 | **Graph-derived power weights** |
| Trial-specific weights | No (global tau) | Yes (posterior) | Yes (posterior) | Yes (data-driven) | No (global) | **Yes (graph-derived)** |
| Uses external structure | No | No | No | No | No | **Yes (knowledge graph)** |
| Computation | MCMC | MCMC | 2^K enumeration | MCMC | MCMC or closed-form | **Closed-form** |
| Robustness mechanism | Vague mixture | Implicit | Model averaging | Implicit in tau | a_0 -> 0 | **Explicit vague component w_0** |
| Theoretical results | Asymptotic properties | Exchangeability consistency | Model selection consistency | Posterior consistency | Propriety conditions | **Conjugate closure, bounded ESS, monotonic borrowing** |
| Evaluation types | Simulation + real data | Simulation + real data | Simulation + real data | Simulation + real data | Simulation + real data | **LOO-CV + sensitivity + contamination + ablation** |

## Gap Analysis

### Gap 1: Simulation Study with Known Ground Truth
- **Gap**: Phase 2 validated on real data LOO-CV only; reference papers include simulation studies with known parameters
- **Evidence**: Schmidli et al. (2014), Alt et al. (2024) both have extensive simulation sections
- **Priority**: Critical (Biometrics expects simulation evidence)
- **Resolution**: Paper Modeler will implement simulation study (3+ scenarios)

### Gap 2: Operating Characteristics (Type I Error, Power)
- **Gap**: No formal operating characteristics (type I error control, power) reported
- **Evidence**: Schmidli et al. (2014) reports operating characteristics under various scenarios
- **Priority**: Important (but can be addressed via the simulation study)
- **Resolution**: Include in simulation study or discuss as future work

### Gap 3: Knowledge Graph Visualization
- **Gap**: No visual representation of the knowledge graph structure
- **Evidence**: Graph-based methods typically show the graph
- **Priority**: Nice-to-have (can include as supplementary figure)
- **Resolution**: Could generate from build_kg.py if time permits

### Gap 4: Comparison with LEAP/MEM on Same Dataset
- **Gap**: Only compare to rMAP, Uniform, Pooled, EqualWeight — not to LEAP or MEM directly
- **Evidence**: LEAP and MEM require patient-level data or MCMC infrastructure not available
- **Priority**: Important (discuss as limitation — cannot run LEAP without patient-level data for all trials)
- **Resolution**: Discuss in limitations; LEAP requires IPD, MEM is O(2^K) infeasible for K=34

## Quality Benchmarks (from Biometrics standards)

Based on the methodology rationale's analysis of reference papers:
- **Theorems/Propositions**: At least 3-4 formal results with proofs (we have 5 properties)
- **Evaluation types**: At least 3 (simulation, real data, sensitivity) — we have 5 (LOO-CV, sensitivity, contamination, ablation, simulation)
- **Tables**: 4-6 in main text, additional in supplementary
- **Figures**: 4-6 publication quality
- **Section structure**: Introduction, Methods, Theoretical Properties, Evaluation/Application, Discussion
- **Page count**: 20-25 pages double-spaced main text

## Related Work Draft

The problem of borrowing information from historical clinical trials has generated a rich methodological literature spanning Bayesian hierarchical models, mixture priors, and adaptive borrowing mechanisms.

**Meta-Analytic-Predictive Priors.** The robust meta-analytic-predictive (rMAP) prior \citep{schmidli2014robust} has become the standard approach in pharmaceutical applications. rMAP fits a normal-normal hierarchical model to logit-transformed historical rates, derives a predictive distribution for the new trial, and robustifies by mixing with a vague component. While widely adopted and regulatory-accepted, rMAP treats all historical trials as exchangeable up to a global heterogeneity parameter $\tau$ --- it cannot assign differential borrowing weights based on trial-specific structural similarity.

**Exchangeability-Based Methods.** The Multisource Exchangeability Model (MEM; \citealt{kaizer2018bayesian}) enumerates all $2^K$ binary exchangeability partitions and computes posterior-weighted mixtures. Extensions include the independent MEM (iMEM; \citealt{kaizer2021bayesian}) and dynamic MEM (dMEM; \citealt{ji2022dynamic}). These methods produce trial-specific borrowing decisions but at $O(2^K)$ computational cost, making them intractable for $K > 20$. Moreover, exchangeability decisions are posterior quantities that cannot be pre-specified.

**Latent Exchangeability Methods.** LEAP (Latent Exchangeability with Aggregation of Priors; \citealt{alt2024leap}) constructs trial-specific priors via latent exchangeability indicators, representing the closest conceptual predecessor to KG-DAP. LEAP's weights adapt to observed outcome data, providing conflict detection that KG-DAP lacks. However, LEAP requires MCMC computation and cannot be fully pre-specified --- a critical limitation for regulatory Statistical Analysis Plans.

**Commensurate and Power Priors.** Commensurate priors \citep{hobbs2011hierarchical} place continuous shrinkage on the difference between current and historical trial parameters. Power priors \citep{ibrahim2000power} discount historical likelihoods by a scalar $a_0$. The propensity-score-integrated power prior \citep{wang2020propensity} extends this with covariate adjustment. All these methods adapt borrowing from outcome data, precluding full pre-specification.

**Graph-Based Methods in Biostatistics.** Knowledge graphs and network structures have been used for drug repositioning \citep{himmelstein2017systematic} and clinical trial matching \citep{zhang2020deepenroll}, but their application to informative prior construction is novel. The closest work is the conditional autoregressive (CAR) prior approach, which uses adjacency structure to induce spatial correlation among trial parameters but requires MCMC and estimates the spatial correlation from data.

**KG-DAP's Position.** KG-DAP uniquely occupies the intersection of three desiderata: (i) trial-specific continuous weights derived from external structural information, (ii) full pre-specifiability of all prior components, and (iii) closed-form computation via the conjugate Beta mixture. No existing method satisfies all three simultaneously.

## Action Items for Writer

### Narrative to Build
The story is about the **direction** of borrowing mattering more than the **amount**: "Even modest ESS (~6) produces large MAE improvements (26%) when the borrowing is properly directed by graph structure." Frame KG-DAP as importing graph kernel ideas from ML into Bayesian clinical trial statistics.

### Claims Needing Citations
1. rMAP is the most widely used method in pharma → cite Schmidli 2014, FDA guidance 2019
2. Pre-specifiability matters for regulatory SAPs → cite FDA adaptive design guidance
3. MEM is computationally intractable for large K → cite Kaizer 2018 (O(2^K))
4. Diffusion kernels on graphs → cite Kondor & Lafferty 2002, Smola & Kondor 2003
5. ESS computation via moment-matching → cite Morita et al. 2008
6. Conjugacy of finite mixtures → cite Diaconis & Ylvisaker 1979
7. DerSimonian-Laird estimator for rMAP → cite DerSimonian & Laird 1986
8. Forecast combination principles → cite Clemen 1989

### Sections to Emphasize
1. **Methods**: Full mathematical specification with enough detail to reproduce
2. **Theoretical Properties**: All 5 properties as formal propositions
3. **LOO-CV Results**: Primary evidence table with ALL 5 comparators
4. **Ablation Study**: Critical for demonstrating that both graph structure AND power weights contribute
5. **Sensitivity Analysis**: Broad plateau shows robustness to hyperparameter choice

### Quality Benchmarks to Meet
- At least 4 formal propositions with proofs
- LOO-CV + simulation + sensitivity + contamination + ablation evaluations
- All 5 comparators in ALL results tables
- Honest discussion of limitations (at least 5 specific limitations)
- At least 10 bibliography entries from established journals
