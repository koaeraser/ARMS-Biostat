# Design Rationale for KG-DAP (Knowledge-Graph-Driven Adaptive Prior)

## 1. Why a Beta Mixture?

### The choice

KG-DAP represents the informative prior as a finite (K+1)-component Beta mixture on the probability scale, rather than using a hierarchical normal model on the logit scale (as in rMAP), a conditional autoregressive (CAR) prior on the logit scale (as in the predecessor KG-CAR), or a nonparametric prior.

### Alternatives considered

**Hierarchical normal on logit scale (rMAP-style).** The standard approach in pharma (Schmidli et al. 2014) places a normal-normal hierarchical model on logit-transformed historical rates: logit(p_k) ~ N(mu, tau^2). This produces a MAP predictive distribution for the new trial, which is then mixed with a vague component. The approach is well-understood and has regulatory acceptance. However, it treats all historical trials as exchangeable up to a random heterogeneity parameter tau — there is no mechanism for graded, trial-specific borrowing based on structural similarity. The between-trial heterogeneity tau is estimated from the data, making the prior partially data-dependent and harder to pre-specify in a Statistical Analysis Plan.

**Conditional autoregressive (CAR) prior (KG-CAR, predecessor method).** The previous version of this project used a joint multivariate normal prior on logit-response rates with CAR spatial structure derived from the knowledge graph: phi | sigma ~ N(0, sigma^2 (D_w - rho W)^{-1}). This incorporated graph structure but required MCMC or multivariate Laplace approximation, and critically, the spatial correlation parameter rho was estimated from the data. This made the prior data-dependent — the borrowing weights were not fully determined before observing outcomes. Additionally, extending the CAR prior to a new trial not yet in the graph required augmenting the joint distribution, which was algebraically awkward.

**Dirichlet process or nonparametric mixture.** A Dirichlet process mixture of Beta distributions would allow the number of mixture components to be learned from data, providing flexibility without specifying K. This was rejected for two reasons: (i) it introduces MCMC dependence that eliminates the closed-form computational advantage; (ii) the number of components is already known — it equals the number of historical trial arms plus one vague component, with the graph determining the weights rather than the data determining the clustering.

**Direct Beta prior via moment-matching.** One could collapse all historical information into a single Beta(a, b) distribution via moment-matching, skipping the mixture entirely. This was rejected because it destroys the trial-level resolution: if historical trials have heterogeneous rates (ranging from 0.07 to 0.79 in our NDMM data), a single Beta cannot represent the bimodal or multimodal structure. The mixture preserves the ability to identify which historical trials are concordant with the new trial's data after the posterior update.

### Why the Beta mixture wins

The Beta mixture was chosen because it uniquely satisfies four requirements simultaneously:

1. **Conjugate closure.** Each Beta component is conjugate to the Binomial likelihood. The posterior is therefore also a Beta mixture — with updated parameters within each component and re-normalized weights across components. No MCMC, no Laplace approximation, no numerical integration. This is the single most important computational property: it makes the entire LOO cross-validation over 35 arms trivially fast (seconds, not minutes).

2. **Fixed, pre-specifiable weights.** The mixture weights omega_k are derived from the knowledge graph via the diffusion kernel, without reference to outcome data. This means the complete prior — all K+1 component parameters and all K+1 weights — can be written into a regulatory Statistical Analysis Plan before unblinding. The hierarchical normal and CAR approaches cannot make this claim because they estimate between-trial variance from the data.

3. **Trial-level interpretability.** Each mixture component corresponds to a specific historical trial arm. After the Bayesian update, the posterior mixture weights reveal which historical trials were "most believed" given the observed data. This component-level interpretability is absent from the hierarchical normal model (where individual trial contributions are absorbed into the hierarchical parameters) and from the CAR model (where the correlation structure is global, not arm-specific).

4. **Natural robustness scaling.** Adding the vague Beta(1,1) component with fixed weight w_0 provides a floor of protection that is transparent and tunable. In the hierarchical normal model, robustness depends on the tail behavior of the MAP predictive distribution; in the CAR model, there is no explicit robustness mechanism.

### What was sacrificed

The Beta mixture does not model borrowing of the between-trial variance (heterogeneity), which the hierarchical normal model does naturally. In the hierarchical normal model, tau provides a data-driven estimate of how much the historical trials disagree, and the MAP predictive distribution automatically widens when tau is large. In KG-DAP, the width of the prior mixture is determined by the fixed spread of the historical Beta components and their weights — there is no separate "heterogeneity" parameter. This is partially compensated by the vague component (which widens the mixture when w_0 is non-negligible) and by the diffusion kernel (which downweights dissimilar trials), but the compensation is structural rather than learned.

---

## 2. Why Diffusion Kernel?

### The choice

KG-DAP converts the composite similarity matrix S into borrowing weights via a graph diffusion kernel: K_beta = exp(-beta * L), where L is the normalized graph Laplacian derived from S. The raw affinity of the current trial to each historical trial is read from the corresponding row of K_beta and then normalized to form the weight vector omega.

### Alternatives considered

**Raw similarity as weights.** The simplest approach: set omega_k proportional to S_{c,k} (the composite similarity between the current trial c and historical trial k). This was rejected because it ignores the global structure of the graph. Two trials may have identical pairwise similarity to c but very different positions in the graph — one may be in a dense cluster of similar trials (providing redundant information) while the other may be an isolated bridge to a different region (providing unique information). Raw similarity treats them identically; diffusion does not.

**k-nearest-neighbor selection.** Select the k most similar trials and weight them equally (or proportionally to similarity), discarding the rest. This was rejected because it introduces a discrete threshold — the k-th and (k+1)-th most similar trials may have nearly identical similarity scores but receive completely different treatment (full weight vs. zero). The diffusion kernel produces a smooth, continuous weight function with no discontinuity.

**Graph neural network (GNN) embeddings.** Learn trial-level embeddings via a GNN (e.g., GraphSAGE, GCN) and use embedding distance to derive weights. This was rejected for three reasons: (i) with only 35 trial arms and a 56-node graph, there is far too little training data to learn meaningful embeddings — the GNN would overfit severely; (ii) a learned embedding makes the weights data-dependent, violating the pre-specifiability requirement; (iii) the embedding is a black box, eliminating the interpretability of the weight derivation.

**Personalized PageRank.** Use the PageRank vector with restart at node c as the weight vector. This is conceptually similar to diffusion but assumes a random walk model (multiplicative propagation along edges) rather than a heat equation model (additive diffusion). PageRank was considered a viable alternative but was not selected because: (i) the diffusion kernel has a cleaner spectral interpretation (exp(-beta * lambda_i) for each eigenvalue), making the effect of beta on borrowing concentration transparent; (ii) the diffusion kernel is symmetric positive definite, ensuring that K_beta(i,j) = K_beta(j,i), which is a natural requirement for similarity-based borrowing; (iii) the literature on diffusion kernels for graphs (Kondor and Lafferty 2002; Smola and Kondor 2003) provides established theoretical properties.

### Why diffusion wins

The diffusion kernel was chosen because it provides three specific capabilities that the alternatives lack:

1. **Higher-order similarity propagation.** Consider two trials that share no drugs but whose drugs target the same molecular pathway via different mechanisms. Their drug Jaccard J_drug is zero, yet they may be biologically similar. In the knowledge graph, they are connected through a two-hop path: trial_A -> drug_X -> target_Z <- drug_Y <- trial_B. The diffusion kernel propagates similarity through these intermediate nodes, producing a nonzero kernel value K_beta(A, B) even though the direct similarity is low. This is the core advantage of diffusion over raw similarity: it captures transitive, higher-order relationships encoded in the graph topology.

2. **Single interpretable hyperparameter.** The diffusion parameter beta controls a continuous spectrum from uniform borrowing (beta -> 0, all trials receive equal weight) to nearest-neighbor borrowing (beta -> infinity, only the most similar trial receives weight). Moderate values of beta produce graded borrowing that reflects the graph's multi-scale structure. This single parameter replaces what would otherwise be multiple discrete choices (how many neighbors? what similarity threshold? how to handle indirect connections?) with a single continuous tuning parameter whose effect is monotonic and interpretable.

3. **Spectral foundation.** The diffusion kernel is computed via the eigendecomposition of the graph Laplacian: K_beta = U diag(exp(-beta * lambda_i)) U^T. This spectral representation means that beta acts as a frequency filter: low eigenvalues (smooth graph modes) are preserved at all beta values, while high eigenvalues (rough graph modes corresponding to local fluctuations) are increasingly suppressed as beta grows. Borrowing weights derived from this kernel therefore favor structurally smooth patterns of similarity — exactly what one expects from a biomedical knowledge graph where drug-target-population relationships vary smoothly across trials.

### What was sacrificed

The diffusion kernel treats the graph as undirected (the Laplacian is symmetric), even though the underlying knowledge graph is a directed multigraph (MultiDiGraph). Information about edge direction (e.g., "trial uses drug" vs. "drug is used by trial") is lost when constructing the symmetric similarity matrix S. This is an acceptable simplification for our application because the relevant relationships are inherently symmetric: if trial A and trial B share a drug, the biological relevance is symmetric regardless of edge direction. However, for knowledge graphs with asymmetric relationships (e.g., causal graphs where direction matters), the diffusion kernel would need to be replaced with an asymmetric propagation mechanism.

---

## 3. Why Pre-Specified Weights?

### The choice

The borrowing weights omega_k in KG-DAP are fully determined by the knowledge graph structure and the hyperparameters (alpha, beta_w, gamma, beta), with no dependence on the observed outcome data (n_k, y_k) or (n_c, y_c). This is a deliberate design decision, not an oversight.

### The regulatory argument

In regulated clinical trials, the statistical analysis plan (SAP) must be finalized and registered before the trial data is unblinded. Any aspect of the analysis that depends on the outcome data — including the choice of prior — raises concerns about data dredging, post-hoc rationalization, and inflated type I error. Pre-specifiable priors avoid this concern entirely: the prior can be written into the SAP with specific numerical values for all mixture weights and component parameters, locked before the first patient is enrolled.

This is not merely an academic point. Regulatory agencies (FDA, EMA) have expressed increasing openness to Bayesian designs with informative priors, but only when the prior can be justified prospectively. The FDA's 2019 guidance on adaptive designs explicitly distinguishes between pre-specified informative priors (acceptable) and data-adaptive borrowing mechanisms (requiring additional justification for type I error control). KG-DAP's pre-specifiability places it firmly in the acceptable category.

### The practical argument

Data-adaptive methods (commensurate priors, elastic priors, LEAP) have a more subtle problem: they require the practitioner to trust a black-box adaptation mechanism. When the commensurate prior shrinks the historical data away because it detects "conflict," the practitioner may not understand why — is it genuine heterogeneity, or a chance fluctuation in a small current trial? When LEAP's exchangeability probabilities shift, the practitioner cannot explain to a clinical team why the prior changed. KG-DAP's weights, by contrast, are derived from the knowledge graph and can be presented to a clinical team as: "We borrow 35% from Trial X because it uses the same drugs in the same population, 20% from Trial Y because it shares molecular targets, and 10% from Trial Z because it has a similar patient population." This transparency supports clinical buy-in, which is essential for actual adoption.

### What was sacrificed

The most significant cost of pre-specification is the inability to downweight a historical trial that turns out to have a very different outcome from the current trial, even when the graph says they should be similar. If the knowledge graph is wrong — if two trials that share drugs, targets, and population characteristics nonetheless have very different response rates due to unmeasured factors (dosing schedule, study conduct quality, geographic differences) — KG-DAP will borrow from the misleading trial at the pre-specified weight. The vague component w_0 provides partial protection (bounding the damage), but it cannot provide the targeted protection of a data-adaptive method that detects and responds to the specific conflict.

This trade-off is fundamental and cannot be resolved within the pre-specification framework. The rationale for accepting it is that the knowledge graph encodes genuine biomedical structure — drugs really do target specific proteins, and trials in the same disease with the same drugs and populations really do tend to have similar outcomes. When this assumption holds (which we expect in the majority of cases), pre-specified weights perform well. When it fails, the robustness component limits the damage. The alternative — data-adaptive borrowing — performs better in the conflict case but cannot be pre-specified, which limits its regulatory acceptability and practitioner trust.

### Comparison with data-adaptive alternatives

| Property | KG-DAP (pre-specified) | Commensurate prior | LEAP | Power prior |
|----------|----------------------|-------------------|------|-------------|
| Weights depend on outcomes? | No | Yes (tau_k adapted) | Yes (exchangeability prob.) | Yes (a_0 adapted) |
| Pre-specifiable in SAP? | Yes | No | No | Only if a_0 fixed |
| Conflict detection? | No (vague component only) | Yes (per-trial) | Yes (pairwise) | Yes (global) |
| Computational cost | Closed-form | MCMC | MCMC | MCMC or closed-form |
| Interpretability of weights | Direct (graph-derived) | Indirect (posterior of tau) | Indirect (exchangeability) | Direct (but scalar) |

---

## 4. Cross-Disciplinary Inspirations

KG-DAP combines ideas from at least four fields outside traditional clinical trial statistics:

### 4.1 Graph kernels from machine learning

The diffusion kernel K_beta = exp(-beta * L) was introduced by Kondor and Lafferty (2002, "Diffusion kernels on graphs and other discrete structures," ICML) as a positive definite kernel on graphs for use in SVMs and kernel methods. Smola and Kondor (2003, "Kernels and regularization on graphs," COLT) established its regularization-theoretic properties. The key insight imported from this literature is that the matrix exponential of the graph Laplacian provides a principled way to propagate local similarity information through the graph's global topology — exactly what is needed to convert a knowledge graph into borrowing weights.

In the machine learning context, diffusion kernels are used for semi-supervised classification (propagating labels from labeled to unlabeled nodes) and graph-based clustering. KG-DAP repurposes this same mathematical object for a different task: propagating borrowing affinity from a target trial to historical trials based on the knowledge graph's structure. The mathematical machinery is identical; the interpretation is new.

### 4.2 Sensor fusion with heterogeneous reliability

The abstract structure of KG-DAP — combining multiple information sources with source-specific weights, plus a "no information" fallback — mirrors the architecture of multi-sensor fusion systems. In sensor fusion (e.g., Kalman filtering with multiple sensors), each sensor provides a noisy estimate of the same quantity, and the fusion algorithm assigns weights based on each sensor's known reliability. When a sensor's reliability is unknown or suspected to be degraded, a robust fusion algorithm includes a "default" estimate (analogous to the vague component) that prevents catastrophic failure.

The specific mapping is: each historical trial arm is a "sensor" reporting an estimate of the response rate that the current trial would observe. The knowledge graph provides a structural reliability estimate for each sensor (how similar its context is to the current trial). The vague component is the "default sensor" — always available, always uninformative, providing a safety net.

This sensor fusion framing influenced the decision to use a weighted mixture (the standard fusion architecture) rather than a hierarchical model (which would be more natural from a purely Bayesian perspective). The mixture architecture makes the source-specific contributions explicit and modular, which is the hallmark of fusion systems.

### 4.3 Heat diffusion on manifolds from physics and applied mathematics

The diffusion kernel exp(-beta * L) is the fundamental solution of the heat equation on a graph. This provides a physical interpretation: if each historical trial is a point source of "heat" (information), then the kernel value K_beta(c, k) represents the amount of heat that reaches the current trial c from source k after time beta. Trials that are nearby on the graph (in the diffusion sense) contribute more heat; trials that are far away contribute less. The diffusion time beta controls how far information propagates: short times produce localized borrowing, long times produce diffuse borrowing.

This physical interpretation influenced the choice of beta as the primary tuning parameter: it has a natural scale (related to the graph's spectral gap) and a monotonic effect (increasing beta increases concentration). It also provides intuition for practitioners: "we let information diffuse through the knowledge graph for time beta, and see how much reaches the new trial from each historical source."

### 4.4 Expert opinion aggregation and forecast combination

The literature on combining expert forecasts (e.g., Clemen 1989, "Combining forecasts: A review and annotated bibliography," International Journal of Forecasting) addresses the same abstract problem: how to weight multiple information sources when combining them. Key principles from this literature that influenced KG-DAP:

- **Diversity matters.** Simply averaging all experts (uniform weights) is often surprisingly competitive because it benefits from error cancellation. This motivated including the equal-weight mixture as a comparator — if KG-DAP cannot beat uniform weighting, the knowledge graph is not providing useful information.

- **Extreme weights are dangerous.** Putting too much weight on any single source is risky because the weight estimate may itself be wrong. This motivated the sample-size cap n_cap (preventing a single large trial from dominating the mixture) and the robustness weight w_0 (ensuring minimum mass on the vague component).

- **Source similarity is a useful proxy for source reliability.** When direct reliability estimates are unavailable, the correlation structure among sources can be used to infer which sources are informative. This is essentially what the knowledge graph provides: a structural proxy for which historical trials are likely to be informative for the current trial.

---

## 5. Design Trade-Offs

### 5.1 Pre-specifiability vs. conflict adaptation

**Gained:** Full pre-specifiability. The entire prior can be locked before data collection begins, satisfying regulatory requirements and enabling transparent communication with clinical teams.

**Sacrificed:** The ability to detect and respond to prior-data conflict at the individual trial level. If the knowledge graph incorrectly identifies a historical trial as similar, KG-DAP will borrow from it at the pre-specified weight regardless of how different its outcome is from the current trial's data. The vague component provides only global protection (not targeted).

**Mitigation:** The vague component w_0 = 0.20 bounds the maximum influence of any misleading historical information. Additionally, the posterior mixture weights automatically downweight components whose predictive distributions are inconsistent with the observed data — this is a natural Bayesian "data filtering" effect that operates within the fixed-weight framework.

### 5.2 Closed-form computation vs. model richness

**Gained:** All quantities (prior mean, variance, quantiles, CDF, posterior update, predictive distribution, ESS) are available in closed form as weighted sums of Beta-function evaluations. No MCMC, no convergence diagnostics, no chain mixing concerns. Full LOO cross-validation over 35 arms completes in seconds.

**Sacrificed:** The ability to model more complex structures. A hierarchical model could estimate between-trial heterogeneity tau from the data; a latent class model could discover clusters of similar trials; a Gaussian process on the graph could provide continuous uncertainty quantification over the graph surface. All of these richer models require MCMC or variational inference, which is slower, less reproducible, and harder to pre-specify.

**Justification:** The computational efficiency is not merely a convenience — it enables the sensitivity analyses that are essential for credible application. A method that takes 10 minutes per run cannot be explored across a grid of hyperparameters; a method that takes 2 seconds can be explored across a 5-dimensional grid with thousands of configurations.

### 5.3 Summary statistics vs. patient-level data

**Gained:** KG-DAP requires only (n_k, y_k) summaries from each historical trial, which are almost always available from published reports. This maximizes the pool of usable historical data.

**Sacrificed:** The ability to adjust for covariate imbalance between trials. When patient-level data is available, methods like LEAP (Alt et al. 2024) or propensity score-based approaches can adjust for differences in patient characteristics that the aggregate summaries miss. KG-DAP cannot exploit patient-level covariates, even when they are available.

**Justification:** In our application domain (newly diagnosed multiple myeloma MRD negativity trials), patient-level data is available for only a small fraction of the historical trials. The knowledge graph's population similarity component S_pop partially compensates by encoding population-level differences (transplant eligibility, risk profile, MRD threshold), but this is a coarser adjustment than what patient-level methods can achieve.

### 5.4 Graph quality dependence vs. autonomous learning

**Gained:** The knowledge graph provides a principled, interpretable source of structural information that can be curated by domain experts and audited for correctness.

**Sacrificed:** The ability to learn or correct the graph from outcome data. If the knowledge graph contains errors or omits relevant relationships, KG-DAP has no mechanism to detect this. The weights will reflect the graph as given, not the graph as it should be.

**Justification:** Learning graph structure from outcomes would require either (i) a much larger dataset than 35 trial arms, or (ii) strong parametric assumptions about how graph structure relates to outcomes. Both are problematic in our setting. The curated knowledge graph, while imperfect, encodes genuine biomedical knowledge (drug-target relationships, patient populations) that has been verified by domain experts.

### 5.5 Continuous weights vs. binary exchangeability

**Gained:** The diffusion kernel produces continuous, graded weights for every historical trial, reflecting the continuous nature of biomedical similarity. Two trials that are 80% similar and 70% similar receive different weights — not the same binary "exchangeable" or "not exchangeable" label.

**Sacrificed:** The discrete exchangeability framework has a cleaner theoretical interpretation: either two trials are drawn from the same population (exchangeable) or they are not. The continuous weighting of KG-DAP blurs this distinction, making it harder to state formal exchangeability assumptions.

**Justification:** Binary exchangeability is an idealization that rarely holds in practice. The MEM (Mutual Exchangeability Model; Kaizer et al. 2018) approach, which enumerates all possible binary exchangeability partitions, requires 2^K evaluations and produces posterior probabilities of exchangeability that are themselves uncertain. The continuous weighting of KG-DAP is a pragmatic alternative that avoids the combinatorial explosion and produces a single, deterministic weight for each trial.

---

## 6. Relationship to Existing Methods

### 6.1 Robust Meta-Analytic-Predictive Prior (rMAP)

**Reference:** Schmidli et al. (2014), *Biometrics* 70:1023-1032.

**Relationship:** KG-DAP and rMAP share the architecture of "informative component + vague robustness component," but differ in how the informative component is constructed. rMAP derives its informative component from a hierarchical normal-normal model fitted to all historical trials (treating them as exchangeable with random heterogeneity). KG-DAP derives its informative component from graph-weighted individual trial posteriors (treating them as differentially relevant based on knowledge graph structure).

**Key difference:** rMAP estimates the between-trial heterogeneity tau from the historical data, making the MAP predictive distribution data-dependent. KG-DAP's weights are graph-derived and fully pre-specified. In a sense, KG-DAP replaces the statistical estimation of "how different are the historical trials from each other?" (which tau answers) with the structural question "how similar is each historical trial to the new trial?" (which the knowledge graph answers). This substitution trades statistical adaptivity for structural information.

**When KG-DAP should outperform rMAP:** When the historical trials are heterogeneous but the knowledge graph correctly identifies which subset is most relevant to the new trial. In this case, rMAP's exchangeability assumption causes it to borrow equally from all historical trials (modulated only by the global tau), while KG-DAP concentrates borrowing on the most relevant subset.

**When rMAP should outperform KG-DAP:** When the knowledge graph provides poor discrimination (all trials appear equally similar) or when the graph structure is misleading. In this case, KG-DAP's weights are essentially uniform (no worse than rMAP) or actively harmful (borrowing from the wrong trials), while rMAP's data-driven tau can at least widen the predictive distribution to account for the heterogeneity.

### 6.2 LEAP (Latent Exchangeability with Aggregation of Priors)

**Reference:** Alt et al. (2024), *Biometrics*.

**Relationship:** LEAP is conceptually the closest existing method to KG-DAP. Both construct priors as weighted combinations of historical trial contributions. However, LEAP determines its exchangeability weights from the outcome data (using latent exchangeability indicators), while KG-DAP determines its weights from the knowledge graph. LEAP also operates at the patient level when individual data is available, while KG-DAP works with summary statistics.

**What KG-DAP borrows from LEAP:** The idea that borrowing should be trial-specific rather than global — that different historical trials should contribute different amounts to the prior. KG-DAP's original project name referenced LEAP explicitly (KG-LEAP).

**What KG-DAP changes:** The source of the borrowing weights. LEAP's data-driven weights are more responsive to actual outcome similarity but less pre-specifiable. KG-DAP's graph-driven weights are less responsive but fully pre-specifiable.

### 6.3 MEM (Mutual Exchangeability Model) and iMEM/dMEM

**References:** Kaizer et al. (2018); Kaizer et al. (2021, iMEM); Ji et al. (2022, dMEM).

**Relationship:** MEM enumerates all 2^K possible binary exchangeability partitions of the historical trials and computes the posterior probability of each partition given the data. This produces posterior-weighted mixtures of priors, where each partition defines a different prior for the new trial. KG-DAP can be seen as replacing MEM's exponential enumeration with a single, deterministic weight vector derived from the knowledge graph.

**Key difference:** MEM's weights are posterior probabilities of exchangeability — they adapt to the data. KG-DAP's weights are structural affinities — they are fixed before data collection. MEM's computational cost is O(2^K), which becomes prohibitive for K > 20; KG-DAP's cost is O(K^3) for the eigendecomposition (and O(K) per evaluation thereafter), which is tractable for any K.

### 6.4 Commensurate Priors

**Reference:** Hobbs et al. (2011).

**Relationship:** The commensurate prior places a spike-and-slab or continuous shrinkage prior on the difference between the current trial's parameter and each historical trial's parameter, with the degree of shrinkage adapted from the data. KG-DAP achieves a similar effect (differential borrowing from different historical trials) but uses graph structure rather than data-driven shrinkage.

**Key difference:** The commensurate prior's shrinkage parameter tau_k controls how much to borrow from trial k, and is estimated from the data. KG-DAP's weight omega_k serves the same role but is derived from the knowledge graph. The commensurate prior is more flexible (it can fully disconnect from a misleading trial) but less pre-specifiable.

### 6.5 Power Priors

**References:** Ibrahim and Chen (2000); Wang et al. (2020, propensity-score-integrated power prior).

**Relationship:** The power prior raises each historical trial's likelihood to a power a_0k in [0, 1], effectively discounting the historical data. KG-DAP's sample-size cap n_cap achieves a similar effect: when n_k > n_cap, the effective contribution of trial k is scaled down by a factor min(1, n_cap/n_k), which is equivalent to a power prior with a_0k = min(1, n_cap/n_k). The difference is that KG-DAP's discount is applied uniformly based on sample size (to prevent any single trial from dominating), while the power prior's a_0k is typically adapted based on compatibility between historical and current data.

### 6.6 Summary: positioning KG-DAP

KG-DAP occupies a specific niche in the landscape: it is the method that uses external structural information (the knowledge graph) to derive trial-specific borrowing weights, while maintaining full pre-specifiability and closed-form computation. No existing method occupies this exact position:

| Method | Trial-specific weights? | Pre-specifiable? | Uses external structure? | Closed-form? |
|--------|------------------------|-----------------|------------------------|-------------|
| rMAP | No (global tau) | Partially | No | No (MCMC) |
| LEAP | Yes | No | No | No (MCMC) |
| MEM | Yes (posterior) | No | No | No (2^K enumeration) |
| Commensurate | Yes | No | No | No (MCMC) |
| Power prior | No (global a_0) | If a_0 fixed | No | Depends |
| **KG-DAP** | **Yes** | **Yes** | **Yes (KG)** | **Yes** |

---

## 7. Why These Default Hyperparameters

### Composite similarity weights: alpha = 0.20, beta_w = 0.20, gamma = 0.60

**Population similarity (gamma = 0.60) dominates.** In the NDMM MRD negativity setting, response rates are most directly driven by patient characteristics: transplant eligibility (TE patients consistently have higher MRD negativity rates than transplant-ineligible patients), cytogenetic risk profile (high-risk patients have lower response rates), disease stage (ISS stage 3 is prognostic), and MRD sensitivity threshold (10^{-5} vs. 10^{-4}). These patient-level factors explain more variance in response rates across trials than the drug identity alone, because the same drug used in different populations produces different rates. Therefore, population similarity receives the highest weight.

**Drug and target similarity share the remaining weight equally (alpha = beta_w = 0.20).** Drug identity and molecular target identity capture overlapping but distinct information. Two trials using the same drug automatically share the same molecular targets, but two trials using different drugs may still share molecular targets (e.g., carfilzomib and bortezomib are both proteasome inhibitors targeting the 20S proteasome subunit). Splitting the remaining 0.40 equally between drugs (direct mechanistic overlap) and targets (indirect mechanistic overlap via 2-hop paths) avoids privileging one level of abstraction over the other.

**Why not learn the weights?** The composite similarity weights could, in principle, be estimated via cross-validation (e.g., minimize LOO MAE over the weight simplex). This was considered but rejected for the default specification because: (i) with 35 trial arms, cross-validation estimates of optimal weights are noisy — different random seeds in the data split would give different optima; (ii) learned weights make the method less pre-specifiable (the weights depend on the data used for cross-validation); (iii) the diffusion kernel smooths over fine-grained differences in the similarity matrix, so moderate changes in the composite weights have a second-order effect on the final borrowing weights. Sensitivity analysis over the weight simplex is recommended, but defaults should be defensible on domain grounds without data-driven optimization.

### Diffusion parameter: beta = 2.0

**Scale calibration.** The normalized graph Laplacian L has eigenvalues in [0, 2]. For beta = 2.0, the diffusion kernel attenuates the highest-frequency eigenmode (eigenvalue 2.0) by a factor of exp(-2.0 * 2.0) = exp(-4) approximately equal to 0.018 — nearly complete suppression of local fluctuations. The lowest non-trivial eigenmode (which captures the graph's dominant partition) is attenuated by exp(-2.0 * lambda_2), where lambda_2 is the spectral gap. For a well-connected graph, lambda_2 is typically between 0.1 and 0.5, giving an attenuation of exp(-0.2) to exp(-1.0), i.e., 37% to 82% preservation. This places beta = 2.0 in the "moderate concentration" regime — neither uniform borrowing (beta -> 0) nor nearest-neighbor borrowing (beta -> infinity).

**Empirical rationale.** In preliminary experiments with the NDMM dataset, beta = 2.0 produced weight vectors that are concentrated on the 3-5 most similar trials for each left-out arm, with the remaining trials receiving small but nonzero weights. This matches the clinical intuition that for any given trial, there are a handful of "most relevant" historical comparators plus a long tail of "somewhat relevant" trials.

**Sensitivity.** Beta is the most important hyperparameter for sensitivity analysis. The specification recommends exploring beta in {0.5, 1.0, 2.0, 5.0, 10.0} as a standard sensitivity grid.

### Robustness weight: w_0 = 0.20

**Protection floor.** With w_0 = 0.20, the vague Beta(1,1) component always retains at least 20% of the prior mass. This means that even in the worst case (all historical trials are completely misleading), the prior is no more informative than a prior with ESS of approximately 0.80 * max_k(ESS_k) + 0.20 * 2. In practice, this bounds the maximum information contribution and prevents the prior from being so concentrated that it overwhelms reasonable amounts of current trial data.

**Convention alignment.** The rMAP prior (Schmidli et al. 2014) uses a robustness weight of w_rob = 0.20 as its default. Using the same value for KG-DAP ensures a fair comparison: any difference in performance is due to the informative component, not the robustness tuning.

**Not too conservative.** A larger w_0 (e.g., 0.50) would make KG-DAP closer to the non-informative prior, diluting the value of the knowledge graph. A smaller w_0 (e.g., 0.05) would make the prior more informative but more vulnerable to graph errors. The value 0.20 represents a standard compromise in the Bayesian historical borrowing literature.

### Base prior: alpha_0 = beta_0 = 1.0 (uniform)

**Minimally informative.** Using Beta(1,1) = Uniform(0,1) as the base prior for each historical arm means that pi_k(theta) = Beta(1 + y_k, 1 + n_k - y_k) is simply the posterior under maximum ignorance. This ensures that each historical arm's contribution is driven entirely by its own data, not by a subjective base prior that might introduce additional bias.

**Jeffreys alternative.** Using Beta(0.5, 0.5) (the Jeffreys prior) instead of Beta(1,1) would produce slightly different posteriors for small n_k but has negligible effect for the sample sizes in our data (n_k ranging from approximately 20 to 700). The uniform prior was chosen for simplicity and transparency.

### Sample-size cap: n_cap = 200

**Preventing single-trial dominance.** The historical trial arms range from n approximately 20 to n approximately 700 in the NDMM dataset. Without a cap, a trial with n = 700 would contribute a Beta component with ESS approximately 700, which would dominate the mixture regardless of its weight omega_k. The cap n_cap = 200 ensures that no single trial contributes more than ESS approximately 200 to the mixture, so the final prior reflects the graph-derived weighting structure rather than raw data volume.

**Scale.** The cap of 200 was chosen to be large enough that typical Phase 2 trials (n = 50-150) are unaffected, while Phase 3 trials (n = 300-700) are modestly discounted. It is equivalent to saying: "we trust each historical trial's data, but we do not let any single trial dominate the prior just because it happened to enroll more patients." This aligns with the principle that the knowledge graph, not the sample size, should determine the relative importance of different historical trials.

---

## 8. Literature Basis

The following papers contributed specific ideas to KG-DAP:

| Paper | Contribution to KG-DAP |
|-------|----------------------|
| Schmidli et al. (2014), *Biometrics* | Architecture of informative + vague component; the rMAP framework that KG-DAP extends |
| Kondor and Lafferty (2002), ICML | Diffusion kernel on graphs: K_beta = exp(-beta * L) |
| Smola and Kondor (2003), COLT | Regularization-theoretic justification of graph diffusion kernels |
| Alt et al. (2024), *Biometrics* (LEAP) | Trial-specific borrowing weights via latent exchangeability; the conceptual predecessor to KG-DAP's graph-derived weights |
| Kaizer et al. (2018) (MEM) | Exchangeability-based mixture priors for historical borrowing |
| Hobbs et al. (2011) | Commensurate priors: the idea of borrowing proportional to compatibility |
| Ibrahim and Chen (2000) | Power priors: discounting historical data by a scalar weight |
| Morita et al. (2008) | ESS calculation for informative priors via moment-matching |
| Diaconis and Ylvisaker (1979) | Conjugacy of finite mixtures of exponential family distributions |
| Clemen (1989) | Forecast combination principles: diversity, extreme-weight avoidance |

---

## 9. Stress Test Results

### Strongest argument against KG-DAP

"The knowledge graph is a curation artifact, not a natural object. Different curators would construct different graphs with different edges, producing different similarity matrices and different borrowing weights. The method's output is therefore sensitive to subjective curation choices that are not part of the formal statistical model."

**Response:** This is a genuine concern. However, the same argument applies to any method that uses covariates or structural information — the choice of which covariates to include, how to encode them, and how to weight them is always partly subjective. KG-DAP makes this subjectivity transparent (the graph is an explicit, inspectable object) and controllable (the composite similarity weights alpha, beta_w, gamma can be varied in sensitivity analyses). The alternative — ignoring structural information entirely (rMAP) or learning it from small-sample outcomes (commensurate prior) — has its own vulnerabilities: rMAP wastes available information, and data-learned structure is noisy with K = 35.

### Second-strongest argument

"With only 35 trial arms, the knowledge graph's 56 nodes and their connections may encode spurious patterns. Two trials may appear similar in the graph simply because multiple myeloma is a single disease with a limited pharmacopoeia — most trials use combinations drawn from the same 11 drugs. High similarity may reflect lack of diversity rather than genuine relevance."

**Response:** This is mitigated by the population similarity component (gamma = 0.60), which captures trial-level differences in eligibility criteria, risk profile, and endpoint sensitivity — factors that vary meaningfully across trials even when the drug combinations overlap. The LOO cross-validation is specifically designed to test whether graph-derived similarity actually predicts outcome similarity: if it does not, KG-DAP will not outperform the equal-weight mixture comparator.

---

## 10. Open Questions

These questions can only be answered through implementation and evaluation:

1. **Does graph diffusion actually outperform raw similarity?** It is possible that the higher-order propagation in the diffusion kernel adds noise rather than signal for this particular graph. The equal-weight mixture comparator will help answer this.

2. **Is beta = 2.0 near the optimal value?** The sensitivity analysis will reveal the MAE-minimizing beta and whether the performance surface has a sharp or flat optimum.

3. **Is the ESS bound tight?** The conjectured ESS bound (Property 2) depends on the mixture ESS being bounded by the maximum component ESS. For well-separated mixture components, the actual ESS can be much smaller than this bound. Empirical ESS values will reveal how tight the bound is in practice.

4. **Does the robustness component activate appropriately?** In the posterior update, the vague component's weight should increase when the observed data conflicts with all historical components. Whether this "automatic robustness" is sufficient in practice (compared to explicit conflict-detection mechanisms) is an empirical question.

5. **Can the method handle prospective validation?** The LOO cross-validation is retrospective. Prospective validation — using the method to predict a newly published trial's results before they are observed — would provide stronger evidence of practical value. The Elo-KRd patient-level dataset is a candidate for this test.

6. **Sensitivity to graph curation choices.** How much do the results change if edges are added or removed from the knowledge graph? A perturbation analysis (randomly dropping 10-20% of edges and re-running LOO) would quantify this.
