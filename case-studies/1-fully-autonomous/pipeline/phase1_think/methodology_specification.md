# Methodology Specification: KG-DAP (Knowledge-Graph-Driven Adaptive Prior)

## One-Sentence Summary

KG-DAP constructs an informative Bayesian prior for a new clinical trial's response rate by forming a weighted mixture of historical trial posteriors, where the mixing weights are derived from a biomedical knowledge graph via a graph diffusion kernel, with a vague robustness component that automatically limits borrowing when the graph structure is uninformative or misleading.

---

## Core Idea

Given a current (new) trial arm c and a set of K historical trial arms {h_1, ..., h_K}, each with summary data (n_k, y_k), the KG-DAP prior for the current trial's response rate theta_c is:

    pi_{DAP}(theta_c) = (1 - w_0) * sum_{k=1}^{K} omega_k * pi_k(theta_c) + w_0 * pi_vague(theta_c)

where:

- **omega_k** are knowledge-graph-derived weights satisfying omega_k >= 0 and sum_k omega_k = 1, computed from a diffusion kernel on the KG
- **pi_k(theta_c)** is the contribution from historical trial k, derived as a Beta posterior given its summary data and a minimally informative base prior
- **w_0** is the robustness weight for the vague component pi_vague
- **pi_vague(theta_c) = Beta(1, 1)** is a uniform (maximally uninformative) prior

The key innovation is that the weights omega_k are derived by applying a **graph diffusion kernel** to the knowledge graph, which captures not only direct structural similarity (shared drugs) but also higher-order similarity (shared molecular targets via different drugs, similar patient populations connected through treatment pathways). The diffusion parameter beta controls the "radius" of information propagation on the graph, providing a single, interpretable hyperparameter that governs the sharpness of borrowing.

---

## Mathematical Specification

### Notation

| Symbol | Definition | Domain |
|--------|-----------|--------|
| theta_c | Response rate (MRD negativity rate) for the current trial arm | [0, 1] |
| (n_k, y_k) | Sample size and number of responders for historical arm k | n_k in N, y_k in {0,...,n_k} |
| K | Number of historical trial arms | positive integer |
| G = (V, E) | Biomedical knowledge graph | MultiDiGraph |
| S | K x K composite similarity matrix derived from G | S_ij in [0, 1] |
| beta | Power sharpness parameter | (0, infty) |
| c | Index identifying the current trial in the graph | integer |
| omega | Normalized borrowing weight vector | K-simplex |
| w_0 | Robustness weight for the vague component | [0, 1] |
| delta | Similarity threshold for weight concentration | [0, 1] |
| alpha_0, beta_0 | Base prior parameters for each historical arm | positive reals |
| ESS | Effective sample size of the KG-DAP prior | positive real |

### Assumptions

1. **Conditional independence**: Given theta_c, the new trial data y_c | theta_c ~ Binomial(n_c, theta_c)
2. **Historical summary sufficiency**: Each historical arm k is summarized by (n_k, y_k), and we treat the posterior Beta(alpha_0 + y_k, beta_0 + n_k - y_k) as a sufficient representation
3. **Knowledge graph faithfulness**: The KG structure reflects genuine biomedical relationships that are relevant to treatment efficacy (drugs target specific proteins, trials enroll specific populations)
4. **Positive definiteness of similarity**: The composite similarity matrix S has entries in [0,1] with S_ii = 1

### Step 1: Composite Similarity Matrix

The composite similarity between any two trial arms i and j is:

    S_ij = alpha * J_drug(i,j) + beta_w * J_target(i,j) + gamma * S_pop(i,j)

where:
- J_drug(i,j) = |D_i intersect D_j| / |D_i union D_j| is the Jaccard similarity of drug sets
- J_target(i,j) = |T_i intersect T_j| / |T_i union T_j| is the Jaccard similarity of 2-hop molecular target sets
- S_pop(i,j) = 1 - ||x_i - x_j||_1 / K_pop is the normalized L1 population similarity (K_pop = number of population features)
- alpha + beta_w + gamma = 1 are component weights (default: 0.20, 0.20, 0.60)

**Justification for default weights**: Population similarity (gamma = 0.60) is weighted highest because patient characteristics (transplant eligibility, risk profile, MRD threshold) have the strongest direct effect on response rates. Drug and target similarities (alpha = beta_w = 0.20 each) capture mechanistic overlap at two levels of abstraction.

### Step 2: Power-Based Similarity Weights

Given the composite similarity vector s_c = (S_{c,1}, ..., S_{c,K}) from the current trial c to each historical arm, compute the borrowing weights via power sharpening:

    omega_k = S_{c,k}^beta / sum_{j=1}^{K} S_{c,j}^beta

where beta > 0 is the sharpness parameter controlling the concentration of borrowing.

**Why power-based weights instead of the diffusion kernel**: The originally proposed diffusion kernel exp(-beta * L) requires constructing a graph Laplacian. On dense similarity graphs (where most pairwise similarities are in [0.3, 0.9], as is typical for trials within a single therapeutic area), the normalized Laplacian's eigenvalues compress into a narrow band, making the kernel unable to discriminate between similar and dissimilar trials regardless of beta. Power-based weights directly amplify similarity differences without the Laplacian bottleneck. This was discovered and validated in Phase 2 (see validation_report.md).

**Properties of the power-based weights**:
- As beta -> 0: omega_k -> 1/K (uniform borrowing, ignoring graph structure)
- As beta -> infty: omega concentrates on the trial with highest S_{c,k} (nearest-neighbor borrowing)
- For moderate beta (validated range [8, 30]): weights smoothly concentrate on structurally similar trials while suppressing dissimilar ones
- **Monotonicity**: If S_{c,j} > S_{c,k}, then omega_j > omega_k for all beta > 0 (directly verifiable)
- **Pre-specifiability preserved**: weights depend only on the KG structure, not on outcome data

### Step 4: Historical Arm Contributions

For each historical arm k, form the posterior distribution under a minimally informative base prior Beta(alpha_0, beta_0):

    pi_k(theta) = Beta(theta | alpha_0 + y_k, beta_0 + n_k - y_k)

Default: alpha_0 = beta_0 = 1 (uniform base prior), so pi_k = Beta(1 + y_k, 1 + n_k - y_k).

**Sample-size scaling** (optional refinement): When historical trials have very large sample sizes, their posteriors are very concentrated. To prevent a single large trial from dominating, we can cap the effective contribution of each historical arm:

    pi_k^{scaled}(theta) = Beta(theta | alpha_0 + y_k * min(1, n_cap / n_k), beta_0 + (n_k - y_k) * min(1, n_cap / n_k))

where n_cap is a cap on the effective sample size (default: n_cap = 200). This is equivalent to using a power prior with a_0k = min(1, n_cap/n_k) for each arm individually.

### Step 5: Robustness Component

The vague component provides protection against misleading KG structure:

    pi_vague(theta) = Beta(theta | 1, 1)    [uniform on (0,1)]

The robustness weight w_0 controls the minimum probability mass on the vague component. Default: w_0 = 0.20.

**Interpretation**: Even if the KG structure is entirely misleading, at least 20% of the prior mass comes from an uninformative source, bounding the potential for harm.

### Step 6: The KG-DAP Prior

The complete prior is the finite mixture:

    pi_{DAP}(theta_c) = (1 - w_0) * sum_{k=1}^{K} omega_k * Beta(theta_c | a_k', b_k') + w_0 * Beta(theta_c | 1, 1)

where a_k' = alpha_0 + y_k (or the scaled version) and b_k' = beta_0 + n_k - y_k.

This is a **(K+1)-component Beta mixture**, which is:
- Analytically tractable (moments, quantiles, and CDF computed in closed form)
- Conjugate-like: the posterior after observing (n_c, y_c) from the current trial is also a Beta mixture with updated parameters
- Interpretable: each component corresponds to a specific historical trial or the vague fallback

### Step 7: Posterior After Observing Current Trial Data

Given current trial data y_c | theta_c ~ Binomial(n_c, theta_c), the posterior is:

    pi(theta_c | y_c) propto L(y_c | theta_c) * pi_{DAP}(theta_c)
                      = sum_{k=0}^{K} v_k * Beta(theta_c | a_k' + y_c, b_k' + n_c - y_c)

where the updated mixture weights are:

    v_k propto w_k^{prior} * B(a_k' + y_c, b_k' + n_c - y_c) / B(a_k', b_k')

with w_k^{prior} = (1-w_0)*omega_k for k=1,...,K and w_0^{prior} = w_0, and B(.,.) is the Beta function. This is the standard Bayesian update for a finite mixture of conjugate priors.

### Step 8: Effective Sample Size (ESS)

The prior ESS quantifies how much information the KG-DAP prior contributes, measured as the number of hypothetical observations from a Beta distribution that would produce the same information content (Morita et al. 2008):

    ESS = a_prior + b_prior

where a_prior and b_prior are the parameters of the "closest" single Beta distribution to the mixture, found by moment-matching:

    mu_prior = E[theta_c]     (mixture mean)
    sigma^2_prior = Var[theta_c]   (mixture variance)
    a_prior = mu_prior * (mu_prior * (1 - mu_prior) / sigma^2_prior - 1)
    b_prior = (1 - mu_prior) * (mu_prior * (1 - mu_prior) / sigma^2_prior - 1)

### Step 9: Predictive Distribution

The prior predictive distribution for a new observation y_new ~ Binomial(n_new, theta_c) under the KG-DAP prior is:

    P(y_new = y | n_new) = sum_{k=0}^{K} w_k * BetaBin(y | n_new, a_k', b_k')

where BetaBin is the Beta-Binomial probability mass function. This enables pre-trial design operating characteristics.

---

## Algorithm

```
ALGORITHM: KG-DAP Prior Construction
======================================

Input:
  - Knowledge graph G = (V, E)
  - Historical trial data {(n_k, y_k)}_{k=1}^K
  - Current trial identifier c (node in G)
  - Parameters: alpha, beta_w, gamma (similarity weights)
                beta (diffusion parameter)
                w_0 (robustness weight)
                alpha_0, beta_0 (base prior)
                n_cap (sample size cap, optional)

Output:
  - KG-DAP prior: mixture weights {w_k}_{k=0}^K and Beta parameters {(a_k', b_k')}_{k=0}^K

Step 1: Compute composite similarity matrix
  For each pair (i, j) of trial arms (including c):
    S_ij = alpha * J_drug(i,j) + beta_w * J_target(i,j) + gamma * S_pop(i,j)

Step 2: Compute power-based borrowing weights
  For k = 1, ..., K:
    raw_k = S[c, k]^beta
  omega = raw / sum(raw)    [normalized to sum to 1]

Step 4: Form historical arm posteriors
  For k = 1, ..., K:
    scale_k = min(1, n_cap / n_k)    [skip if n_cap = infinity]
    a_k' = alpha_0 + y_k * scale_k
    b_k' = beta_0 + (n_k - y_k) * scale_k

Step 5: Assemble the KG-DAP prior
  Component 0 (vague): weight = w_0, parameters = (1, 1)
  Component k (k=1..K): weight = (1 - w_0) * omega_k, parameters = (a_k', b_k')

Step 6: (Optional) Compute ESS
  mu = sum_k w_k * a_k' / (a_k' + b_k')
  mu2 = sum_k w_k * [a_k'*(a_k'+1)] / [(a_k'+b_k')*(a_k'+b_k'+1)]
  var = mu2 - mu^2
  ESS = mu*(1-mu)/var - 1

Return {(w_k, a_k', b_k')}_{k=0}^K, ESS
```

```
ALGORITHM: KG-DAP Posterior Update
====================================

Input:
  - KG-DAP prior: {(w_k, a_k', b_k')}_{k=0}^K
  - Current trial data: (n_c, y_c)

Output:
  - Posterior mixture: {(v_k, a_k'', b_k'')}_{k=0}^K

For k = 0, ..., K:
  a_k'' = a_k' + y_c
  b_k'' = b_k' + n_c - y_c
  log_marginal_k = lbeta(a_k'', b_k'') - lbeta(a_k', b_k')
  log_v_k = log(w_k) + log_marginal_k

Normalize: v_k = exp(log_v_k) / sum_j exp(log_v_j)

Return {(v_k, a_k'', b_k'')}_{k=0}^K
```

---

## Parameters

| Parameter | Symbol | Default | Valid Range | Interpretation |
|-----------|--------|---------|-------------|----------------|
| Drug similarity weight | alpha | 0.20 | [0, 1] | Weight of drug-level Jaccard in composite similarity |
| Target similarity weight | beta_w | 0.20 | [0, 1] | Weight of molecular target Jaccard in composite similarity |
| Population similarity weight | gamma | 0.60 | [0, 1] | Weight of population similarity in composite similarity |
| Power sharpness parameter | beta | 15.0 | (0, infty) | Controls sharpness of graph-based weight concentration; higher = more concentrated on similar trials. Validated safe range: [8, 30] |
| Robustness weight | w_0 | 0.20 | [0, 1) | Minimum weight on vague (uniform) component |
| Base prior alpha | alpha_0 | 1.0 | (0, infty) | First parameter of base prior for each historical arm |
| Base prior beta | beta_0 | 1.0 | (0, infty) | Second parameter of base prior for each historical arm |
| Sample size cap | n_cap | 200 | (0, infty] | Maximum effective sample size contributed by any single historical arm; infty = no cap |

**Constraint**: alpha + beta_w + gamma = 1.

**Hyperparameter sensitivity**: The primary hyperparameters requiring sensitivity analysis are beta (diffusion sharpness) and w_0 (robustness weight). The composite similarity weights (alpha, beta_w, gamma) can be explored but have a secondary effect given that the diffusion kernel smooths over fine-grained differences in the similarity matrix.

---

## Comparators

### 1. Robust Meta-Analytic-Predictive Prior (rMAP)

**Reference**: Schmidli et al. (2014), *Biometrics* 70:1023--1032.

**Description**: A Bayesian hierarchical normal-normal model is fit to the logit-transformed historical rates, producing a MAP predictive distribution for the new trial. A vague component is mixed in for robustness:

    pi_{rMAP}(theta) = (1 - w_rob) * pi_{MAP}(theta) + w_rob * pi_vague(theta)

**Implementation**:
1. Fit: logit(p_hat_k) ~ N(mu, tau^2), mu ~ N(0, 10^2), tau ~ HalfNormal(1)
2. Derive: pi_{MAP}(logit(theta)) = t-distribution (marginalizing mu, tau)
3. Approximate: Laplace approximation or numerical integration on the probability scale
4. Mix: add vague component with weight w_rob = 0.20

**Why this comparator**: rMAP is the most widely used informative prior method in pharma and has regulatory acceptance. It represents the state of the art for "graph-free" borrowing.

### 2. Uniform (Non-Informative) Prior

**Description**: pi(theta) = Beta(1, 1). This is the baseline representing no borrowing at all.

**Why this comparator**: Establishes the value of any borrowing. If KG-DAP cannot beat this, borrowing is not worthwhile.

### 3. Pooled Beta Prior (Full Borrowing)

**Description**: Pool all historical data into a single Beta posterior:

    pi_{pool}(theta) = Beta(1 + sum_k y_k, 1 + sum_k (n_k - y_k))

**Why this comparator**: Represents maximum (naive) borrowing. Shows the cost of ignoring heterogeneity.

### 4. Equal-Weight Mixture (Unweighted KG-DAP)

**Description**: Same as KG-DAP but with omega_k = 1/K for all k (ignoring the knowledge graph).

**Why this comparator**: Isolates the contribution of the KG-derived weights vs. the mixture architecture itself.

---

## Minimum Viable Evaluation

### Experiment: Leave-One-Out Cross-Validation on 35 NDMM Trial Arms

**Protocol**: For each arm h in {1, ..., 35}:
1. Hold out arm h's data (n_h, y_h)
2. Construct KG-DAP prior from the remaining 34 arms using the knowledge graph
3. Compute the prior predictive distribution for arm h
4. Compute metrics comparing the prior prediction to the observed (n_h, y_h)

**Primary Metric**: Mean Absolute Error (MAE) of the prior predictive mean vs. observed rate:

    MAE = (1/35) * sum_h |E_{prior_h}[theta_h] - y_h/n_h|

**Secondary Metrics**:
- **Coverage**: Proportion of arms where the observed rate falls within the 95% prior predictive interval
- **Interval Width**: Mean width of the 95% prior predictive interval (measures efficiency)
- **Log Predictive Score**: (1/35) * sum_h log P(y_h | n_h, prior_h) (proper scoring rule)

**Success Threshold**:
- MAE(KG-DAP) < MAE(rMAP) by at least 5% relative improvement
- Coverage(KG-DAP) in [0.88, 0.98] (calibrated, not over- or under-covering)
- MAE(KG-DAP) < MAE(Uniform) (borrowing provides value)

---

## Success Criterion

| Criterion | Metric | Minimum Required | Comparison Target | Notes |
|-----------|--------|------------------|-------------------|-------|
| **Primary** | LOO MAE reduction | >= 5% relative to rMAP | rMAP (Schmidli 2014) | Averaged over all 35 arms |
| **Calibration** | 95% interval coverage | 0.88 -- 0.98 | rMAP, Uniform | Both prior and posterior |
| **Robustness** | MAE degradation under contamination | < 5% absolute increase | KG-DAP without contamination | Add 10 adjacent + 10 foreign trials to pool |
| **Robustness** | Coverage under contamination | > 0.85 | Maintain calibration | After adding irrelevant trials |
| **Efficiency** | Mean interval width | Narrower than rMAP | rMAP | Especially for small-n arms (n < 100) |
| **Computation** | Wall-clock time for full LOO | < 10 seconds | -- | On standard laptop |

---

## Theoretical Properties (Expected)

### Property 1: Conjugate Closure (PROVEN)

**Statement**: If the KG-DAP prior is a finite Beta mixture and the likelihood is Binomial, the posterior is also a finite Beta mixture (with re-weighted components).

**Status**: PROVEN (standard result for conjugate mixtures; see Diaconis & Ylvisaker 1979).

**Proof sketch**: Each component Beta(a_k, b_k) is conjugate to Binomial. The posterior mixture weight for component k is proportional to w_k * B(a_k + y, b_k + n - y) / B(a_k, b_k), where B is the Beta function. The posterior is therefore a (K+1)-component Beta mixture with updated parameters and re-normalized weights.

### Property 2: Bounded ESS (CONJECTURED)

**Statement**: The effective sample size of the KG-DAP prior satisfies:

    ESS(pi_{DAP}) <= (1 - w_0) * max_k(n_k * min(1, n_cap/n_k) + alpha_0 + beta_0) + 2

**Status**: CONJECTURED. The bound follows from the fact that the mixture ESS is bounded by the maximum component ESS (the most informative single component), not their sum. The "+2" comes from the vague Beta(1,1) component's ESS.

**Significance**: This guarantees that the KG-DAP prior cannot contribute more effective information than the single most informative historical arm (scaled), even when many historical arms are included. This is a crucial safety property.

### Property 3: Asymptotic Dominance of Data (PROVEN)

**Statement**: As n_c -> infty, the posterior pi(theta_c | y_c) under the KG-DAP prior converges to the likelihood-driven estimate regardless of the prior:

    pi(theta_c | y_c) -> delta(y_c / n_c) as n_c -> infty

**Status**: PROVEN (Bernstein-von Mises theorem applies since the prior has full support on (0,1) via the vague component with w_0 > 0).

### Property 4: Monotonic Borrowing (PROVEN)

**Statement**: If S_{c,j} > S_{c,k}, then omega_j > omega_k for all beta > 0.

**Status**: PROVEN. With power-based weights omega_k = S_{c,k}^beta / Z, if S_{c,j} > S_{c,k} then S_{c,j}^beta > S_{c,k}^beta for all beta > 0 (since x -> x^beta is strictly monotone on (0, infty)). The normalization constant Z is the same for both, so omega_j > omega_k.

### Property 5: Power Parameter Limits (PROVEN)

**Statement**:
- As beta -> 0: omega_k -> 1/K for all k (uniform weights; no graph discrimination)
- As beta -> infty: omega concentrates on the trial with highest S_{c,k} (nearest-neighbor borrowing)

**Status**: PROVEN. As beta -> 0, S_{c,k}^beta -> 1 for all k, so omega_k -> 1/K. As beta -> infty, the term with largest S_{c,k} dominates the sum, so omega concentrates on argmax_k S_{c,k}.

---

## Expected Strengths

1. **Principled graph-based borrowing**: Unlike rMAP (which treats all historical trials as exchangeable up to random heterogeneity) or MEM (which uses binary exchangeability), KG-DAP uses the continuous graph structure to derive graded, trial-specific borrowing weights. Trials that share more drugs, molecular targets, and patient characteristics receive higher weight.

2. **Pre-specifiability**: The borrowing weights omega_k are fully determined by the knowledge graph and the diffusion parameter beta, without reference to the outcome data of the current trial. This means the prior can be fully specified in a regulatory Statistical Analysis Plan before unblinding. This is a major practical advantage over data-driven methods (commensurate prior, elastic prior, LEAP) that adapt borrowing based on observed outcomes.

3. **Computational efficiency**: The prior is a finite Beta mixture, so all key quantities (moments, quantiles, CDF, predictive distribution, posterior update) are available in closed form or as weighted sums of closed-form expressions. No MCMC is required. Full LOO cross-validation over 35 arms takes seconds, not minutes.

4. **Built-in robustness**: The vague component with weight w_0 ensures that the prior always has full support on (0,1), providing a floor of protection against misleading graph structure. Even if the knowledge graph is entirely wrong, the posterior is no worse than using a prior with ESS = w_0 * 2 (essentially uninformative).

5. **Higher-order structural similarity via diffusion**: Raw Jaccard similarity only captures direct overlap (shared drugs). The diffusion kernel propagates similarity through the graph, capturing indirect relationships: two trials that use different drugs targeting the same molecular pathway will have elevated diffusion-kernel similarity even if their drug Jaccard is zero.

6. **Handles heterogeneous sample sizes**: The optional sample-size cap n_cap prevents large Phase 3 trials (n > 500) from dominating the prior, ensuring that the mixture reflects structural similarity rather than raw data volume.

7. **Natural uncertainty quantification**: The Beta mixture prior naturally produces wider credible intervals when the historical evidence is conflicting (components are dispersed) and narrower intervals when it is concordant (components cluster together).

---

## Expected Weaknesses

1. **Hyperparameter sensitivity**: The method has 5 primary hyperparameters (alpha, beta_w, gamma, beta, w_0). While defaults are provided, the optimal values may be context-dependent. The diffusion parameter beta is particularly important: too small yields uniform borrowing (wasting the graph), too large yields nearest-neighbor borrowing (ignoring useful moderate-similarity trials). Sensitivity analysis is essential.

2. **Knowledge graph quality dependence**: KG-DAP is only as good as the underlying knowledge graph. If the graph encodes relationships that are not relevant to treatment efficacy (e.g., two drugs target the same protein but one is effective and the other is not due to pharmacokinetic differences), the graph-derived weights will be suboptimal. The method provides no mechanism to learn that the graph is wrong from outcome data (by design, for pre-specifiability).

3. **Mixture identifiability**: With K+1 components, the Beta mixture is over-parameterized relative to the information content. In practice this is handled by the fact that the mixture weights are fixed (not estimated), but the posterior mixture can have many small-weight components that contribute little. Computational overhead is O(K) per evaluation.

4. **Limited to summary statistics**: KG-DAP uses (n_k, y_k) summaries. It cannot exploit patient-level covariates from historical trials (unlike LEAP). When patient-level data is available, methods that adjust for covariate imbalance may be preferable.

5. **Assumption that graph structure predicts outcome similarity**: The fundamental assumption that structurally similar trials (in the KG sense) have similar response rates is not always true. Trials using the same drugs in different disease stages, or with different dosing schedules, may have very different outcomes despite high graph similarity.

---

## Appendix A: Composite Similarity Details

### Drug Jaccard (J_drug)
    J_drug(i,j) = |drugs(i) intersect drugs(j)| / |drugs(i) union drugs(j)|

Uses the drug nodes directly connected to each trial arm in the KG.

### Target Jaccard (J_target)
    J_target(i,j) = |targets(i) intersect targets(j)| / |targets(i) union targets(j)|

Uses 2-hop neighbors: trial -> drug -> molecular_target.

### Population Similarity (S_pop)
    S_pop(i,j) = 1 - ||x_i - x_j||_1 / K_pop

where x is a K_pop-dimensional feature vector with components:
- TE indicator: 1.0 (TE), 0.5 (mixed), 0.0 (TIE)
- frac_high_risk (proportion of high-risk cytogenetics patients)
- median_age / 100 (normalized age)
- frac_iss3 (proportion with ISS stage 3)
- MRD threshold indicator: 1.0 (10^-5), 0.0 (10^-4)
- ASCT required: 0 or 1

K_pop = 6 (number of features, also maximum possible L1 distance).

## Appendix B: Comparison with Previous KG-CAR Approach

The previous approach (KG-CAR, from the hand-crafted paper) used a conditional autoregressive (CAR) prior on the logit-transformed response rates:

    phi | sigma_phi ~ N(0, sigma_phi^2 * (D_w - rho * W)^{-1})
    logit(theta_k) = mu + phi_k + epsilon_k

Key differences from KG-DAP:

| Aspect | KG-CAR | KG-DAP |
|--------|--------|--------|
| Working scale | Logit (Gaussian) | Probability (Beta) |
| Prior structure | Joint multivariate normal | Finite Beta mixture |
| Computation | MCMC or Laplace approximation | Closed-form |
| Parameters estimated from data | mu, sigma_phi, sigma_eps, rho | None (fully pre-specified) |
| Robustness mechanism | None explicit | Vague component w_0 |
| Graph utilization | Adjacency matrix W | Diffusion kernel exp(-beta*L) |
| Handles new trial not in graph | Awkward (must extend joint prior) | Natural (compute similarity to all historical, form mixture) |

KG-DAP is designed to address two key limitations of KG-CAR:
1. **Computational complexity**: KG-CAR requires MCMC or multivariate Laplace approximation; KG-DAP is closed-form
2. **Pre-specifiability**: KG-CAR estimates the spatial correlation parameter rho from data, making the prior data-dependent; KG-DAP's weights are entirely pre-specified from the graph

## Appendix C: Extension to Randomized Controlled Trials

For a randomized trial with treatment arm T and control arm C, KG-DAP can be applied independently to each arm:

    pi_{DAP}(theta_T) = mixture based on historical treatment arms
    pi_{DAP}(theta_C) = mixture based on historical control arms

The treatment effect delta = theta_T - theta_C can then be estimated by:
1. Computing the full posterior for each arm
2. Deriving P(theta_T > theta_C | data) via numerical integration over the mixture posteriors

This factored approach is natural when the KG contains arm-level (not trial-level) nodes.
