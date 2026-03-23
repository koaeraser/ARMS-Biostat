# Literature Briefing: KG-Informed Bayesian Borrowing for Clinical Trials

**Prepared for**: Phase 1 THINK (methodology-architect)
**Date**: 2026-03-23

---

## 1. The Borrowing Landscape in Clinical Trial Statistics

### 1.1 Power Prior (Ibrahim & Chen 2000)

The power prior raises the historical likelihood to a discount parameter a0 in [0,1]:

    pi(theta | D0, a0) propto L(theta | D0)^{a0} * pi_0(theta)

- a0 = 0: ignore history; a0 = 1: full borrowing
- Extension to multiple historical datasets: separate a0k per dataset (Ibrahim et al. 2015)
- Limitation: a0 is typically fixed or given a Beta prior; it does not leverage structural knowledge about trial similarity

**Key reference**: Ibrahim MH, Chen MH, Gwon Y, Chen F. The power prior: theory and applications. *Statist. Sci.* 2015;30:46--60.

### 1.2 Robust Meta-Analytic-Predictive Prior (Schmidli et al. 2014)

Derives a MAP prior from a Bayesian hierarchical meta-analysis of historical data, then adds a vague component:

    pi_rMAP(theta) = (1 - w_rob) * pi_MAP(theta) + w_rob * pi_vague(theta)

- pi_MAP is a (possibly approximated) marginal predictive from a normal-normal hierarchical model
- w_rob (typically 0.1--0.5) provides robustness against prior-data conflict
- Limitation: ALL historical trials are pooled through a single between-trial heterogeneity parameter tau; no trial-specific borrowing weights; treats all sources as equidistant

**Key reference**: Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B. Robust meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics* 2014;70:1023--1032.

### 1.3 Commensurate Prior (Hobbs et al. 2011, 2012)

Places a spike-and-slab-like prior on the commensurability between current and historical parameters:

    theta_current | theta_hist ~ N(theta_hist, 1/tau)

- tau controls borrowing: large tau = strong borrowing, small tau = weak borrowing
- Empirical Bayes or fully Bayesian estimation of tau
- Limitation: designed for one-to-one (single historical study) borrowing; extension to multiple sources requires K commensurability parameters, losing interpretability

**Key reference**: Hobbs BP, Carlin BP, Mandrekar SJ, Sargent DJ. Hierarchical commensurate and power prior models for adaptive incorporation of historical information in clinical trials. *Biometrics* 2011;67:1047--1056.

### 1.4 LEAP: Latent Exchangeability Prior (Alt et al. 2024)

Models individual-level exchangeability: each subject in the historical data is classified as exchangeable or non-exchangeable with the current data:

    z_i ~ Bernoulli(pi_ex)  (latent exchangeability indicator)

- If z_i = 1: subject i's data comes from the same model as the current trial
- If z_i = 0: subject i is non-exchangeable (different distribution)
- Requires patient-level historical data (not just summary statistics)
- Limitation: fundamentally patient-level; does not apply to trial-level summary data without modification; does not incorporate structural knowledge about WHY trials may be exchangeable

**Key reference**: Alt EM, Chang X, Jiang X, Liu Q, Mo M, Xia HA, Ibrahim JG. LEAP: the latent exchangeability prior for borrowing information from historical data. *Biometrics* 2024;80:ujae083.

### 1.5 Elastic Prior (Jiang, Nie & Yuan 2023)

Defines an elastic function that monotonically maps a congruence measure between historical and current data to a borrowing weight:

    a0 = f(d(D_current, D_hist))

where f satisfies pre-specified information-borrowing constraints. The method can be pre-specified but the congruence function is computed from observed data (data-driven discounting).

**Key reference**: Jiang Y, Nie L, Yuan Y. Elastic priors to dynamically borrow information from historical data in clinical trials. *Biometrics* 2023;79:49--60.

### 1.6 Multisource Exchangeability Models (MEM) (Kaizer et al. 2018; Hobbs & Landin 2018)

Enumerates all 2^{K-1} exchangeability configurations for K baskets and model-averages:

    pi(theta_k | data) = sum_{gamma} w_gamma * pi(theta_k | data, gamma)

where gamma is a binary vector indicating which baskets are exchangeable.

- iMEM (Kaizer et al. 2021): iterative version for sequential platform trials
- dMEM (Ji et al. 2022): dynamic version with temporal weighting
- Limitation: combinatorial explosion (2^{K-1} configurations); exchangeability is binary (yes/no), not graded; does not use external structural information

**Key reference**: Kaizer AM, Koopmeiners JS, Hobbs BP. Bayesian hierarchical modeling based on multisource exchangeability. *Biostatistics* 2018;19:169--184.

### 1.7 Local Power Prior for Basket Trials (Zhou et al. 2025)

Recent proposal combining power priors with pairwise similarity assessments:

    pi(theta_k | data) propto prod_j L(theta_k | D_j)^{a_{kj}}

where a_{kj} reflects similarity between baskets k and j.

**Key reference**: Zhou T, Ji Y, Bhatt DL. A Bayesian basket trial design using local power prior. *Biometrical J.* 2025;67:e70069.

### 1.8 Bayesian Local Exchangeability (BLE) (Liu et al. 2022)

Only allows borrowing among "locally exchangeable" baskets -- those with similar response rates:

    Borrowing restricted to baskets where |p_hat_k - p_hat_j| < delta

**Key reference**: Liu R, Rong A, Bhatt DL. Bayesian local exchangeability design for phase II basket trials. *Statist. Med.* 2022;41:4367--4384.

---

## 2. Cross-Disciplinary Insights

### 2.1 Graph Diffusion Kernels (Kondor & Lafferty 2002)

The heat (diffusion) kernel on a graph defines a family of similarity measures parameterized by diffusion time t:

    K(t) = exp(-t * L)

where L is the graph Laplacian. For small t, only immediate neighbors contribute; for large t, information diffuses globally. This provides a principled, single-parameter family of graph-based similarity measures that smoothly interpolates between local and global structure.

**Key insight for our problem**: Rather than using raw Jaccard similarities, we can use the diffusion kernel to derive "effective similarity" that accounts for higher-order graph structure (e.g., two trials that share no drugs but target the same molecular pathway via different drugs).

### 2.2 Multi-Source Domain Adaptation (ML literature)

In transfer learning, the problem of combining information from K heterogeneous source domains to improve prediction on a target domain is structurally identical to our borrowing problem:

    w_k* = argmin_w E_target[L(f_w(x), y)]  subject to f_w = sum_k w_k * f_k

Key results:
- Source weights should be inversely proportional to the "distance" between source and target distributions
- Using uniform weights is provably suboptimal when sources differ in relevance
- The A-distance (Ben-David et al. 2010) between source and target distributions bounds the transfer risk

**Key insight**: The optimal borrowing weight for historical trial k should be a function of the structural distance between trial k and the current trial, not a single global parameter.

### 2.3 Graph Signal Processing: Spectral Filtering

In GSP, signals on graphs (e.g., response rates as signals on the trial-graph) can be filtered in the spectral domain:

    f_filtered = U * h(Lambda) * U^T * f

where U, Lambda are the eigenvectors/eigenvalues of the graph Laplacian. Low-pass filtering (h suppresses high eigenvalues) produces smooth signals over the graph -- analogous to borrowing from similar neighbors.

**Key insight**: The borrowing mechanism can be viewed as a spectral filter on the graph of trials. This suggests that the prior mean for a new trial should be a spectrally filtered version of the observed historical rates.

### 2.4 Expert Aggregation and Forecast Combination

The linear opinion pool combines K expert forecasts with weights:

    p_combined = sum_k w_k * p_k

Optimal weights (Ranjan & Gneiting 2010) depend on:
- Each expert's historical calibration
- Correlations between experts
- Relevance of each expert to the current prediction target

**Key insight**: The "experts" in our setting are historical trials. Their "relevance" is encoded by the KG. But unlike expert aggregation, we also have outcome data from the experts, so we can learn relevance empirically while using the KG as a structural prior on relevance.

### 2.5 Conditional Autoregressive (CAR) Spatial Models

CAR models (Besag 1974; Leroux et al. 1999) place a joint prior on spatially-indexed random effects where each unit's prior mean is the weighted average of its neighbors:

    phi_i | phi_{-i} ~ N(sum_j w_ij * phi_j / sum_j w_ij, sigma^2 / sum_j w_ij)

The joint distribution has precision Q = D_w - rho * W (Leroux CAR), where D_w = diag(rowSums(W)), W is the adjacency/weight matrix, and rho controls the strength of spatial correlation.

**Key insight**: This directly models our setting -- trial parameters as "spatial" random effects on a knowledge graph, with the KG providing the adjacency structure. The CAR prior naturally borrows more from structurally similar trials.

---

## 3. Gap Analysis

| Requirement | Power Prior | rMAP | Commensurate | LEAP | MEM | CAR |
|-------------|:-----------:|:----:|:------------:|:----:|:---:|:---:|
| Trial-specific borrowing weights | Partial | No | No | Yes* | Yes | Yes |
| Uses structural knowledge (KG) | No | No | No | No | No | Possible |
| Works with summary statistics | Yes | Yes | Yes | No | Yes | Yes |
| Robustness to irrelevant trials | Manual a0 | w_rob | tau | Automatic | Model avg | rho |
| Pre-specifiable | If a0 fixed | Yes | No | No | No | If W fixed |
| Computational tractability | Fast | Fast | MCMC | MCMC | 2^K | Fast |
| Interpretability | High | High | Medium | Medium | Low | High |

*LEAP provides patient-level, not trial-level, specific weights.

### Critical Gap
No existing method combines:
1. **Graph-informed trial-specific borrowing** (what the KG offers)
2. **Robustness to prior-data conflict** (what rMAP offers)
3. **Pre-specifiability** (what regulatory practice demands)
4. **Computational tractability** (what practical use demands)
5. **Summary-statistic compatibility** (what available data provides)

---

## 4. Synthesis: Design Principles for the New Method

1. **The KG should determine the STRUCTURE of borrowing** (which trials, how much) -- not just provide a single pooled estimate
2. **Robustness must be built in** -- a vague component or dynamic discounting that protects against misleading graph structure
3. **The borrowing weights should be pre-specifiable** from the KG alone (no peeking at outcome data) for regulatory acceptability
4. **The method should work with summary statistics** (n, y) since most historical trials only report aggregates
5. **Computation should be closed-form or near-closed-form** (conjugate or Laplace approximation) for practical adoption
6. **The framework should produce a proper prior distribution** that can be combined with the likelihood of the new trial via standard Bayesian updating

---

## 5. References

1. Ibrahim MH, Chen MH, Gwon Y, Chen F. The power prior: theory and applications. *Statist. Sci.* 2015;30:46--60.
2. Schmidli H, et al. Robust meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics* 2014;70:1023--1032.
3. Hobbs BP, et al. Hierarchical commensurate and power prior models. *Biometrics* 2011;67:1047--1056.
4. Alt EM, et al. LEAP: the latent exchangeability prior. *Biometrics* 2024;80:ujae083.
5. Jiang Y, Nie L, Yuan Y. Elastic priors to dynamically borrow information. *Biometrics* 2023;79:49--60.
6. Kaizer AM, Koopmeiners JS, Hobbs BP. Bayesian hierarchical modeling based on multisource exchangeability. *Biostatistics* 2018;19:169--184.
7. Zhou T, Ji Y, Bhatt DL. A Bayesian basket trial design using local power prior. *Biometrical J.* 2025;67:e70069.
8. Liu R, Rong A, Bhatt DL. Bayesian local exchangeability design for phase II basket trials. *Statist. Med.* 2022;41:4367--4384.
9. Kondor RI, Lafferty J. Diffusion kernels on graphs and other discrete structures. *Proc. ICML* 2002.
10. Ben-David S, et al. A theory of learning from different domains. *Machine Learning* 2010;79:151--175.
11. Besag J. Spatial interaction and the statistical analysis of lattice systems. *JRSS-B* 1974;36:192--236.
12. Leroux BG, Lei X, Bhesania N. Estimation of disease rates in small areas: a new mixed model for spatial dependence. In *Statistical Models in Epidemiology* 1999.
13. Morita S, Thall PF, Mueller P. Determining the effective sample size of a parametric prior. *Biometrics* 2008;64:595--602.
14. Ranjan R, Gneiting T. Combining probability forecasts. *JRSS-B* 2010;72:71--91.
