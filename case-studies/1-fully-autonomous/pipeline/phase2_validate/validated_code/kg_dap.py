"""
KG-DAP (Knowledge-Graph-Driven Adaptive Prior) Implementation
==============================================================

Implements the KG-DAP methodology as specified in methodology_specification.md:
  - Composite similarity matrix (reusing build_kg.py infrastructure)
  - Normalized graph Laplacian
  - Diffusion kernel weights
  - Beta mixture prior construction
  - Posterior update (conjugate)
  - Prior/posterior predictive moments and intervals
  - ESS computation
"""

import os
import sys
import pickle

import numpy as np
import pandas as pd
from scipy import linalg as la
from scipy.special import betaln, gammaln
from scipy.stats import beta as beta_dist

# ---------------------------------------------------------------------------
# Add project root to path so we can import build_kg
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

import build_kg


# ===================================================================
# Input validation
# ===================================================================
def validate_inputs(n_arr, y_arr, alpha, beta_w, gamma, beta_diff, w0,
                    alpha0, beta0, n_cap):
    """Validate all inputs before computation. Raises ValueError on failure."""
    # Sample sizes positive
    if np.any(n_arr <= 0):
        raise ValueError(f"All sample sizes must be positive. Got min n={n_arr.min()}")
    # Responses in valid range
    if np.any(y_arr < 0) or np.any(y_arr > n_arr):
        raise ValueError("Responses must satisfy 0 <= y <= n for all arms")
    # Similarity weights sum to 1
    if abs(alpha + beta_w + gamma - 1.0) > 1e-9:
        raise ValueError(f"Similarity weights must sum to 1. Got {alpha + beta_w + gamma}")
    # Diffusion parameter positive
    if beta_diff <= 0:
        raise ValueError(f"Diffusion parameter beta must be > 0. Got {beta_diff}")
    # Robustness weight in [0, 1)
    if w0 < 0 or w0 >= 1:
        raise ValueError(f"Robustness weight w0 must be in [0, 1). Got {w0}")
    # Base prior parameters positive
    if alpha0 <= 0 or beta0 <= 0:
        raise ValueError(f"Base prior params must be > 0. Got alpha0={alpha0}, beta0={beta0}")
    # n_cap positive
    if n_cap <= 0:
        raise ValueError(f"Sample size cap must be > 0. Got {n_cap}")


# ===================================================================
# Step 1: Composite similarity (delegating to build_kg.py)
# ===================================================================
def compute_similarity_matrix(G, trial_ids, alpha=0.20, beta_w=0.20, gamma=0.60):
    """Compute the composite similarity matrix using build_kg infrastructure.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Knowledge graph from build_kg.
    trial_ids : list
        Sorted trial identifiers.
    alpha, beta_w, gamma : float
        Similarity component weights (must sum to 1).

    Returns
    -------
    S : np.ndarray, shape (K, K)
        Composite similarity matrix with S_ii = 1, S_ij in [0, 1].
    """
    J_drug = build_kg.compute_drug_jaccard(G, trial_ids)
    J_target = build_kg.compute_target_jaccard(G, trial_ids)
    S_pop = build_kg.compute_population_similarity(G, trial_ids)
    S = build_kg.compute_composite_similarity(J_drug, J_target, S_pop,
                                              alpha=alpha, beta=beta_w,
                                              gamma=gamma)
    # Ensure diagonal is exactly 1
    np.fill_diagonal(S, 1.0)
    # Clip to [0, 1]
    S = np.clip(S, 0.0, 1.0)
    return S


# ===================================================================
# Step 2: Normalized graph Laplacian
# ===================================================================
def compute_normalized_laplacian(S):
    """Compute the normalized graph Laplacian L = I - D^{-1/2} S D^{-1/2}.

    Parameters
    ----------
    S : np.ndarray, shape (K, K)
        Symmetric similarity matrix with positive entries.

    Returns
    -------
    L : np.ndarray, shape (K, K)
        Normalized Laplacian with eigenvalues in [0, 2].
    """
    K = S.shape[0]
    d = S.sum(axis=1)  # degree vector
    # Guard against zero degree (isolated node)
    d_inv_sqrt = np.where(d > 0, 1.0 / np.sqrt(d), 0.0)
    D_inv_sqrt = np.diag(d_inv_sqrt)
    L = np.eye(K) - D_inv_sqrt @ S @ D_inv_sqrt
    # Ensure symmetry (floating point)
    L = 0.5 * (L + L.T)
    return L


# ===================================================================
# Step 3: Borrowing weight computation
# ===================================================================
def compute_diffusion_weights(S, current_idx, beta_diff=2.0):
    """Compute diffusion-kernel borrowing weights for a current trial.

    Parameters
    ----------
    S : np.ndarray, shape (K+1, K+1)
        Similarity matrix including the current trial.
    current_idx : int
        Index of the current trial in S.
    beta_diff : float
        Diffusion time parameter (higher = more concentrated).

    Returns
    -------
    omega : np.ndarray, shape (K,)
        Normalized borrowing weights for historical arms (sums to 1).
    """
    L = compute_normalized_laplacian(S)
    # Eigendecomposition (L is symmetric)
    eigenvalues, U = la.eigh(L)
    # Clip eigenvalues to [0, 2] for numerical safety
    eigenvalues = np.clip(eigenvalues, 0.0, 2.0)
    # Diffusion kernel: K_beta = U diag(exp(-beta * lambda)) U^T
    exp_eigvals = np.exp(-beta_diff * eigenvalues)
    K_beta = U @ np.diag(exp_eigvals) @ U.T

    # Raw affinities from current trial to all others
    K_total = S.shape[0]
    hist_indices = [i for i in range(K_total) if i != current_idx]
    raw_affinities = K_beta[current_idx, hist_indices]

    # Ensure non-negative (numerical safety)
    raw_affinities = np.maximum(raw_affinities, 0.0)

    # Normalize
    total = raw_affinities.sum()
    if total <= 0:
        # Fallback to uniform if all affinities are zero
        omega = np.ones(len(hist_indices)) / len(hist_indices)
    else:
        omega = raw_affinities / total

    # Assert simplex
    assert abs(omega.sum() - 1.0) < 1e-10, f"Weights don't sum to 1: {omega.sum()}"
    assert np.all(omega >= 0), "Negative weights detected"

    return omega


def compute_power_weights(S, current_idx, beta_power=10.0):
    """Compute power-sharpened borrowing weights from raw similarity.

    omega_k = S_{c,k}^beta / sum_j S_{c,j}^beta

    This is equivalent to a softmax over beta * log(S_{c,k}), providing
    a direct, monotonic mapping from similarity to borrowing weight.
    The parameter beta controls concentration: beta=0 gives uniform,
    beta->inf gives nearest-neighbor.

    Parameters
    ----------
    S : np.ndarray, shape (K+1, K+1)
        Similarity matrix including the current trial.
    current_idx : int
        Index of the current trial in S.
    beta_power : float
        Sharpness parameter (higher = more concentrated on similar trials).

    Returns
    -------
    omega : np.ndarray, shape (K,)
        Normalized borrowing weights (sums to 1).
    """
    K_total = S.shape[0]
    hist_indices = [i for i in range(K_total) if i != current_idx]
    sims = S[current_idx, hist_indices]

    # Clip small similarities to avoid 0^beta issues
    sims = np.maximum(sims, 1e-10)

    # Power sharpening
    raw = sims ** beta_power

    # Normalize
    total = raw.sum()
    if total <= 0:
        omega = np.ones(len(hist_indices)) / len(hist_indices)
    else:
        omega = raw / total

    # Assert simplex
    assert abs(omega.sum() - 1.0) < 1e-10, f"Weights don't sum to 1: {omega.sum()}"
    assert np.all(omega >= 0), "Negative weights detected"

    return omega


# ===================================================================
# Step 4 & 5 & 6: KG-DAP Prior Construction
# ===================================================================
def construct_prior(omega, n_hist, y_hist, w0=0.20, alpha0=1.0, beta0=1.0,
                    n_cap=200):
    """Construct the KG-DAP Beta mixture prior.

    Parameters
    ----------
    omega : np.ndarray, shape (K,)
        Graph-derived borrowing weights (sum to 1).
    n_hist : np.ndarray, shape (K,)
        Historical sample sizes.
    y_hist : np.ndarray, shape (K,)
        Historical responders.
    w0 : float
        Robustness weight for vague component.
    alpha0, beta0 : float
        Base prior parameters.
    n_cap : float
        Sample size cap.

    Returns
    -------
    prior : dict
        Keys: 'weights' (K+1,), 'alphas' (K+1,), 'betas' (K+1,)
        Component 0 is the vague component.
    """
    K = len(omega)

    # Mixture weights: [w0, (1-w0)*omega_1, ..., (1-w0)*omega_K]
    weights = np.zeros(K + 1)
    weights[0] = w0
    weights[1:] = (1.0 - w0) * omega

    # Beta parameters
    alphas = np.zeros(K + 1)
    betas_arr = np.zeros(K + 1)

    # Component 0: vague Beta(1,1)
    alphas[0] = 1.0
    betas_arr[0] = 1.0

    # Components 1..K: historical arms with optional sample-size cap
    for k in range(K):
        scale_k = min(1.0, n_cap / n_hist[k]) if n_cap < np.inf else 1.0
        alphas[k + 1] = alpha0 + y_hist[k] * scale_k
        betas_arr[k + 1] = beta0 + (n_hist[k] - y_hist[k]) * scale_k

    # Validate
    assert abs(weights.sum() - 1.0) < 1e-10, f"Prior weights sum={weights.sum()}"
    assert np.all(alphas > 0), "Non-positive alpha parameter"
    assert np.all(betas_arr > 0), "Non-positive beta parameter"

    return {'weights': weights, 'alphas': alphas, 'betas': betas_arr}


# ===================================================================
# Step 7: Posterior Update
# ===================================================================
def posterior_update(prior, n_c, y_c):
    """Bayesian update of the Beta mixture given observed data.

    Parameters
    ----------
    prior : dict
        From construct_prior().
    n_c : int
        Current trial sample size.
    y_c : int
        Current trial responders.

    Returns
    -------
    posterior : dict
        Same structure as prior, with updated weights and parameters.
    """
    w_prior = prior['weights']
    a_prior = prior['alphas']
    b_prior = prior['betas']
    K_plus_1 = len(w_prior)

    # Updated parameters
    a_post = a_prior + y_c
    b_post = b_prior + (n_c - y_c)

    # Updated weights via marginal likelihood
    # log(v_k) = log(w_k) + lbeta(a_k+y, b_k+n-y) - lbeta(a_k, b_k)
    log_v = np.zeros(K_plus_1)
    for k in range(K_plus_1):
        log_v[k] = (np.log(w_prior[k] + 1e-300)
                     + betaln(a_post[k], b_post[k])
                     - betaln(a_prior[k], b_prior[k]))

    # Log-sum-exp normalization
    max_log_v = log_v.max()
    log_v_shifted = log_v - max_log_v
    v_post = np.exp(log_v_shifted)
    v_post /= v_post.sum()

    return {'weights': v_post, 'alphas': a_post, 'betas': b_post}


# ===================================================================
# Mixture moments and quantiles
# ===================================================================
def mixture_mean(mixture):
    """Compute E[theta] for a Beta mixture."""
    w = mixture['weights']
    a = mixture['alphas']
    b = mixture['betas']
    component_means = a / (a + b)
    return np.dot(w, component_means)


def mixture_variance(mixture):
    """Compute Var[theta] for a Beta mixture."""
    w = mixture['weights']
    a = mixture['alphas']
    b = mixture['betas']
    ab = a + b
    component_means = a / ab
    component_vars = (a * b) / (ab**2 * (ab + 1))
    # Var = E[E[theta|k]^2] + E[Var[theta|k]] - (E[theta])^2
    mu = np.dot(w, component_means)
    var = np.dot(w, component_means**2 + component_vars) - mu**2
    return var


def mixture_cdf(mixture, x):
    """Compute CDF of Beta mixture at point x."""
    w = mixture['weights']
    a = mixture['alphas']
    b = mixture['betas']
    cdf_val = 0.0
    for k in range(len(w)):
        cdf_val += w[k] * beta_dist.cdf(x, a[k], b[k])
    return cdf_val


def mixture_quantile(mixture, q, tol=1e-8):
    """Compute quantile of Beta mixture via bisection."""
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mixture_cdf(mixture, mid) < q:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def mixture_interval(mixture, level=0.95):
    """Compute equal-tailed credible interval."""
    alpha_level = (1 - level) / 2
    lower = mixture_quantile(mixture, alpha_level)
    upper = mixture_quantile(mixture, 1 - alpha_level)
    return lower, upper


# ===================================================================
# Step 8: ESS
# ===================================================================
def compute_ess(mixture):
    """Compute effective sample size via moment-matching to a Beta."""
    mu = mixture_mean(mixture)
    var = mixture_variance(mixture)
    if var <= 0 or mu <= 0 or mu >= 1:
        return 2.0  # degenerate case
    kappa = mu * (1 - mu) / var - 1
    if kappa <= 0:
        return 2.0  # mixture is more dispersed than Beta(1,1)
    ess = kappa  # ESS = a + b = kappa + 2 for Beta; using the formula from spec: ESS = mu*(1-mu)/var - 1
    return max(ess, 0.0)


# ===================================================================
# Step 9: Predictive distribution (log-scale for numerical stability)
# ===================================================================
def log_beta_binomial_pmf(y, n, a, b):
    """Log of BetaBinomial(y | n, a, b) PMF.

    P(y|n,a,b) = C(n,y) * B(y+a, n-y+b) / B(a, b)
    """
    log_comb = gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1)
    log_bb = log_comb + betaln(y + a, n - y + b) - betaln(a, b)
    return log_bb


def log_predictive_pmf(mixture, y, n):
    """Log P(y | n, prior) for the Beta mixture prior.

    P(y|n) = sum_k w_k * BetaBin(y|n, a_k, b_k)
    """
    w = mixture['weights']
    a = mixture['alphas']
    b = mixture['betas']
    log_terms = np.zeros(len(w))
    for k in range(len(w)):
        log_terms[k] = np.log(w[k] + 1e-300) + log_beta_binomial_pmf(y, n, a[k], b[k])
    # Log-sum-exp
    max_lt = log_terms.max()
    return max_lt + np.log(np.exp(log_terms - max_lt).sum())


# ===================================================================
# Full KG-DAP pipeline for LOO-CV
# ===================================================================
def build_kg_dap_prior_loocv(G, trial_ids, trials_df, held_out_idx,
                              alpha=0.20, beta_w=0.20, gamma=0.60,
                              beta_diff=2.0, w0=0.20, alpha0=1.0, beta0=1.0,
                              n_cap=200, weight_method='power'):
    """Build KG-DAP prior for one held-out trial arm in LOO-CV.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Full knowledge graph.
    trial_ids : list
        All trial IDs (length K_total = 35).
    trials_df : pd.DataFrame
        Trial data with columns: trial_id, n, y, p_hat.
    held_out_idx : int
        Index (0-based) of the held-out arm in trial_ids.
    weight_method : str
        'diffusion' for original Laplacian diffusion kernel,
        'power' for power-sharpened raw similarity weights.
    [other hyperparameters]

    Returns
    -------
    prior : dict
        KG-DAP Beta mixture prior for the held-out trial.
    """
    K_total = len(trial_ids)

    # Get data for all arms
    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']))

    # All trial IDs including held-out
    all_ids = trial_ids

    # Compute similarity matrix for ALL arms (including held-out)
    S_full = compute_similarity_matrix(G, all_ids, alpha=alpha, beta_w=beta_w, gamma=gamma)

    # Compute borrowing weights
    if weight_method == 'diffusion':
        omega = compute_diffusion_weights(S_full, held_out_idx, beta_diff=beta_diff)
    elif weight_method == 'power':
        omega = compute_power_weights(S_full, held_out_idx, beta_power=beta_diff)
    else:
        raise ValueError(f"Unknown weight_method: {weight_method}")

    # Historical arms = all except held-out
    hist_indices = [i for i in range(K_total) if i != held_out_idx]
    n_hist = np.array([tid_to_data[all_ids[i]][0] for i in hist_indices], dtype=float)
    y_hist = np.array([tid_to_data[all_ids[i]][1] for i in hist_indices], dtype=float)

    # Validate
    validate_inputs(n_hist, y_hist, alpha, beta_w, gamma, beta_diff, w0,
                    alpha0, beta0, n_cap)

    # Construct prior
    prior = construct_prior(omega, n_hist, y_hist, w0=w0, alpha0=alpha0,
                           beta0=beta0, n_cap=n_cap)

    return prior


def build_kg_dap_prior_with_external(G, trial_ids, trials_df, held_out_idx,
                                      external_df,
                                      alpha=0.20, beta_w=0.20, gamma=0.60,
                                      beta_diff=2.0, w0=0.20, alpha0=1.0,
                                      beta0=1.0, n_cap=200,
                                      weight_method='power'):
    """Build KG-DAP prior including external (adjacent/foreign) trials.

    External trials don't have KG nodes. We assign them similarity based on
    drug overlap and population features from their CSV data, then include
    them in the mixture with diffusion-derived weights.

    For simplicity with external trials that are NOT in the KG:
    - We compute their similarity to internal trials using only population
      similarity and drug overlap as described in the CSV columns
    - External trials get a fixed low similarity to the current trial
      (0.3 for adjacent, 0.1 for foreign) since they measure different
      endpoints/indications

    Returns
    -------
    prior : dict
        KG-DAP Beta mixture prior incorporating external trials.
    """
    K_total = len(trial_ids)

    # Get data for internal arms
    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']))

    # Compute internal similarity matrix
    S_internal = compute_similarity_matrix(G, trial_ids, alpha=alpha,
                                           beta_w=beta_w, gamma=gamma)

    # Number of external trials
    n_ext = len(external_df)
    K_aug = K_total + n_ext

    # Build augmented similarity matrix
    S_aug = np.zeros((K_aug, K_aug))
    S_aug[:K_total, :K_total] = S_internal

    # Assign similarity for external trials
    # Adjacent trials (same disease, different endpoint/drug): moderate similarity
    # Foreign trials (different disease): low similarity
    for ext_idx, (_, ext_row) in enumerate(external_df.iterrows()):
        aug_idx = K_total + ext_idx
        S_aug[aug_idx, aug_idx] = 1.0  # self-similarity

        indication = ext_row.get('indication', '')
        # Determine base similarity
        if indication == 'NDMM':
            base_sim = 0.40  # adjacent: same disease
        elif indication in ('RRMM',):
            base_sim = 0.30  # adjacent: related disease
        else:
            base_sim = 0.10  # foreign: different disease

        # Set similarity to all internal trials
        for int_idx in range(K_total):
            S_aug[int_idx, aug_idx] = base_sim
            S_aug[aug_idx, int_idx] = base_sim

        # Set similarity between external trials
        for ext_idx2 in range(ext_idx):
            aug_idx2 = K_total + ext_idx2
            S_aug[aug_idx, aug_idx2] = base_sim * 0.5
            S_aug[aug_idx2, aug_idx] = base_sim * 0.5

    # Compute borrowing weights on augmented matrix
    if weight_method == 'diffusion':
        omega_aug = compute_diffusion_weights(S_aug, held_out_idx, beta_diff=beta_diff)
    elif weight_method == 'power':
        omega_aug = compute_power_weights(S_aug, held_out_idx, beta_power=beta_diff)
    else:
        raise ValueError(f"Unknown weight_method: {weight_method}")

    # Gather all historical data (internal + external)
    hist_indices_internal = [i for i in range(K_total) if i != held_out_idx]
    n_all = []
    y_all = []

    for i in hist_indices_internal:
        n_all.append(tid_to_data[trial_ids[i]][0])
        y_all.append(tid_to_data[trial_ids[i]][1])

    for _, ext_row in external_df.iterrows():
        n_all.append(int(ext_row['n']))
        y_all.append(int(ext_row['y']))

    n_hist = np.array(n_all, dtype=float)
    y_hist = np.array(y_all, dtype=float)

    # Construct prior
    prior = construct_prior(omega_aug, n_hist, y_hist, w0=w0, alpha0=alpha0,
                           beta0=beta0, n_cap=n_cap)

    return prior
