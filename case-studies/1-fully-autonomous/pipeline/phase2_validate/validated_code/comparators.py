"""
Comparator Methods for KG-DAP Validation
=========================================

Implements all four comparators from the methodology specification:
  1. rMAP (Robust Meta-Analytic-Predictive) — Schmidli et al. 2014
  2. Uniform prior Beta(1,1)
  3. Pooled Beta prior
  4. Equal-Weight Mixture (KG-DAP without graph weights)
"""

import numpy as np
from scipy.special import betaln, gammaln
from scipy.stats import beta as beta_dist
from scipy.optimize import minimize_scalar

# Import KG-DAP utilities for mixture operations
from kg_dap import (mixture_mean, mixture_variance, mixture_cdf,
                    mixture_quantile, mixture_interval, compute_ess,
                    log_predictive_pmf, posterior_update, log_beta_binomial_pmf,
                    construct_prior)


# ===================================================================
# Comparator 1: Uniform prior Beta(1,1)
# ===================================================================
def uniform_prior():
    """Return the uniform Beta(1,1) prior as a 1-component mixture.

    Returns
    -------
    prior : dict
        Single-component mixture: Beta(1,1).
    """
    return {
        'weights': np.array([1.0]),
        'alphas': np.array([1.0]),
        'betas': np.array([1.0])
    }


# ===================================================================
# Comparator 2: Pooled Beta prior
# ===================================================================
def pooled_prior(n_hist, y_hist, alpha0=1.0, beta0=1.0):
    """Pool all historical data into a single Beta posterior.

    pi_pool(theta) = Beta(alpha0 + sum(y_k), beta0 + sum(n_k - y_k))

    Parameters
    ----------
    n_hist : array-like
        Historical sample sizes.
    y_hist : array-like
        Historical responders.
    alpha0, beta0 : float
        Base prior parameters.

    Returns
    -------
    prior : dict
        Single-component mixture.
    """
    a = alpha0 + np.sum(y_hist)
    b = beta0 + np.sum(n_hist - y_hist)
    return {
        'weights': np.array([1.0]),
        'alphas': np.array([a]),
        'betas': np.array([b])
    }


# ===================================================================
# Comparator 3: Equal-Weight Mixture
# ===================================================================
def equal_weight_prior(n_hist, y_hist, w0=0.20, alpha0=1.0, beta0=1.0,
                       n_cap=200):
    """KG-DAP mixture but with uniform omega_k = 1/K.

    Same architecture as KG-DAP but without graph-derived weights.

    Parameters
    ----------
    n_hist, y_hist : array-like
        Historical data.
    w0 : float
        Robustness weight.
    alpha0, beta0 : float
        Base prior parameters.
    n_cap : float
        Sample size cap.

    Returns
    -------
    prior : dict
        (K+1)-component Beta mixture with equal borrowing weights.
    """
    K = len(n_hist)
    omega_equal = np.ones(K) / K
    return construct_prior(omega_equal, np.asarray(n_hist, dtype=float),
                          np.asarray(y_hist, dtype=float),
                          w0=w0, alpha0=alpha0, beta0=beta0, n_cap=n_cap)


# ===================================================================
# Comparator 4: rMAP (Robust Meta-Analytic-Predictive)
# ===================================================================
def rmap_prior(n_hist, y_hist, w_rob=0.20, tau_prior_scale=1.0):
    """Robust Meta-Analytic-Predictive (rMAP) prior.

    Implementation following Schmidli et al. (2014):
    1. Logit-transform historical rates
    2. Fit normal-normal hierarchical model: logit(p_k) ~ N(mu, tau^2)
    3. Derive MAP predictive: logit(theta) ~ N(mu_hat, sigma^2_pred)
    4. Mix with vague component

    We use the analytical Laplace/moment-matching approach:
    - mu_hat = weighted mean of logit(p_k)
    - tau_hat^2 = DerSimonian-Laird estimator
    - sigma^2_pred = tau_hat^2 + 1/K * sum(1/V_k) where V_k = within-study variance

    Then we approximate the logit-normal MAP predictive with a Beta distribution
    via moment-matching, and robustify by mixing with Beta(1,1).

    Parameters
    ----------
    n_hist : array-like
        Historical sample sizes.
    y_hist : array-like
        Historical responders.
    w_rob : float
        Weight for vague (robustness) component.
    tau_prior_scale : float
        Scale for half-normal prior on tau (not used in DL estimator).

    Returns
    -------
    prior : dict
        2-component Beta mixture: (1-w_rob)*MAP + w_rob*Uniform.
    """
    n_hist = np.asarray(n_hist, dtype=float)
    y_hist = np.asarray(y_hist, dtype=float)
    K = len(n_hist)

    # Compute observed rates (add 0.5 continuity correction for 0/1 rates)
    p_hat = (y_hist + 0.5) / (n_hist + 1.0)

    # Logit transform
    logit_p = np.log(p_hat / (1.0 - p_hat))

    # Within-study variances on logit scale (delta method)
    V_k = 1.0 / (n_hist * p_hat * (1.0 - p_hat))

    # Fixed-effect estimate (inverse-variance weighted)
    w_fe = 1.0 / V_k
    mu_fe = np.sum(w_fe * logit_p) / np.sum(w_fe)

    # DerSimonian-Laird between-study variance estimator
    Q = np.sum(w_fe * (logit_p - mu_fe)**2)
    c = np.sum(w_fe) - np.sum(w_fe**2) / np.sum(w_fe)
    tau2_dl = max(0.0, (Q - (K - 1)) / c)

    # Random-effects weights
    w_re = 1.0 / (V_k + tau2_dl)
    mu_re = np.sum(w_re * logit_p) / np.sum(w_re)

    # Predictive variance on logit scale
    # Var_pred = tau^2 + 1/sum(w_re)  [uncertainty in mu]
    sigma2_pred = tau2_dl + 1.0 / np.sum(w_re)

    # Map logit-normal predictive to Beta via moment-matching
    # E[theta] and Var[theta] from logit-normal(mu_re, sigma2_pred)
    # Using the approximation: E[expit(X)] ≈ expit(mu / sqrt(1 + c^2 * sigma^2))
    # where c = 16*sqrt(3)/(15*pi) ≈ 0.5878
    c_approx = 16.0 * np.sqrt(3.0) / (15.0 * np.pi)  # ≈ 0.5878
    kappa = np.sqrt(1.0 + c_approx**2 * sigma2_pred)
    mu_pred = 1.0 / (1.0 + np.exp(-mu_re / kappa))

    # For variance, use the full moment formula:
    # Var[theta] ≈ mu_pred * (1-mu_pred) * (1 - 1/kappa2)
    # where kappa2 = sqrt(1 + c^2 * 2*sigma^2)
    kappa2 = np.sqrt(1.0 + c_approx**2 * 2.0 * sigma2_pred)
    mu2_pred = 1.0 / (1.0 + np.exp(-mu_re / kappa2))
    # E[theta^2] ≈ E[theta^2] ≈ mu2_pred (for the second moment)
    # A simpler, more robust variance approximation:
    var_pred = mu_pred * (1.0 - mu_pred) * (1.0 - 1.0 / kappa)

    # Moment-match to Beta(a, b)
    if var_pred <= 0 or var_pred >= mu_pred * (1.0 - mu_pred):
        # Very dispersed — default to near-uniform
        a_map = 1.0
        b_map = 1.0
    else:
        nu = mu_pred * (1.0 - mu_pred) / var_pred - 1.0
        if nu <= 0:
            a_map = 1.0
            b_map = 1.0
        else:
            a_map = max(mu_pred * nu, 0.5)
            b_map = max((1.0 - mu_pred) * nu, 0.5)

    # Robustify: (1-w_rob)*Beta(a_map, b_map) + w_rob*Beta(1,1)
    return {
        'weights': np.array([w_rob, 1.0 - w_rob]),
        'alphas': np.array([1.0, a_map]),
        'betas': np.array([1.0, b_map])
    }
