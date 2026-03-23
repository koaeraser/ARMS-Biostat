"""
Comparator methods for KG-CAR validation.

Implements:
1. No borrowing (independent Beta-Binomial)
2. Standard BHM (normal hierarchical model, homogeneous tau)
3. rMAP (robust meta-analytic predictive prior)
4. KG-PP (power prior with alpha = KG similarity)
5. Full pooling

All methods use analytical or empirical Bayes approximations for speed.
All predict_loo_samples methods return consistent format with:
- p_mean, theta_mean, theta_sd, ci, p_samples, theta_samples
"""

import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.special import expit, logit as sp_logit
from scipy.stats import norm, beta as beta_dist


class NoBorrowing:
    """
    Independent Beta-Binomial model for each arm.

    LOO prediction: use Jeffreys prior Beta(0.5, 0.5) — no information
    from other arms.
    """

    def predict_loo_samples(self, h_target, y, n, W=None, B=200, rng=None, **kwargs):
        if rng is None:
            rng = np.random.default_rng(42)

        # Jeffreys prior — no borrowing at all
        a_prior = 0.5
        b_prior = 0.5

        p_samples = beta_dist.rvs(a_prior, b_prior, size=B, random_state=rng)
        p_samples = np.clip(p_samples, 1e-6, 1 - 1e-6)

        p_mean = np.mean(p_samples)
        p_lo = np.percentile(p_samples, 2.5)
        p_hi = np.percentile(p_samples, 97.5)

        theta_samples = sp_logit(p_samples)
        theta_mean = np.mean(theta_samples)
        theta_sd = np.std(theta_samples)

        return {
            'p_mean': p_mean,
            'theta_mean': theta_mean,
            'theta_sd': theta_sd,
            'ci': (p_lo, p_hi),
            'p_samples': p_samples,
            'theta_samples': theta_samples,
        }


class StandardBHM:
    """
    Standard Bayesian Hierarchical Model with homogeneous tau.

    theta_h | mu, tau ~ N(mu, tau^2)
    Uses empirical Bayes to estimate (mu, tau), then predicts theta_new ~ N(mu, tau^2).
    """

    def _eb_estimates(self, theta_obs, theta_var_obs):
        """EB estimates via DerSimonian-Laird style moment matching."""
        # Weighted mean
        weights = 1.0 / (theta_var_obs + 0.01)
        mu_hat = np.average(theta_obs, weights=weights)

        # Between-study variance
        total_var = np.var(theta_obs, ddof=1) if len(theta_obs) > 1 else 0.5
        mean_samp_var = np.mean(theta_var_obs)
        tau2_hat = max(total_var - mean_samp_var, 0.05)
        tau_hat = np.sqrt(tau2_hat)

        return mu_hat, tau_hat

    def predict_loo_samples(self, h_target, y, n, W=None, B=200, rng=None, **kwargs):
        if rng is None:
            rng = np.random.default_rng(42)

        H = len(y)
        mask = np.ones(H, dtype=bool)
        mask[h_target] = False

        p_obs = (y[mask] + 0.5) / (n[mask] + 1.0)
        p_obs = np.clip(p_obs, 0.01, 0.99)
        theta_obs = sp_logit(p_obs)
        theta_var_obs = 1.0 / (n[mask] * p_obs * (1 - p_obs) + 1e-6)

        mu_hat, tau_hat = self._eb_estimates(theta_obs, theta_var_obs)

        # Predictive: theta_new ~ N(mu_hat, tau_hat^2)
        theta_samples = rng.normal(mu_hat, tau_hat, size=B)
        p_samples = expit(theta_samples)

        p_mean = np.mean(p_samples)
        p_lo = np.percentile(p_samples, 2.5)
        p_hi = np.percentile(p_samples, 97.5)

        return {
            'p_mean': p_mean,
            'theta_mean': np.mean(theta_samples),
            'theta_sd': np.std(theta_samples),
            'ci': (p_lo, p_hi),
            'p_samples': p_samples,
            'theta_samples': theta_samples,
        }


class RobustMAP:
    """
    Robust Meta-Analytic Predictive (rMAP) prior.

    Following Schmidli et al. 2014:
    1. Fit hierarchical model to historical data -> MAP prior
    2. Mix with vague: pi_rMAP = (1-w)*pi_MAP + w*pi_vague

    The vague component is N(mu_hat, sigma_vague^2) — centered on the grand
    mean but with much wider spread, following the spirit of the original paper.
    """

    def __init__(self, w_rob=0.10, vague_sd=2.0):
        self.w_rob = w_rob
        self.vague_sd = vague_sd
        self.bhm = StandardBHM()

    def predict_loo_samples(self, h_target, y, n, W=None, B=200, rng=None, **kwargs):
        if rng is None:
            rng = np.random.default_rng(42)

        H = len(y)
        mask = np.ones(H, dtype=bool)
        mask[h_target] = False

        p_obs = (y[mask] + 0.5) / (n[mask] + 1.0)
        p_obs = np.clip(p_obs, 0.01, 0.99)
        theta_obs = sp_logit(p_obs)
        theta_var_obs = 1.0 / (n[mask] * p_obs * (1 - p_obs) + 1e-6)

        mu_hat, tau_hat = self.bhm._eb_estimates(theta_obs, theta_var_obs)

        # Mixture samples
        w = self.w_rob
        n_map = int(B * (1 - w))
        n_vague = B - n_map

        theta_map = rng.normal(mu_hat, tau_hat, size=n_map)
        theta_vague = rng.normal(mu_hat, self.vague_sd, size=n_vague)

        theta_samples = np.concatenate([theta_map, theta_vague])
        rng.shuffle(theta_samples)
        p_samples = expit(theta_samples)

        p_mean = np.mean(p_samples)
        p_lo = np.percentile(p_samples, 2.5)
        p_hi = np.percentile(p_samples, 97.5)

        return {
            'p_mean': p_mean,
            'theta_mean': np.mean(theta_samples),
            'theta_sd': np.std(theta_samples),
            'ci': (p_lo, p_hi),
            'p_samples': p_samples,
            'theta_samples': theta_samples,
        }


class KGPP:
    """
    Knowledge Graph Power Prior.

    For each historical arm j, discount its contribution by alpha_j = S(target, j).
    Conjugate Beta approximation: Beta(0.5 + sum(alpha*y), 0.5 + sum(alpha*(n-y))).
    """

    def __init__(self, W):
        self.W = W
        self.H = W.shape[0]

    def predict_loo_samples(self, h_target, y, n, B=200, rng=None, **kwargs):
        if rng is None:
            rng = np.random.default_rng(42)

        mask = np.ones(self.H, dtype=bool)
        mask[h_target] = False

        alpha = self.W[h_target, mask]
        y_hist = y[mask]
        n_hist = n[mask]

        a_post = 0.5 + np.sum(alpha * y_hist)
        b_post = 0.5 + np.sum(alpha * (n_hist - y_hist))

        p_samples = beta_dist.rvs(a_post, b_post, size=B, random_state=rng)
        p_samples = np.clip(p_samples, 1e-6, 1 - 1e-6)

        p_mean = np.mean(p_samples)
        p_lo = np.percentile(p_samples, 2.5)
        p_hi = np.percentile(p_samples, 97.5)

        theta_samples = sp_logit(p_samples)

        return {
            'p_mean': p_mean,
            'theta_mean': np.mean(theta_samples),
            'theta_sd': np.std(theta_samples),
            'ci': (p_lo, p_hi),
            'p_samples': p_samples,
            'theta_samples': theta_samples,
        }


class FullPooling:
    """
    Full pooling: single parameter for all trials.
    Pool all other arms equally.
    """

    def predict_loo_samples(self, h_target, y, n, W=None, B=200, rng=None, **kwargs):
        if rng is None:
            rng = np.random.default_rng(42)

        H = len(y)
        mask = np.ones(H, dtype=bool)
        mask[h_target] = False

        a_post = 0.5 + y[mask].sum()
        b_post = 0.5 + (n[mask] - y[mask]).sum()

        p_samples = beta_dist.rvs(a_post, b_post, size=B, random_state=rng)
        p_samples = np.clip(p_samples, 1e-6, 1 - 1e-6)

        p_mean = np.mean(p_samples)
        p_lo = np.percentile(p_samples, 2.5)
        p_hi = np.percentile(p_samples, 97.5)

        theta_samples = sp_logit(p_samples)

        return {
            'p_mean': p_mean,
            'theta_mean': np.mean(theta_samples),
            'theta_sd': np.std(theta_samples),
            'ci': (p_lo, p_hi),
            'p_samples': p_samples,
            'theta_samples': theta_samples,
        }
