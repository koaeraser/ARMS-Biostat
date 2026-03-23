"""
KG-CAR: Knowledge Graph Conditional Autoregressive model for clinical trial borrowing.

Implements:
- Leroux CAR parameterization with BYM decomposition
- Robustification mixture (KG-CAR + vague)
- Laplace approximation for fast posterior inference
- Leave-one-out prediction

Mathematical model:
  y_h | p_h ~ Binomial(n_h, p_h)
  theta_h = logit(p_h) = mu + phi_h + eps_h
  phi | sigma_phi, rho ~ N(0, sigma_phi^2 * [rho*L + (1-rho)*I]^{-1})
  eps_h | sigma_eps ~ N(0, sigma_eps^2)  iid
  Robustified: (1 - w_rob) * pi_KG_CAR + w_rob * pi_vague

For LOO prediction:
  - Fit BYM model on H-1 arms
  - Use CAR conditional to predict left-out arm's phi
  - theta_target ~ N(mu + phi_pred, var_phi + var_eps)
  - Robustify by mixing with weakly informative prior centered on grand mean
"""

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit as sp_logit
from scipy.stats import norm
import warnings


def build_adjacency_and_laplacian(S_composite, trial_ids=None, delta_min=0.0):
    """
    Build adjacency matrix W and graph Laplacian L from composite similarity.

    Parameters
    ----------
    S_composite : ndarray (H, H)
        Composite similarity matrix.
    trial_ids : list, optional
        Trial identifiers.
    delta_min : float
        Minimum similarity threshold. Entries below this are set to 0.

    Returns
    -------
    W : ndarray (H, H)
        Adjacency matrix (zero diagonal).
    L : ndarray (H, H)
        Graph Laplacian D_W - W.
    """
    H = S_composite.shape[0]
    W = S_composite.copy()
    np.fill_diagonal(W, 0.0)
    W[W < delta_min] = 0.0
    W = 0.5 * (W + W.T)
    D = np.diag(W.sum(axis=1))
    L = D - W
    return W, L


def leroux_precision(L, rho, H):
    """
    Compute Leroux precision matrix: Q(rho) = rho * L + (1 - rho) * I
    """
    return rho * L + (1.0 - rho) * np.eye(H)


class KGCAR:
    """
    KG-CAR model with Laplace approximation.

    Uses marginalized BYM: eta = phi + eps, where
    eta ~ N(0, sigma_phi^2 * Q^{-1} + sigma_eps^2 * I)
    and theta_h = mu + eta_h.

    For LOO prediction:
    1. Fit model on H-1 arms to get (mu, eta_{-h}, sigma_phi, sigma_eps, rho)
    2. Predict eta_h using the CAR conditional distribution
    3. Robustify by mixing samples from CAR predictive with weakly informative prior
    """

    def __init__(self, W, L, w_rob=0.10, rho=None, sigma_phi_init=0.5,
                 sigma_eps_init=0.5, vague_sd=2.0):
        """
        Parameters
        ----------
        W : ndarray (H, H)
            Adjacency matrix.
        L : ndarray (H, H)
            Graph Laplacian.
        w_rob : float
            Robustification weight (0 to 1).
        rho : float or None
            If None, rho is estimated. If float, rho is fixed.
        vague_sd : float
            SD of the vague component on logit scale (default 2.0, which gives
            a reasonably wide but not absurd range on probability scale).
        """
        self.W = W
        self.L = L
        self.H = W.shape[0]
        self.w_rob = w_rob
        self.rho_fixed = rho
        self.sigma_phi_init = sigma_phi_init
        self.sigma_eps_init = sigma_eps_init
        self.vague_sd = vague_sd

        # Results stored after fitting
        self.mu_hat = None
        self.eta_hat = None      # posterior mode of eta (= theta - mu)
        self.sigma_phi_hat = None
        self.sigma_eps_hat = None
        self.rho_hat = None

    def _neg_log_posterior(self, params, y_sub, n_sub, L_sub, H_eff):
        """
        Negative log posterior for the marginalized BYM model.

        params: [mu, eta_1, ..., eta_H, log_sigma_phi, log_sigma_eps, logit_rho]
        """
        mu = params[0]
        eta = params[1:1 + H_eff]
        log_sigma_phi = params[1 + H_eff]
        log_sigma_eps = params[2 + H_eff]
        sigma_phi = np.exp(log_sigma_phi)
        sigma_eps = np.exp(log_sigma_eps)

        if self.rho_fixed is not None:
            rho = self.rho_fixed
        else:
            logit_rho = params[3 + H_eff]
            rho = expit(logit_rho)

        # Leroux precision
        Q_sub = leroux_precision(L_sub, rho, H_eff)

        # Marginal covariance of eta = phi + eps:
        # Sigma = sigma_phi^2 * Q^{-1} + sigma_eps^2 * I
        # Precision: P = Sigma^{-1}
        try:
            Q_inv = np.linalg.solve(Q_sub, np.eye(H_eff))
        except np.linalg.LinAlgError:
            return 1e10

        Sigma_eta = sigma_phi**2 * Q_inv + sigma_eps**2 * np.eye(H_eff)

        try:
            # Use Cholesky for stability
            cho = np.linalg.cholesky(Sigma_eta)
            logdet_Sigma = 2.0 * np.sum(np.log(np.diag(cho)))
            z = np.linalg.solve(cho, eta)
            quad_form = np.dot(z, z)
        except np.linalg.LinAlgError:
            return 1e10

        # Log prior on eta: -0.5 * (logdet_Sigma + eta' P eta + H*log(2pi))
        log_prior_eta = -0.5 * (logdet_Sigma + quad_form + H_eff * np.log(2 * np.pi))

        # Binomial log likelihood
        theta = mu + eta
        # Use numerically stable formula
        log_lik = np.sum(y_sub * theta - n_sub * np.logaddexp(0, theta))

        # Prior on mu: N(0, 25)
        log_prior_mu = -0.5 * mu**2 / 25.0

        # Prior on sigma_phi: Half-N(0, 1)
        # On log scale: log p(log_sigma) = log(2) - sigma^2/2 + log_sigma
        # Jacobian: |d(sigma)/d(log_sigma)| = sigma = exp(log_sigma)
        log_prior_sigma_phi = np.log(2) - 0.5 * sigma_phi**2 + log_sigma_phi
        log_prior_sigma_eps = np.log(2) - 0.5 * sigma_eps**2 + log_sigma_eps

        # Prior on rho: Beta(1,1) = Uniform
        log_prior_rho = 0.0
        if self.rho_fixed is None:
            # Jacobian of logit transform
            log_prior_rho = np.log(rho * (1 - rho) + 1e-20)

        total = log_lik + log_prior_eta + log_prior_mu + log_prior_sigma_phi + log_prior_sigma_eps + log_prior_rho
        return -total

    def fit(self, y, n, mask=None):
        """
        Fit the KG-CAR model using Laplace approximation (mode finding).
        """
        if mask is None:
            mask = np.ones(self.H, dtype=bool)
        H_eff = mask.sum()

        y_sub = y[mask]
        n_sub = n[mask]
        L_sub = self.L[np.ix_(mask, mask)]

        # Initialize
        p_hat = (y_sub + 0.5) / (n_sub + 1.0)
        theta_init = sp_logit(np.clip(p_hat, 0.02, 0.98))
        mu_init = np.mean(theta_init)
        eta_init = theta_init - mu_init

        x0 = np.concatenate([
            [mu_init],
            eta_init,
            [np.log(self.sigma_phi_init)],
            [np.log(self.sigma_eps_init)],
        ])
        if self.rho_fixed is None:
            x0 = np.concatenate([x0, [0.0]])

        result = minimize(
            self._neg_log_posterior,
            x0,
            args=(y_sub, n_sub, L_sub, H_eff),
            method='L-BFGS-B',
            options={'maxiter': 3000, 'ftol': 1e-12, 'gtol': 1e-8}
        )

        if not result.success:
            # Fallback
            result2 = minimize(
                self._neg_log_posterior,
                x0,
                args=(y_sub, n_sub, L_sub, H_eff),
                method='Nelder-Mead',
                options={'maxiter': 20000, 'xatol': 1e-8, 'fatol': 1e-12}
            )
            if result2.fun < result.fun:
                result = result2

        params = result.x
        self.mu_hat = params[0]

        self.eta_hat_sub = params[1:1 + H_eff]
        self.eta_hat = np.full(self.H, np.nan)
        self.eta_hat[mask] = self.eta_hat_sub

        self.sigma_phi_hat = np.exp(params[1 + H_eff])
        self.sigma_eps_hat = np.exp(params[2 + H_eff])

        if self.rho_fixed is not None:
            self.rho_hat = self.rho_fixed
        else:
            self.rho_hat = expit(params[3 + H_eff])

        # Store theta hat for fitted arms
        self.theta_hat = np.full(self.H, np.nan)
        self.theta_hat[mask] = self.mu_hat + self.eta_hat_sub

        return self

    def predict_loo(self, h_target, y, n):
        """
        Leave-one-out prediction for arm h_target.

        1. Fit on remaining arms
        2. Predict using CAR conditional for phi + iid for eps
        3. Robustify with weakly informative mixture component
        """
        mask = np.ones(self.H, dtype=bool)
        mask[h_target] = False

        self.fit(y, n, mask=mask)

        rho = self.rho_hat
        W_row = self.W[h_target, mask]
        w_sum = W_row.sum()
        denom = rho * w_sum + (1 - rho)

        # The marginal of eta_h is eta_h = phi_h + eps_h
        # But the CAR conditional is on phi, and eps is independent.
        # For the BYM2 marginal eta, we need to think about the marginal
        # covariance structure.
        #
        # However, the fitted eta values for historical arms combine both
        # phi and eps. We approximate: use the fitted theta values to
        # extract the "effective" spatial contribution.
        #
        # The CAR conditional for the marginal BYM model:
        # For the target arm, given the fitted theta values of neighbors,
        # the predictive is approximately:
        #   theta_target ~ N(mu + phi_pred, var_phi + var_eps)
        # where phi_pred is the CAR conditional mean based on the
        # fitted eta values (treating them as phi estimates).

        eta_hist = self.eta_hat[mask]

        if denom < 1e-10:
            phi_pred = 0.0
        else:
            phi_pred = rho * np.sum(W_row * eta_hist) / denom

        # Predictive variance for structured component
        var_phi_pred = self.sigma_phi_hat**2 / max(denom, 1e-10)
        # Add unstructured variance
        var_theta_pred = var_phi_pred + self.sigma_eps_hat**2

        # KG-CAR predictive
        theta_mean_car = self.mu_hat + phi_pred
        theta_sd_car = np.sqrt(var_theta_pred)

        return theta_mean_car, theta_sd_car

    def predict_loo_samples(self, h_target, y, n, B=200, rng=None, **kwargs):
        """
        LOO prediction with Monte Carlo samples from the robustified predictive.
        """
        if rng is None:
            rng = np.random.default_rng(42)

        theta_mean_car, theta_sd_car = self.predict_loo(h_target, y, n)

        # Robustified mixture sampling
        w = self.w_rob
        n_car = int(B * (1 - w))
        n_vague = B - n_car

        # CAR component
        theta_samples_car = rng.normal(theta_mean_car, theta_sd_car, size=n_car)

        # Vague component: centered on grand mean with wider SD
        # This is the "weakly informative" prior that protects against
        # KG misspecification while still being centered appropriately.
        theta_mean_vague = self.mu_hat  # center on fitted grand mean
        theta_sd_vague = self.vague_sd   # wider spread

        theta_samples_vague = rng.normal(theta_mean_vague, theta_sd_vague, size=n_vague)

        theta_samples = np.concatenate([theta_samples_car, theta_samples_vague])
        rng.shuffle(theta_samples)
        p_samples = expit(theta_samples)

        # Point estimates from samples
        p_mean = np.mean(p_samples)
        theta_mean = np.mean(theta_samples)
        theta_sd = np.std(theta_samples)

        # CI from sample quantiles
        p_lo = np.percentile(p_samples, 2.5)
        p_hi = np.percentile(p_samples, 97.5)

        return {
            'p_mean': p_mean,
            'theta_mean': theta_mean,
            'theta_sd': theta_sd,
            'ci': (p_lo, p_hi),
            'p_samples': p_samples,
            'theta_samples': theta_samples,
            'rho_hat': self.rho_hat,
            'sigma_phi_hat': self.sigma_phi_hat,
            'sigma_eps_hat': self.sigma_eps_hat,
        }


class KGCARNoRob(KGCAR):
    """KG-CAR without robustification (ablation)."""
    def __init__(self, W, L, **kwargs):
        kwargs['w_rob'] = 0.0
        super().__init__(W, L, **kwargs)


class KGCARNoBYM(KGCAR):
    """KG-CAR without BYM decomposition (no eps, only structured).

    Forces sigma_eps to be very small so only the structured component remains.
    """
    def __init__(self, W, L, **kwargs):
        super().__init__(W, L, **kwargs)
        self.sigma_eps_init = 1e-6

    def _neg_log_posterior(self, params, y_sub, n_sub, L_sub, H_eff):
        """Override to fix sigma_eps near zero."""
        params = params.copy()
        params[2 + H_eff] = np.log(1e-6)
        return super()._neg_log_posterior(params, y_sub, n_sub, L_sub, H_eff)
