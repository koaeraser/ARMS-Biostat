"""
Compute ELIR ESS (Expected Local Information Ratio) for KG-CAR.

Implements the ESS definition from Neuenschwander, Weber, Schmidli, O'Hagan (2020),
"Predictively consistent prior effective sample sizes," Biometrics 76(2):578-587.

For a prior pi(theta) on the logit scale with binomial data:
  ESS_ELIR = E_{theta ~ pi} [ i_pi(theta) / i_1(theta) ]

where:
  i_pi(theta) = -d^2/dtheta^2 log pi(theta)  (prior information)
  i_1(theta) = expit(theta)(1 - expit(theta))  (Fisher info for one Bernoulli obs)

For our Normal prior N(mu_pred, sigma^2_pred):
  i_pi(theta) = 1/sigma^2_pred  (constant)

So: ESS_ELIR = (1/sigma^2_pred) * E_{theta ~ N(mu, sigma^2)} [1 / (expit(theta)(1-expit(theta)))]

The expectation is computed by Monte Carlo.
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.special import expit

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent.parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR.parent / "validated_results"

sys.path.insert(0, str(SCRIPT_DIR))

from method import KGCAR, build_adjacency_and_laplacian
from run_validation import load_data


def elir_ess_normal_binomial(mu_pred, sigma2_pred, M=100_000, rng=None):
    """
    Compute ELIR ESS for a Normal(mu_pred, sigma2_pred) prior on logit(p)
    with Binomial(1, p) data (single Bernoulli observation).

    Parameters
    ----------
    mu_pred : float
        Prior mean on the logit scale.
    sigma2_pred : float
        Prior variance on the logit scale.
    M : int
        Number of Monte Carlo samples for the expectation.
    rng : numpy.random.Generator

    Returns
    -------
    ess_elir : float
        ELIR effective sample size.
    """
    if rng is None:
        rng = np.random.default_rng(12345)

    sigma_pred = np.sqrt(sigma2_pred)
    theta_samples = rng.normal(mu_pred, sigma_pred, size=M)

    # Fisher information for one Bernoulli observation: p(1-p)
    p = expit(theta_samples)
    fisher_info = p * (1.0 - p)

    # Avoid division by zero at extremes
    fisher_info = np.clip(fisher_info, 1e-15, None)

    # Prior information (constant for Normal): 1/sigma^2
    prior_info = 1.0 / sigma2_pred

    # Information ratio r(theta) = i_pi / i_1
    r_theta = prior_info / fisher_info

    # ELIR ESS = E[r(theta)]
    ess_elir = np.mean(r_theta)

    return ess_elir


def compute_elir_for_all_folds(y, n, p_hat, trial_ids, W, L):
    """
    For each LOO fold, fit KG-CAR and compute ELIR ESS.

    Returns DataFrame with per-fold results.
    """
    H = len(y)
    model = KGCAR(W, L, w_rob=0.10)
    rng = np.random.default_rng(42)
    records = []

    print("  Computing ELIR ESS for each of 40 LOO folds...")
    for h in range(H):
        # Fit model on H-1 arms
        mask = np.ones(H, dtype=bool)
        mask[h] = False

        model_copy = KGCAR(W, L, w_rob=model.w_rob)
        model_copy.fit(y, n, mask=mask)

        rho = model_copy.rho_hat
        sigma_phi = model_copy.sigma_phi_hat
        sigma_eps = model_copy.sigma_eps_hat

        W_row = W[h, mask]
        w_sum = W_row.sum()
        denom = rho * w_sum + (1.0 - rho)

        # CAR predictive mean: mu + phi_pred
        eta_hist = model_copy.eta_hat[mask]
        if denom < 1e-10:
            phi_pred = 0.0
        else:
            phi_pred = rho * np.sum(W_row * eta_hist) / denom
        mu_pred = model_copy.mu_hat + phi_pred

        # Predictive variance
        var_phi_pred = sigma_phi**2 / max(denom, 1e-10)
        sigma2_pred = var_phi_pred + sigma_eps**2

        # ELIR ESS
        ess_elir = elir_ess_normal_binomial(mu_pred, sigma2_pred, M=200_000, rng=rng)

        # Variance decomposition
        frac_structured = sigma_phi**2 / (sigma_phi**2 + sigma_eps**2)

        records.append({
            'trial_id': trial_ids[h],
            'arm_idx': h,
            'n_h': n[h],
            'p_hat': p_hat[h],
            'ess_elir': ess_elir,
            'mu_pred': mu_pred,
            'sigma2_pred': sigma2_pred,
            'rho': rho,
            'sigma_phi': sigma_phi,
            'sigma_eps': sigma_eps,
            'frac_structured': frac_structured,
            'w_sum': w_sum,
        })

        if (h + 1) % 10 == 0:
            print(f"    Fold {h+1}/{H} done (ESS_ELIR = {ess_elir:.1f})")

    df = pd.DataFrame(records)
    return df


def main():
    print("=" * 60)
    print("ELIR ESS COMPUTATION FOR KG-CAR")
    print("=" * 60)

    # Load data
    trials, y, n, p_hat, trial_ids, S_composite, J_drug, J_target, S_pop = load_data()
    W, L = build_adjacency_and_laplacian(S_composite)

    # Compute ELIR ESS
    ess_df = compute_elir_for_all_folds(y, n, p_hat, trial_ids, W, L)

    # Save
    out_path = RESULTS_DIR / 'ess_elir_analysis.csv'
    ess_df.to_csv(out_path, index=False)

    # Summary
    print(f"\n  === ELIR ESS Summary ===")
    print(f"  ESS_ELIR: mean={ess_df['ess_elir'].mean():.1f}, "
          f"min={ess_df['ess_elir'].min():.1f}, "
          f"max={ess_df['ess_elir'].max():.1f}")
    print(f"  sigma2_pred: mean={ess_df['sigma2_pred'].mean():.4f}")
    print(f"  rho: mean={ess_df['rho'].mean():.3f}")
    print(f"  sigma_phi: mean={ess_df['sigma_phi'].mean():.3f}")
    print(f"  sigma_eps: mean={ess_df['sigma_eps'].mean():.3f}")
    print(f"  frac_structured: mean={ess_df['frac_structured'].mean():.3f}")
    print(f"  w_sum: mean={ess_df['w_sum'].mean():.1f}")
    print(f"\n  Results saved to: {out_path}")

    # Per-fold details
    print(f"\n  Per-fold ESS_ELIR:")
    for _, row in ess_df.iterrows():
        print(f"    {row['trial_id']:>20s}: ESS_ELIR={row['ess_elir']:6.1f}, "
              f"sigma2_pred={row['sigma2_pred']:.4f}, rho={row['rho']:.3f}")


if __name__ == "__main__":
    main()
