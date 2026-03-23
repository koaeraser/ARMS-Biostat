"""
Structured Design OC v2: KG-CAR vs rMAP with shuffled similarities.

Correct protocol (fixes v1 flaws):
1. Fix one arm h as "current trial," all 40 as historical
2. Per replicate: shuffle similarities, compute KG-CAR center, generate data
3. Calibrate ONE γ* per method such that simulated T1E = 0.05
4. Power = rejection rate at that γ* under H1

Key efficiency: analytical grid quadrature + Beta CDF (no MC sampling).
"""

import argparse
import time
import pickle
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.special import expit, logit as sp_logit
from scipy.stats import beta as beta_dist, norm

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

DATA_DIR = SCRIPT_DIR.parent.parent.parent / "data"
RESULTS_DIR = SCRIPT_DIR.parent / "validated_results"


# ── Analytical posterior probability ──────────────────────────────────────


def compute_posterior_prob(y_T, n_T, y_C, n_C, mu_prior, sigma2_prior,
                           w_rob=0.10, vague_sd=2.0, n_grid=500):
    """
    P(p_T > p_C | y_T, y_C) via grid quadrature + Beta CDF.

    Treatment: Beta(0.5 + y_T, 0.5 + n_T - y_T)  [Jeffreys]
    Control:   mixture-Normal prior on logit scale, grid quadrature posterior.

    Returns scalar probability, exact up to grid resolution.
    """
    # Control posterior on logit-scale grid
    sigma = np.sqrt(sigma2_prior)
    max_sd = max(sigma, vague_sd)
    lo = mu_prior - 6 * max_sd
    hi = mu_prior + 6 * max_sd
    theta_grid = np.linspace(lo, hi, n_grid)

    log_map = norm.logpdf(theta_grid, mu_prior, sigma)
    log_vague = norm.logpdf(theta_grid, mu_prior, vague_sd)
    log_prior = np.logaddexp(
        np.log(1.0 - w_rob) + log_map,
        np.log(w_rob) + log_vague,
    )

    p_grid = expit(theta_grid)
    p_grid = np.clip(p_grid, 1e-15, 1 - 1e-15)
    log_lik = y_C * np.log(p_grid) + (n_C - y_C) * np.log(1 - p_grid)

    log_post = log_prior + log_lik
    log_post -= log_post.max()
    w_C = np.exp(log_post)
    w_C /= w_C.sum()

    # Treatment posterior: Beta(a_T, b_T)
    a_T = 0.5 + y_T
    b_T = 0.5 + (n_T - y_T)

    # P(p_T > p_C) = sum_j w_C_j * [1 - F_Beta(p_C_j; a_T, b_T)]
    prob_T_gt = 1.0 - beta_dist.cdf(p_grid, a_T, b_T)
    return float(np.dot(w_C, prob_T_gt))


def compute_prob_matrix(n_T, n_C, mu_prior, sigma2_prior,
                         w_rob=0.10, vague_sd=2.0, n_grid=500):
    """
    Precompute P[y_T, y_C] for all outcomes. Used for rMAP (fixed center).
    """
    P = np.zeros((n_T + 1, n_C + 1))
    for y_T in range(n_T + 1):
        a_T = 0.5 + y_T
        b_T = 0.5 + (n_T - y_T)
        for y_C in range(n_C + 1):
            P[y_T, y_C] = compute_posterior_prob(
                y_T, n_T, y_C, n_C, mu_prior, sigma2_prior,
                w_rob, vague_sd, n_grid,
            )
    return P


# ── Shuffled KG center ───────────────────────────────────────────────────


def shuffle_center(theta_obs, similarities, power_k, top_m, self_idx, rng):
    """
    Shuffle non-self similarities and compute weighted center.

    Returns logit-scale center.
    """
    H = len(theta_obs)
    mask = np.ones(H, dtype=bool)
    mask[self_idx] = False

    # Shuffle non-self similarity values
    non_self_sims = similarities[mask].copy()
    rng.shuffle(non_self_sims)

    # Build weight vector: power transform
    w = np.zeros(H)
    w[mask] = non_self_sims ** power_k

    # Sparsify: keep top-m
    if top_m < H - 1:
        w_active = w[mask]
        if len(w_active) > top_m:
            threshold = np.sort(w_active)[::-1][top_m]
            w[w < threshold] = 0
            w[self_idx] = 0

    denom = w.sum()
    if denom > 0:
        return np.dot(w, theta_obs) / denom
    else:
        return np.mean(theta_obs)


# ── Calibration ──────────────────────────────────────────────────────────


def calibrate_gamma(probs, target=0.05):
    """Find γ such that mean(probs > γ) ≈ target."""
    sorted_p = np.sort(probs)[::-1]
    n = len(sorted_p)
    target_count = int(np.ceil(target * n))
    if target_count <= 0:
        return 1.0
    if target_count >= n:
        return 0.0
    return sorted_p[target_count - 1] - 1e-10


# ── Noise variance from unshuffled config ────────────────────────────────


def compute_unshuffled_noise(S_composite, theta_obs, power_k, top_m):
    """
    Compute sigma2_noise = Var(theta_obs - unshuffled_center) across all arms.
    """
    H = len(theta_obs)
    W = S_composite.copy() ** power_k
    np.fill_diagonal(W, 0)

    if top_m < H - 1:
        for h in range(H):
            row = W[h, :].copy()
            row[h] = -1
            threshold_val = np.sort(row)[::-1][top_m]
            W[h, row < threshold_val] = 0
            W[h, h] = 0

    centers = np.zeros(H)
    for h in range(H):
        w_h = W[h, :]
        denom = w_h.sum()
        if denom > 0:
            centers[h] = np.dot(w_h, theta_obs) / denom
        else:
            centers[h] = np.mean(theta_obs)

    residuals = theta_obs - centers
    return residuals.var(), centers


# ── Main simulation ──────────────────────────────────────────────────────


def run_v2(arm_idx, S_composite, theta_obs, power_k, top_m,
           n_T=50, n_C=25, sigma2_kg=1.00, sigma2_rm=1.06,
           delta_grid=None, R=200, w_rob=0.10, vague_sd=2.0,
           n_grid=500, seed=42):
    """
    v2 structured design OC for one arm, one (k, m) config.
    """
    if delta_grid is None:
        delta_grid = np.array([0.00, 0.10, 0.15, 0.20])

    H = len(theta_obs)
    rng = np.random.default_rng(seed)

    # Similarity vector for this arm
    sims = S_composite[arm_idx, :].copy()

    # Noise variance from unshuffled config
    sigma2_noise, centers_unshuffled = compute_unshuffled_noise(
        S_composite, theta_obs, power_k, top_m)

    # rMAP center (invariant to shuffling)
    mu_rmap = np.mean(theta_obs)

    print(f"\n{'='*60}")
    print(f"v2 STRUCTURED DESIGN OC")
    print(f"  Arm {arm_idx}, config (k={power_k}, m={top_m}), R={R}")
    print(f"  sigma2_noise = {sigma2_noise:.4f}")
    print(f"  Unshuffled KG center for arm {arm_idx}: "
          f"{centers_unshuffled[arm_idx]:.4f} "
          f"(p = {expit(centers_unshuffled[arm_idx]):.3f})")
    print(f"  rMAP center: {mu_rmap:.4f} (p = {expit(mu_rmap):.3f})")
    print(f"  sigma2_kg={sigma2_kg}, sigma2_rm={sigma2_rm}")
    print(f"{'='*60}")

    # ── Precompute rMAP probability matrix ──
    print(f"\n  Precomputing rMAP P matrix ({n_T+1}×{n_C+1})...", end=" ",
          flush=True)
    t0 = time.time()
    P_rmap = compute_prob_matrix(n_T, n_C, mu_rmap, sigma2_rm,
                                  w_rob, vague_sd, n_grid)
    print(f"done ({time.time() - t0:.1f}s)")

    # ── Simulate per delta ──
    results = {}
    kg_centers_record = []

    for delta in delta_grid:
        print(f"\n  Δ = {delta:.2f}: simulating R={R}...", end=" ", flush=True)
        t0 = time.time()

        prob_kg = np.zeros(R)
        prob_rm = np.zeros(R)

        for r in range(R):
            # 1. Shuffle similarities → KG-CAR center
            mu_kg = shuffle_center(theta_obs, sims, power_k, top_m,
                                   arm_idx, rng)

            if delta == 0.0:
                kg_centers_record.append(mu_kg)

            # 2. True control rate
            eps = rng.normal(0, np.sqrt(sigma2_noise))
            true_theta_C = mu_kg + eps
            true_pC = expit(true_theta_C)

            # 3. True treatment rate
            true_pT = np.clip(true_pC + delta, 0.0, 0.999)

            # 4. Generate data
            y_T = rng.binomial(n_T, true_pT)
            y_C = rng.binomial(n_C, true_pC)

            # 5. Posterior probabilities
            prob_kg[r] = compute_posterior_prob(
                y_T, n_T, y_C, n_C, mu_kg, sigma2_kg,
                w_rob, vague_sd, n_grid,
            )
            prob_rm[r] = P_rmap[y_T, y_C]

        results[delta] = {'prob_kg': prob_kg, 'prob_rm': prob_rm}
        print(f"done ({time.time() - t0:.1f}s)")

    # ── Center distribution summary ──
    kg_centers = np.array(kg_centers_record)
    print(f"\n  KG-CAR center distribution (H0 shuffles):")
    print(f"    mean={kg_centers.mean():.4f}, std={kg_centers.std():.4f}")
    print(f"    range=[{kg_centers.min():.4f}, {kg_centers.max():.4f}]")
    print(f"    rMAP center (fixed): {mu_rmap:.4f}")

    # ── Calibrate γ* ──
    prob_kg_h0 = results[0.0]['prob_kg']
    prob_rm_h0 = results[0.0]['prob_rm']

    gamma_kg = calibrate_gamma(prob_kg_h0, target=0.05)
    gamma_rm = calibrate_gamma(prob_rm_h0, target=0.05)

    t1e_kg = np.mean(prob_kg_h0 > gamma_kg)
    t1e_rm = np.mean(prob_rm_h0 > gamma_rm)

    print(f"\n  Calibrated thresholds (target T1E = 0.05):")
    print(f"    KG-CAR: γ* = {gamma_kg:.4f}, verified T1E = {t1e_kg:.4f}")
    print(f"    rMAP:   γ* = {gamma_rm:.4f}, verified T1E = {t1e_rm:.4f}")

    # ── Power ──
    print(f"\n  Size-adjusted power at calibrated γ*:")
    print(f"  {'Δ':>6s}  {'KG-CAR':>8s}  {'rMAP':>8s}  {'Diff':>8s}")
    print(f"  {'-'*36}")

    records = []
    for delta in delta_grid:
        power_kg = np.mean(results[delta]['prob_kg'] > gamma_kg)
        power_rm = np.mean(results[delta]['prob_rm'] > gamma_rm)
        diff = power_kg - power_rm
        label = "T1E" if delta == 0.0 else "Power"
        print(f"  {delta:6.2f}  {power_kg:8.4f}  {power_rm:8.4f}  "
              f"{diff:+8.4f}  ({label})")
        records.append({
            'delta': delta,
            'kg_rate': power_kg,
            'rm_rate': power_rm,
            'diff': diff,
            'gamma_kg': gamma_kg,
            'gamma_rm': gamma_rm,
        })

    # ── MC standard errors ──
    print(f"\n  MC standard errors (R={R}):")
    for delta in delta_grid:
        p_kg = np.mean(results[delta]['prob_kg'] > gamma_kg)
        p_rm = np.mean(results[delta]['prob_rm'] > gamma_rm)
        se_kg = np.sqrt(p_kg * (1 - p_kg) / R)
        se_rm = np.sqrt(p_rm * (1 - p_rm) / R)
        print(f"    Δ={delta:.2f}: SE_KG={se_kg:.4f}, SE_rMAP={se_rm:.4f}")

    return pd.DataFrame(records), results, gamma_kg, gamma_rm


# ── Fixed-truth comparison ────────────────────────────────────────────────


def run_fixed_truth(true_pC_values, mu_rm, sigma2_kg=1.00, sigma2_rm=1.06,
                    n_T=50, n_C=10, threshold=0.90, delta_grid=None,
                    R=2000, w_rob=0.10, vague_sd=2.0, n_grid=500, seed=42):
    """
    Compare KG-CAR vs rMAP at fixed truth values with a FIXED threshold.

    KG-CAR: prior centered at logit(true_pC) — ideal centering.
    rMAP:   prior centered at mu_rm (grand mean) — fixed.

    Both use the same threshold γ = 0.90. No recalibration.
    When true_pC < grand mean, rMAP becomes conservative → loses power.
    """
    if delta_grid is None:
        delta_grid = np.array([0.00, 0.05, 0.10, 0.15, 0.20])

    print(f"\n{'='*60}")
    print(f"FIXED-TRUTH DESIGN OC (threshold = {threshold})")
    print(f"  KG-CAR: prior centered at truth (ideal)")
    print(f"  rMAP:   prior centered at {expit(mu_rm):.3f} "
          f"(logit = {mu_rm:.3f})")
    print(f"  n_T={n_T}, n_C={n_C}, R={R}")
    print(f"{'='*60}")

    # Precompute rMAP P matrix (fixed center)
    print(f"\n  Computing rMAP P matrix...", end=" ", flush=True)
    t0 = time.time()
    P_rm = compute_prob_matrix(n_T, n_C, mu_rm, sigma2_rm,
                                w_rob, vague_sd, n_grid)
    print(f"done ({time.time() - t0:.1f}s)")

    all_records = []

    # Header
    print(f"\n  {'pC':>5s}  {'Δ':>5s}  {'KG T1E/Pw':>10s}  "
          f"{'rMAP T1E/Pw':>12s}  {'Diff':>8s}")
    print(f"  {'-'*50}")

    for true_pC in true_pC_values:
        mu_kg = sp_logit(true_pC)

        # KG-CAR P matrix (centered at truth)
        P_kg = compute_prob_matrix(n_T, n_C, mu_kg, sigma2_kg,
                                    w_rob, vague_sd, n_grid)

        rng = np.random.default_rng(seed)

        for delta in delta_grid:
            true_pT = np.clip(true_pC + delta, 0.0, 0.999)
            y_T = rng.binomial(n_T, true_pT, size=R)
            y_C = rng.binomial(n_C, true_pC, size=R)

            rate_kg = np.mean(P_kg[y_T, y_C] > threshold)
            rate_rm = np.mean(P_rm[y_T, y_C] > threshold)
            diff = rate_kg - rate_rm

            label = "T1E" if delta == 0.0 else ""
            print(f"  {true_pC:5.2f}  {delta:5.2f}  {rate_kg:10.4f}  "
                  f"{rate_rm:12.4f}  {diff:+8.4f}  {label}")

            all_records.append({
                'true_pC': true_pC,
                'delta': delta,
                'kg_rate': rate_kg,
                'rm_rate': rate_rm,
                'diff': diff,
                'threshold': threshold,
                'n_C': n_C,
            })

    return pd.DataFrame(all_records)


# ── Realistic centering (KG-weighted average per arm) ─────────────────────


def run_realistic_centering(true_pC_values, S_composite, theta_obs,
                             power_k, top_m, mu_rm,
                             sigma2_kg=1.00, sigma2_rm=1.06,
                             n_T=50, n_C=25, threshold=0.90,
                             delta_grid=None, R=2000,
                             w_rob=0.10, vague_sd=2.0, n_grid=500, seed=42):
    """
    Compare KG-CAR vs rMAP using ACTUAL KG-computed centers per arm.

    For each arm h, KG-CAR center = weighted average of historical logit-rates
    using arm h's KG similarities. Results averaged across all 40 arms.
    """
    if delta_grid is None:
        delta_grid = np.array([0.00, 0.05, 0.10, 0.15, 0.20])

    H = len(theta_obs)

    # Compute KG-CAR centers for all arms
    _, kg_centers = compute_unshuffled_noise(S_composite, theta_obs,
                                             power_k, top_m)

    print(f"\n{'='*60}")
    print(f"REALISTIC CENTERING DESIGN OC (threshold = {threshold})")
    print(f"  KG-CAR: arm-specific centers from KG (k={power_k}, m={top_m})")
    print(f"    center range: [{expit(kg_centers.min()):.3f}, "
          f"{expit(kg_centers.max()):.3f}]")
    print(f"    center std (logit): {kg_centers.std():.3f}")
    print(f"  rMAP: fixed at {expit(mu_rm):.3f}")
    print(f"  n_T={n_T}, n_C={n_C}, R={R}, H={H} arms")
    print(f"{'='*60}")

    # Precompute rMAP P matrix (one, fixed center)
    print(f"\n  Computing rMAP P matrix...", end=" ", flush=True)
    t0 = time.time()
    P_rm = compute_prob_matrix(n_T, n_C, mu_rm, sigma2_rm,
                                w_rob, vague_sd, n_grid)
    print(f"done ({time.time() - t0:.1f}s)")

    # Precompute KG-CAR P matrices (one per arm)
    print(f"  Computing {H} KG-CAR P matrices...", end=" ", flush=True)
    t0 = time.time()
    P_kg_all = {}
    for h in range(H):
        P_kg_all[h] = compute_prob_matrix(n_T, n_C, kg_centers[h],
                                           sigma2_kg, w_rob, vague_sd, n_grid)
    print(f"done ({time.time() - t0:.1f}s)")

    # Simulate
    all_records = []

    print(f"\n  {'pC':>5s}  {'Δ':>5s}  {'KG avg':>8s}  {'rMAP':>8s}  "
          f"{'Diff':>8s}")
    print(f"  {'-'*45}")

    for true_pC in true_pC_values:
        rng = np.random.default_rng(seed)

        for delta in delta_grid:
            true_pT = np.clip(true_pC + delta, 0.0, 0.999)
            y_T = rng.binomial(n_T, true_pT, size=R)
            y_C = rng.binomial(n_C, true_pC, size=R)

            # rMAP: same for all arms
            rate_rm = np.mean(P_rm[y_T, y_C] > threshold)

            # KG-CAR: average across all arms
            arm_rates = np.zeros(H)
            for h in range(H):
                arm_rates[h] = np.mean(P_kg_all[h][y_T, y_C] > threshold)
            rate_kg = arm_rates.mean()

            diff = rate_kg - rate_rm
            label = "T1E" if delta == 0.0 else ""
            print(f"  {true_pC:5.2f}  {delta:5.2f}  {rate_kg:8.4f}  "
                  f"{rate_rm:8.4f}  {diff:+8.4f}  {label}")

            all_records.append({
                'true_pC': true_pC,
                'delta': delta,
                'kg_rate': rate_kg,
                'rm_rate': rate_rm,
                'diff': diff,
                'kg_rate_std': arm_rates.std(),
                'threshold': threshold,
                'n_C': n_C,
                'power_k': power_k,
                'top_m': top_m,
            })

    return pd.DataFrame(all_records)


# ── Main ─────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description='Structured Design OC v2 (shuffled similarities)')
    parser.add_argument('--arm', type=int, default=0,
                        help='Arm index to use as current trial')
    parser.add_argument('--R', type=int, default=200,
                        help='Number of replicates')
    parser.add_argument('--k', type=int, default=1,
                        help='Power transformation exponent')
    parser.add_argument('--m', type=int, default=10,
                        help='Top-m neighbors')
    parser.add_argument('--nC', type=int, default=25,
                        help='Control arm sample size')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--fixed', action='store_true',
                        help='Run fixed-truth comparison (ideal centering)')
    parser.add_argument('--realistic', action='store_true',
                        help='Run realistic centering (actual KG centers)')
    args = parser.parse_args()

    # Load data
    with open(DATA_DIR / "similarity_matrices.pkl", 'rb') as f:
        sim_data = pickle.load(f)

    trials_df = pd.read_csv(DATA_DIR / "trials_data.csv")

    S_composite = sim_data['S_composite']
    trial_ids = sim_data['trial_ids']

    trials_indexed = trials_df.set_index('trial_id')
    p_obs = np.array([trials_indexed.loc[tid, 'p_hat'] for tid in trial_ids])
    theta_obs = sp_logit(p_obs)

    H = len(trial_ids)
    mu_rmap = np.mean(theta_obs)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.realistic:
        true_pC_values = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
                          0.50, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85]
        df = run_realistic_centering(
            true_pC_values=true_pC_values,
            S_composite=S_composite,
            theta_obs=theta_obs,
            power_k=args.k,
            top_m=args.m,
            mu_rm=mu_rmap,
            n_T=50, n_C=args.nC,
            R=args.R,
            seed=args.seed,
        )
        out_path = RESULTS_DIR / 'design_oc_realistic.csv'
        df.to_csv(out_path, index=False)
        print(f"\n  Results saved: {out_path}")
    elif args.fixed:
        # Fixed-truth mode: KG-CAR (ideal center) vs rMAP at several truths
        true_pC_values = [0.15, 0.20, 0.25, 0.30, 0.35, 0.50, 0.65, 0.75, 0.85]
        df = run_fixed_truth(
            true_pC_values=true_pC_values,
            mu_rm=mu_rmap,
            n_T=50, n_C=args.nC,
            R=args.R,
            seed=args.seed,
        )
        out_path = RESULTS_DIR / 'design_oc_fixed_truth.csv'
        df.to_csv(out_path, index=False)
        print(f"\n  Results saved: {out_path}")
    else:
        # Shuffled mode
        print(f"  {H} arms loaded, p_obs range: "
              f"[{p_obs.min():.3f}, {p_obs.max():.3f}]")
        print(f"  Current arm: {args.arm} "
              f"(trial {trial_ids[args.arm]}, p_obs={p_obs[args.arm]:.3f})")

        df, results, gamma_kg, gamma_rm = run_v2(
            arm_idx=args.arm,
            S_composite=S_composite,
            theta_obs=theta_obs,
            power_k=args.k,
            top_m=args.m,
            n_C=args.nC,
            R=args.R,
            seed=args.seed,
        )
        out_path = RESULTS_DIR / 'design_oc_structured_v2_test.csv'
        df.to_csv(out_path, index=False)
        print(f"\n  Results saved: {out_path}")


if __name__ == "__main__":
    main()
