"""
Simulation study for KG-CAR.

Generates synthetic trial data with known ground truth parameters and KG
similarity matrices, then evaluates KG-CAR against comparators.

5 Scenarios:
1. Homogeneous trials (all similar)
2. Heterogeneous trials (KG correctly identifies clusters)
3. Misleading KG (high similarity but heterogeneous outcomes)
4. Contamination (foreign trials with low KG similarity)
5. Corrupted KG (randomly permuted adjacency matrix)

Usage:
    python simulation_study.py              # Quick test (R=20, B=500)
    python simulation_study.py --full       # Full run (R=1000, B=2000)
    python simulation_study.py --scenario 2 # Single scenario
"""

import sys
import os
import time
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.special import expit, logit as sp_logit
from scipy.stats import norm

# Add validated code to path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
VALIDATED_DIR = PROJECT_DIR / "pipeline" / "phase2_validate" / "validated_code"
DATA_OUT_DIR = SCRIPT_DIR / "data"

sys.path.insert(0, str(VALIDATED_DIR))

from method import KGCAR, build_adjacency_and_laplacian
from comparator import NoBorrowing, StandardBHM, RobustMAP, KGPP, FullPooling


# ---------------------------------------------------------------------------
# CRPS (copied from run_validation.py for self-containedness)
# ---------------------------------------------------------------------------

def crps_sample(y_true, samples):
    """CRPS from empirical samples: E|X-y| - 0.5*E|X-X'|."""
    n = len(samples)
    term1 = np.mean(np.abs(samples - y_true))
    sorted_samples = np.sort(samples)
    indices = np.arange(1, n + 1)
    term2 = 2 * np.sum((2 * indices - n - 1) * sorted_samples) / (n * n)
    return term1 - 0.5 * term2


# ---------------------------------------------------------------------------
# Scenario data-generating processes
# ---------------------------------------------------------------------------

def generate_scenario_1(rng, H=20, n_target=50):
    """
    Homogeneous trials: all similar rates, KG similarity high.

    - H historical + 1 target
    - True p_h ~ Beta(3,7) for all (mean ~0.3)
    - n_h ~ DiscreteUniform(30, 100)
    - W_ij ~ Uniform(0.7, 1.0) for all i != j
    """
    H_total = H + 1

    # True rates (mean ~0.3, away from 0.5 to avoid Jeffreys prior artifact)
    p_true = rng.beta(3, 7, size=H_total)

    # Sample sizes
    n = np.concatenate([
        rng.integers(30, 101, size=H),
        [n_target]
    ]).astype(float)

    # Observed data
    y = rng.binomial(n.astype(int), p_true).astype(float)

    # KG similarity: all high
    W_raw = rng.uniform(0.7, 1.0, size=(H_total, H_total))
    W_raw = 0.5 * (W_raw + W_raw.T)
    np.fill_diagonal(W_raw, 0.0)

    return p_true, y, n, W_raw


def generate_scenario_2(rng, H=20, n_target=50):
    """
    Heterogeneous trials: two clusters, KG correctly separates them.

    - 10 similar (p ~ Beta(7,3), mean ~0.7)
    - 10 dissimilar (p ~ Beta(2,8), mean ~0.2)
    - Target in similar cluster (p ~ Beta(7,3))
    - W ~0.8 within clusters, ~0.1 between
    """
    H_sim = H // 2
    H_dis = H - H_sim
    H_total = H + 1

    # True rates: similar cluster (mean ~0.7), dissimilar cluster (mean ~0.2)
    p_sim = rng.beta(7, 3, size=H_sim)
    p_dis = rng.beta(2, 8, size=H_dis)
    p_target = rng.beta(7, 3, size=1)
    p_true = np.concatenate([p_sim, p_dis, p_target])

    # Sample sizes
    n = np.concatenate([
        rng.integers(30, 101, size=H),
        [n_target]
    ]).astype(float)

    y = rng.binomial(n.astype(int), p_true).astype(float)

    # KG similarity: block structure
    W_raw = np.zeros((H_total, H_total))
    # Similar cluster (indices 0..H_sim-1) + target (index H)
    sim_indices = list(range(H_sim)) + [H]
    dis_indices = list(range(H_sim, H))

    for i in sim_indices:
        for j in sim_indices:
            if i != j:
                W_raw[i, j] = rng.uniform(0.7, 0.9)
    for i in dis_indices:
        for j in dis_indices:
            if i != j:
                W_raw[i, j] = rng.uniform(0.7, 0.9)
    # Between clusters: low similarity
    for i in sim_indices:
        for j in dis_indices:
            val = rng.uniform(0.05, 0.15)
            W_raw[i, j] = val
            W_raw[j, i] = val

    W_raw = 0.5 * (W_raw + W_raw.T)
    np.fill_diagonal(W_raw, 0.0)

    return p_true, y, n, W_raw


def generate_scenario_3(rng, H=20, n_target=50):
    """
    Misleading KG: high similarity everywhere, but outcomes have two hidden clusters.

    - H/2 trials from Beta(7,3) (mean ~0.7), H/2 from Beta(3,7) (mean ~0.3)
    - Target drawn from Beta(3,7) (low cluster)
    - W_ij all high (~0.8) — KG claims all trials are similar
    - KG misses the bimodal structure (realistic: same disease, different subpopulations)
    - Tests graceful degradation: rho should shrink, epsilon absorbs
    """
    H_total = H + 1
    H_high = H // 2
    H_low = H - H_high

    # Bimodal true rates: KG misses this structure
    p_high = rng.beta(7, 3, size=H_high)   # mean ~0.7
    p_low = rng.beta(3, 7, size=H_low)     # mean ~0.3
    p_target = rng.beta(3, 7, size=1)      # target in low cluster
    p_true = np.concatenate([p_high, p_low, p_target])

    n = np.concatenate([
        rng.integers(30, 101, size=H),
        [n_target]
    ]).astype(float)

    y = rng.binomial(n.astype(int), p_true).astype(float)

    # KG says all similar (misleading)
    W_raw = rng.uniform(0.7, 0.9, size=(H_total, H_total))
    W_raw = 0.5 * (W_raw + W_raw.T)
    np.fill_diagonal(W_raw, 0.0)

    return p_true, y, n, W_raw


def generate_scenario_4(rng, H_base=20, H_foreign=5, n_target=50):
    """
    Contamination: base trials + foreign trials with very different rates.

    - 20 base trials: p ~ Beta(3,7) (mean ~0.3)
    - 5 foreign trials: p ~ Beta(9,1) (mean ~0.9, opposite direction)
    - Target in base cluster (p ~ Beta(3,7))
    - W low between base and foreign (~0.05), high within base (~0.8)
    """
    H_total = H_base + H_foreign + 1

    p_base = rng.beta(3, 7, size=H_base)
    p_foreign = rng.beta(9, 1, size=H_foreign)
    p_target = rng.beta(3, 7, size=1)
    p_true = np.concatenate([p_base, p_foreign, p_target])

    n = np.concatenate([
        rng.integers(30, 101, size=H_base + H_foreign),
        [n_target]
    ]).astype(float)

    y = rng.binomial(n.astype(int), p_true).astype(float)

    # KG: base cluster high, foreign low, target with base
    W_raw = np.zeros((H_total, H_total))
    base_idx = list(range(H_base)) + [H_base + H_foreign]  # base + target
    foreign_idx = list(range(H_base, H_base + H_foreign))

    for i in base_idx:
        for j in base_idx:
            if i != j:
                W_raw[i, j] = rng.uniform(0.7, 0.9)
    for i in foreign_idx:
        for j in foreign_idx:
            if i != j:
                W_raw[i, j] = rng.uniform(0.3, 0.5)
    for i in base_idx:
        for j in foreign_idx:
            val = rng.uniform(0.02, 0.08)
            W_raw[i, j] = val
            W_raw[j, i] = val

    W_raw = 0.5 * (W_raw + W_raw.T)
    np.fill_diagonal(W_raw, 0.0)

    return p_true, y, n, W_raw


def generate_scenario_5(rng, H=20, n_target=50):
    """
    Corrupted KG: same DGP as Scenario 1 (homogeneous), but W is randomly permuted.

    - True rates homogeneous as in Scenario 1
    - W constructed correctly, then rows/columns randomly permuted
    - Tests robustness: rho should shrink toward 0
    """
    p_true, y, n, W_good = generate_scenario_1(rng, H, n_target)

    H_total = H + 1
    # Randomly permute rows and columns
    perm = rng.permutation(H_total)
    W_raw = W_good[np.ix_(perm, perm)]

    return p_true, y, n, W_raw


SCENARIO_GENERATORS = {
    1: ('Homogeneous', generate_scenario_1),
    2: ('Heterogeneous (clustered)', generate_scenario_2),
    3: ('Misleading KG', generate_scenario_3),
    4: ('Contamination', generate_scenario_4),
    5: ('Corrupted KG', generate_scenario_5),
}


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate_prediction(result, p_true_target):
    """Compute metrics for a single prediction."""
    p_pred = result['p_mean']
    p_samples = result['p_samples']
    ci = result['ci']

    bias = p_pred - p_true_target
    se = bias ** 2
    ae = abs(bias)
    covered = 1 if ci[0] <= p_true_target <= ci[1] else 0
    crps_val = crps_sample(p_true_target, p_samples)
    ci_width = ci[1] - ci[0]

    return {
        'bias': bias,
        'se': se,
        'ae': ae,
        'covered': covered,
        'crps': crps_val,
        'ci_width': ci_width,
        'p_pred': p_pred,
        'p_lo': ci[0],
        'p_hi': ci[1],
        'theta_mean': result['theta_mean'],
        'theta_sd': result['theta_sd'],
    }


def run_single_replicate(scenario_id, rep, rng_data, B, seed_base):
    """
    Run one replicate of a scenario.

    Returns a list of dicts (one per method).
    """
    # Generate data
    gen_name, gen_func = SCENARIO_GENERATORS[scenario_id]
    p_true, y, n_arr, W_raw = gen_func(rng_data)

    H_total = len(p_true)
    target_idx = H_total - 1
    p_true_target = p_true[target_idx]

    # Build adjacency and Laplacian
    W, L = build_adjacency_and_laplacian(W_raw)

    # Set up methods
    methods = {
        'No Borrowing': NoBorrowing(),
        'Full Pooling': FullPooling(),
        'BHM': StandardBHM(),
        'rMAP': RobustMAP(w_rob=0.10),
        'KG-PP': KGPP(W),
        'KG-CAR': KGCAR(W, L, w_rob=0.10),
    }

    records = []
    for method_name, method in methods.items():
        # Deterministic method offset (hash() is non-deterministic across Python sessions)
        method_offset = sum(ord(c) for c in method_name) % 97
        arm_rng = np.random.default_rng(seed_base + rep * 100 + method_offset)
        try:
            result = method.predict_loo_samples(
                target_idx, y, n_arr, W=W, B=B, rng=arm_rng
            )
            metrics = evaluate_prediction(result, p_true_target)

            # Extra info for KG-CAR
            rho_hat = result.get('rho_hat', np.nan)
            sigma_phi = result.get('sigma_phi_hat', np.nan)
            sigma_eps = result.get('sigma_eps_hat', np.nan)

        except Exception as e:
            metrics = {k: np.nan for k in
                       ['bias', 'se', 'ae', 'covered', 'crps', 'ci_width',
                        'p_pred', 'p_lo', 'p_hi', 'theta_mean', 'theta_sd']}
            rho_hat = sigma_phi = sigma_eps = np.nan

        records.append({
            'scenario': scenario_id,
            'scenario_name': gen_name,
            'replicate': rep,
            'method': method_name,
            'p_true': p_true_target,
            'rho_hat': rho_hat,
            'sigma_phi': sigma_phi,
            'sigma_eps': sigma_eps,
            **metrics,
        })

    return records


def run_scenario(scenario_id, R, B, seed=42, verbose=True):
    """Run all replicates for one scenario."""
    gen_name = SCENARIO_GENERATORS[scenario_id][0]
    if verbose:
        print(f"\n{'='*60}")
        print(f"Scenario {scenario_id}: {gen_name}  (R={R}, B={B})")
        print(f"{'='*60}")

    all_records = []
    t0 = time.time()

    for rep in range(R):
        rng_data = np.random.default_rng(seed + rep)
        records = run_single_replicate(scenario_id, rep, rng_data, B, seed)
        all_records.extend(records)

        if verbose and (rep + 1) % max(1, R // 10) == 0:
            elapsed = time.time() - t0
            rate = (rep + 1) / elapsed
            eta = (R - rep - 1) / rate
            print(f"  Rep {rep+1}/{R}  ({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

    results = pd.DataFrame(all_records)

    # Summary
    summary = results.groupby('method').agg(
        Bias=('bias', 'mean'),
        MSE=('se', 'mean'),
        RMSE=('se', lambda x: np.sqrt(x.mean())),
        MAE=('ae', 'mean'),
        Coverage=('covered', 'mean'),
        CRPS=('crps', 'mean'),
        CI_Width=('ci_width', 'mean'),
        # MC standard errors
        MSE_SE=('se', lambda x: x.std() / np.sqrt(len(x))),
        Coverage_SE=('covered', lambda x: x.std() / np.sqrt(len(x))),
        CRPS_SE=('crps', lambda x: x.std() / np.sqrt(len(x))),
    ).reset_index()

    if verbose:
        elapsed = time.time() - t0
        print(f"\n  Completed in {elapsed:.1f}s")
        print(f"\n  Summary for Scenario {scenario_id} ({gen_name}):")
        # Format with SEs inline
        fmt_rows = []
        for _, row in summary.iterrows():
            fmt_rows.append(
                f"  {row['method']:>14s}  "
                f"Bias={row['Bias']:+.4f}  "
                f"MSE={row['MSE']:.5f}({row['MSE_SE']:.5f})  "
                f"Cov={row['Coverage']:.3f}({row['Coverage_SE']:.3f})  "
                f"CRPS={row['CRPS']:.5f}({row['CRPS_SE']:.5f})  "
                f"CIW={row['CI_Width']:.3f}"
            )
        print('\n'.join(fmt_rows))

        # KG-CAR diagnostics
        kgcar = results[results['method'] == 'KG-CAR']
        print(f"\n  KG-CAR diagnostics:")
        print(f"    rho_hat: mean={kgcar['rho_hat'].mean():.3f}, "
              f"sd={kgcar['rho_hat'].std():.3f}")
        print(f"    sigma_phi: mean={kgcar['sigma_phi'].mean():.3f}")
        print(f"    sigma_eps: mean={kgcar['sigma_eps'].mean():.3f}")

        # Relative improvement vs rMAP
        kgcar_mse = summary.loc[summary['method'] == 'KG-CAR', 'MSE'].values[0]
        rmap_mse = summary.loc[summary['method'] == 'rMAP', 'MSE'].values[0]
        pct = (1 - kgcar_mse / rmap_mse) * 100
        print(f"\n  KG-CAR vs rMAP MSE reduction: {pct:+.1f}%")

    return results, summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="KG-CAR Simulation Study")
    parser.add_argument('--full', action='store_true',
                        help='Full run: R=500, B=1000')
    parser.add_argument('--scenario', type=int, default=None,
                        help='Run single scenario (1-5)')
    parser.add_argument('--R', type=int, default=None,
                        help='Override number of replicates')
    parser.add_argument('--B', type=int, default=None,
                        help='Override number of MC samples')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    if args.full:
        R = args.R or 1000
        B = args.B or 2000
    else:
        R = args.R or 20
        B = args.B or 500

    DATA_OUT_DIR.mkdir(parents=True, exist_ok=True)

    scenarios = [args.scenario] if args.scenario else [1, 2, 3, 4, 5]

    print("=" * 60)
    print("KG-CAR SIMULATION STUDY")
    print(f"  R={R} replicates, B={B} MC samples, seed={args.seed}")
    print(f"  Scenarios: {scenarios}")
    print("=" * 60)

    all_results = []
    all_summaries = []
    t_total = time.time()

    for sc in scenarios:
        results, summary = run_scenario(sc, R, B, seed=args.seed)
        results.to_csv(DATA_OUT_DIR / f'simulation_scenario_{sc}.csv', index=False)
        summary.to_csv(DATA_OUT_DIR / f'simulation_summary_{sc}.csv', index=False)
        all_results.append(results)
        all_summaries.append(summary)

    # Combined results
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(DATA_OUT_DIR / 'simulation_all.csv', index=False)

    combined_summary = pd.concat(all_summaries, ignore_index=True)

    # Print grand summary table
    print("\n" + "=" * 60)
    print("GRAND SUMMARY")
    print("=" * 60)

    for sc in scenarios:
        sc_sum = combined_summary[
            combined_summary.index.isin(
                all_summaries[scenarios.index(sc)].index
            )
        ] if len(scenarios) > 1 else combined_summary

    # Pivot: scenario x method for MSE, Coverage, CRPS
    pivot_mse = combined.groupby(['scenario', 'method'])['se'].mean().unstack('method')
    pivot_cov = combined.groupby(['scenario', 'method'])['covered'].mean().unstack('method')
    pivot_crps = combined.groupby(['scenario', 'method'])['crps'].mean().unstack('method')

    # MC standard errors: SE = sd / sqrt(R)
    pivot_mse_se = combined.groupby(['scenario', 'method'])['se'].apply(
        lambda x: x.std() / np.sqrt(len(x))).unstack('method')
    pivot_cov_se = combined.groupby(['scenario', 'method'])['covered'].apply(
        lambda x: x.std() / np.sqrt(len(x))).unstack('method')
    pivot_crps_se = combined.groupby(['scenario', 'method'])['crps'].apply(
        lambda x: x.std() / np.sqrt(len(x))).unstack('method')

    method_order = ['No Borrowing', 'Full Pooling', 'BHM', 'rMAP', 'KG-PP', 'KG-CAR']
    existing_methods = [m for m in method_order if m in pivot_mse.columns]

    # Format: mean (SE)
    def format_with_se(mean_df, se_df, fmt='.5f'):
        result = mean_df.copy().astype(str)
        for col in result.columns:
            for idx in result.index:
                m = mean_df.loc[idx, col]
                s = se_df.loc[idx, col]
                result.loc[idx, col] = f"{m:{fmt}} ({s:{fmt}})"
        return result

    print("\n  MSE (SE) by scenario and method:")
    print(format_with_se(pivot_mse[existing_methods], pivot_mse_se[existing_methods]).to_string())

    print("\n  Coverage (SE) by scenario and method:")
    print(format_with_se(pivot_cov[existing_methods], pivot_cov_se[existing_methods], '.3f').to_string())

    print("\n  CRPS (SE) by scenario and method:")
    print(format_with_se(pivot_crps[existing_methods], pivot_crps_se[existing_methods]).to_string())

    # KG-CAR vs rMAP reduction by scenario
    print("\n  KG-CAR vs rMAP MSE reduction (%) by scenario:")
    for sc in scenarios:
        sc_name = SCENARIO_GENERATORS[sc][0]
        kgcar_mse = pivot_mse.loc[sc, 'KG-CAR']
        rmap_mse = pivot_mse.loc[sc, 'rMAP']
        kgcar_se = pivot_mse_se.loc[sc, 'KG-CAR']
        rmap_se = pivot_mse_se.loc[sc, 'rMAP']
        pct = (1 - kgcar_mse / rmap_mse) * 100
        # SE of difference via delta method (approximate)
        diff = rmap_mse - kgcar_mse
        diff_se = np.sqrt(kgcar_se**2 + rmap_se**2)
        print(f"    Scenario {sc} ({sc_name}): {pct:+.1f}%  "
              f"(MSE diff = {diff:.5f} +/- {diff_se:.5f})")

    elapsed_total = time.time() - t_total
    print(f"\n  Total time: {elapsed_total:.0f}s ({elapsed_total/60:.1f} min)")
    print(f"  Results saved to: {DATA_OUT_DIR}/")


if __name__ == "__main__":
    main()
