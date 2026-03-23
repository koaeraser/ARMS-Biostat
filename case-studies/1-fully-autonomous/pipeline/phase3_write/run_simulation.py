"""
Simulation Study with Known Ground Truth for KG-DAP
=====================================================
Implements 3 diverse scenarios with known parameters to assess
KG-DAP under controlled conditions.

Uses the validated code as foundation (import, don't rewrite).
"""

import os
import sys
import time
import pickle
import numpy as np
import pandas as pd

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
VALIDATED_CODE = os.path.join(PROJECT_ROOT, "pipeline", "phase2_validate", "validated_code")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")

sys.path.insert(0, VALIDATED_CODE)
sys.path.insert(0, PROJECT_ROOT)

from kg_dap import (
    compute_similarity_matrix, compute_power_weights,
    construct_prior, mixture_mean, mixture_variance,
    mixture_interval, compute_ess, log_predictive_pmf
)
from comparators import uniform_prior, pooled_prior, equal_weight_prior, rmap_prior
import build_kg


def load_real_infrastructure():
    """Load the real KG and trial data for constructing realistic scenarios."""
    trials_csv = os.path.join(PROJECT_ROOT, "data", "trials_data.csv")
    trials_df = pd.read_csv(trials_csv)

    graph_path = os.path.join(PROJECT_ROOT, "data", "kg_graph.gpickle")
    with open(graph_path, "rb") as f:
        G = pickle.load(f)

    trial_ids = build_kg.get_trial_ids(G)
    return trials_df, G, trial_ids


def compute_metrics(prior, true_rate, n_obs, y_obs):
    """Compute evaluation metrics for a prior against known ground truth."""
    pm = mixture_mean(prior)
    mae = abs(pm - true_rate)
    lower, upper = mixture_interval(prior, level=0.95)
    coverage = 1 if (lower <= true_rate <= upper) else 0
    width = upper - lower
    log_pred = log_predictive_pmf(prior, y_obs, n_obs)
    ess = compute_ess(prior)
    return {
        'prior_mean': pm,
        'abs_error': mae,
        'coverage_95': coverage,
        'interval_lower': lower,
        'interval_upper': upper,
        'interval_width': width,
        'log_pred_score': log_pred,
        'ess': ess,
    }


def run_scenario(scenario_name, true_rates, n_arms, trial_ids, G, trials_df,
                 B=2000, seed=42):
    """Run one simulation scenario with known ground truth.

    Parameters
    ----------
    scenario_name : str
    true_rates : dict mapping trial_id -> true_rate
    n_arms : dict mapping trial_id -> sample size
    trial_ids : list of all trial IDs
    G : knowledge graph
    trials_df : DataFrame with original trial data (for comparators)
    B : int, number of Monte Carlo replications
    seed : int

    Returns
    -------
    results : list of dicts, one per (replication, arm, method)
    """
    rng = np.random.default_rng(seed)
    K = len(trial_ids)

    # Precompute similarity matrix once
    S_full = compute_similarity_matrix(G, trial_ids, alpha=0.20, beta_w=0.20, gamma=0.60)

    all_results = []

    for b in range(B):
        # Generate synthetic data from known true rates
        synthetic_data = {}
        for tid in trial_ids:
            n_i = n_arms[tid]
            p_i = true_rates[tid]
            y_i = rng.binomial(n_i, p_i)
            synthetic_data[tid] = (n_i, y_i)

        # For each arm, hold it out and build priors
        for h_idx in range(K):
            held_out_id = trial_ids[h_idx]
            true_rate = true_rates[held_out_id]
            n_h = synthetic_data[held_out_id][0]
            y_h = synthetic_data[held_out_id][1]

            # Historical data (exclude held-out)
            hist_indices = [i for i in range(K) if i != h_idx]
            n_hist = np.array([synthetic_data[trial_ids[i]][0] for i in hist_indices], dtype=float)
            y_hist = np.array([synthetic_data[trial_ids[i]][1] for i in hist_indices], dtype=float)

            # --- KG-DAP ---
            omega = compute_power_weights(S_full, h_idx, beta_power=15.0)
            prior_dap = construct_prior(omega, n_hist, y_hist, w0=0.20, n_cap=200)
            metrics_dap = compute_metrics(prior_dap, true_rate, n_h, y_h)
            metrics_dap.update({
                'scenario': scenario_name, 'rep': b, 'trial_id': held_out_id,
                'true_rate': true_rate, 'n': n_h, 'y': y_h, 'method': 'KG-DAP'
            })
            all_results.append(metrics_dap)

            # --- rMAP ---
            prior_rmap = rmap_prior(n_hist, y_hist, w_rob=0.20)
            metrics_rmap = compute_metrics(prior_rmap, true_rate, n_h, y_h)
            metrics_rmap.update({
                'scenario': scenario_name, 'rep': b, 'trial_id': held_out_id,
                'true_rate': true_rate, 'n': n_h, 'y': y_h, 'method': 'rMAP'
            })
            all_results.append(metrics_rmap)

            # --- Uniform ---
            prior_unif = uniform_prior()
            metrics_unif = compute_metrics(prior_unif, true_rate, n_h, y_h)
            metrics_unif.update({
                'scenario': scenario_name, 'rep': b, 'trial_id': held_out_id,
                'true_rate': true_rate, 'n': n_h, 'y': y_h, 'method': 'Uniform'
            })
            all_results.append(metrics_unif)

            # --- Pooled ---
            prior_pool = pooled_prior(n_hist, y_hist)
            metrics_pool = compute_metrics(prior_pool, true_rate, n_h, y_h)
            metrics_pool.update({
                'scenario': scenario_name, 'rep': b, 'trial_id': held_out_id,
                'true_rate': true_rate, 'n': n_h, 'y': y_h, 'method': 'Pooled'
            })
            all_results.append(metrics_pool)

            # --- EqualWeight ---
            prior_ew = equal_weight_prior(n_hist, y_hist, w0=0.20, n_cap=200)
            metrics_ew = compute_metrics(prior_ew, true_rate, n_h, y_h)
            metrics_ew.update({
                'scenario': scenario_name, 'rep': b, 'trial_id': held_out_id,
                'true_rate': true_rate, 'n': n_h, 'y': y_h, 'method': 'EqualWeight'
            })
            all_results.append(metrics_ew)

        if (b + 1) % max(1, B // 10) == 0:
            print(f"    {scenario_name}: rep {b+1}/{B}")

    return all_results


def main():
    print("=" * 60)
    print("Simulation Study with Known Ground Truth")
    print("=" * 60)

    trials_df, G, trial_ids = load_real_infrastructure()
    K = len(trial_ids)

    # Get real sample sizes for each trial
    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']),
                                         float(row['p_hat']))

    real_n = {tid: tid_to_data[tid][0] for tid in trial_ids}
    real_rates = {tid: tid_to_data[tid][2] for tid in trial_ids}

    B = 500  # Monte Carlo replications (sufficient for SE < 50% of effect)

    # =========================================================
    # Scenario 1: FAVORABLE — True rates match real data
    # KG structure IS informative because real rates correlate
    # with KG similarity
    # =========================================================
    print("\n[Scenario 1] Favorable: true rates = observed rates")
    t0 = time.time()
    true_rates_1 = {tid: tid_to_data[tid][2] for tid in trial_ids}
    results_1 = run_scenario("Favorable", true_rates_1, real_n,
                              trial_ids, G, trials_df, B=B, seed=42)
    t1 = time.time() - t0
    print(f"    Completed in {t1:.1f}s")

    # =========================================================
    # Scenario 2: ADVERSE — True rates are shuffled randomly
    # KG structure is UNINFORMATIVE (rates don't correlate with
    # graph similarity)
    # =========================================================
    print("\n[Scenario 2] Adverse: shuffled true rates")
    rng = np.random.default_rng(123)
    shuffled_rates = list(real_rates.values())
    rng.shuffle(shuffled_rates)
    true_rates_2 = {tid: r for tid, r in zip(trial_ids, shuffled_rates)}
    t0 = time.time()
    results_2 = run_scenario("Adverse", true_rates_2, real_n,
                              trial_ids, G, trials_df, B=B, seed=43)
    t2 = time.time() - t0
    print(f"    Completed in {t2:.1f}s")

    # =========================================================
    # Scenario 3: MIXED — Half the trials keep their real rates,
    # half are perturbed by +/-0.15
    # KG structure is partially informative
    # =========================================================
    print("\n[Scenario 3] Mixed: half true, half perturbed")
    rng3 = np.random.default_rng(456)
    true_rates_3 = {}
    for i, tid in enumerate(trial_ids):
        if i % 2 == 0:
            true_rates_3[tid] = real_rates[tid]
        else:
            perturbation = rng3.choice([-0.15, 0.15])
            true_rates_3[tid] = np.clip(real_rates[tid] + perturbation, 0.02, 0.98)
    t0 = time.time()
    results_3 = run_scenario("Mixed", true_rates_3, real_n,
                              trial_ids, G, trials_df, B=B, seed=44)
    t3 = time.time() - t0
    print(f"    Completed in {t3:.1f}s")

    # =========================================================
    # Aggregate and save results
    # =========================================================
    all_results = results_1 + results_2 + results_3
    df = pd.DataFrame(all_results)

    # Save arm-level results
    df.to_csv(os.path.join(DATA_DIR, "simulation_arm_level.csv"), index=False)

    # Compute summary per scenario x method
    summary_rows = []
    for scenario in ['Favorable', 'Adverse', 'Mixed']:
        for method in ['KG-DAP', 'rMAP', 'Uniform', 'Pooled', 'EqualWeight']:
            sub = df[(df['scenario'] == scenario) & (df['method'] == method)]
            # Average MAE across all (rep x arm) combinations
            # Then compute MC standard error of the mean across reps
            rep_maes = sub.groupby('rep')['abs_error'].mean()
            rep_coverage = sub.groupby('rep')['coverage_95'].mean()
            rep_width = sub.groupby('rep')['interval_width'].mean()

            mae_mean = rep_maes.mean()
            mae_se = rep_maes.std() / np.sqrt(len(rep_maes))
            cov_mean = rep_coverage.mean()
            cov_se = rep_coverage.std() / np.sqrt(len(rep_coverage))
            width_mean = rep_width.mean()

            summary_rows.append({
                'scenario': scenario,
                'method': method,
                'MAE': mae_mean,
                'MAE_SE': mae_se,
                'RMSE': np.sqrt((sub['abs_error']**2).mean()),
                'Coverage_95': cov_mean,
                'Coverage_SE': cov_se,
                'Mean_Width': width_mean,
                'Mean_ESS': sub['ess'].mean(),
                'N_reps': B,
                'N_arms': K,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(DATA_DIR, "simulation_summary.csv"), index=False)

    # Print summary
    print("\n" + "=" * 60)
    print("SIMULATION STUDY SUMMARY")
    print("=" * 60)
    for scenario in ['Favorable', 'Adverse', 'Mixed']:
        print(f"\n--- {scenario} ---")
        sub = summary_df[summary_df['scenario'] == scenario]
        print(sub[['method', 'MAE', 'MAE_SE', 'Coverage_95', 'Mean_Width', 'Mean_ESS']].to_string(index=False))

    # Check MC-SE vs effect size
    print("\n--- MC-SE / Effect Size Check ---")
    for scenario in ['Favorable', 'Adverse', 'Mixed']:
        sub = summary_df[summary_df['scenario'] == scenario]
        dap_mae = sub[sub['method'] == 'KG-DAP']['MAE'].values[0]
        rmap_mae = sub[sub['method'] == 'rMAP']['MAE'].values[0]
        effect = abs(dap_mae - rmap_mae)
        dap_se = sub[sub['method'] == 'KG-DAP']['MAE_SE'].values[0]
        ratio = dap_se / effect if effect > 0 else float('inf')
        flag = "OK" if ratio < 0.50 else "WARNING: SE > 50% of effect"
        print(f"  {scenario}: effect={effect:.4f}, SE={dap_se:.4f}, SE/effect={ratio:.2f} [{flag}]")

    print(f"\nFiles saved to: {DATA_DIR}")


if __name__ == "__main__":
    main()
