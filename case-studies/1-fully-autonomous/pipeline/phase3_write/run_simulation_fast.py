"""
Fast Simulation Study with Known Ground Truth for KG-DAP
=========================================================
Optimized version: precompute similarity matrix once, precompute all
weights once, then only vary the synthetic data.
"""

import os
import sys
import time
import pickle
import numpy as np
import pandas as pd

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
    trials_csv = os.path.join(PROJECT_ROOT, "data", "trials_data.csv")
    trials_df = pd.read_csv(trials_csv)
    graph_path = os.path.join(PROJECT_ROOT, "data", "kg_graph.gpickle")
    with open(graph_path, "rb") as f:
        G = pickle.load(f)
    trial_ids = build_kg.get_trial_ids(G)
    return trials_df, G, trial_ids


def compute_metrics_fast(prior, true_rate, n_obs, y_obs):
    """Compute metrics without interval computation for speed."""
    pm = mixture_mean(prior)
    mae = abs(pm - true_rate)
    # Skip interval for speed during MC reps — compute on summary only
    lower, upper = mixture_interval(prior, level=0.95)
    coverage = 1 if (lower <= true_rate <= upper) else 0
    width = upper - lower
    ess = compute_ess(prior)
    return mae, coverage, width, ess


def main():
    print("=" * 60)
    print("Fast Simulation Study with Known Ground Truth")
    print("=" * 60)

    trials_df, G, trial_ids = load_real_infrastructure()
    K = len(trial_ids)
    B = 200  # Reduced for speed; still gives SE/effect < 50%

    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']),
                                         float(row['p_hat']))

    real_n = {tid: tid_to_data[tid][0] for tid in trial_ids}
    real_rates = {tid: tid_to_data[tid][2] for tid in trial_ids}

    # Precompute similarity matrix ONCE (doesn't change across reps)
    S_full = compute_similarity_matrix(G, trial_ids, alpha=0.20, beta_w=0.20, gamma=0.60)

    # Precompute KG-DAP weights for each held-out arm ONCE
    precomputed_omega = {}
    for h in range(K):
        precomputed_omega[h] = compute_power_weights(S_full, h, beta_power=15.0)

    # Define scenarios
    rng_shuffle = np.random.default_rng(123)
    shuffled = list(real_rates.values())
    rng_shuffle.shuffle(shuffled)

    rng_mix = np.random.default_rng(456)
    mixed_rates = {}
    for i, tid in enumerate(trial_ids):
        if i % 2 == 0:
            mixed_rates[tid] = real_rates[tid]
        else:
            pert = rng_mix.choice([-0.15, 0.15])
            mixed_rates[tid] = np.clip(real_rates[tid] + pert, 0.02, 0.98)

    scenarios = {
        'Favorable': {tid: real_rates[tid] for tid in trial_ids},
        'Adverse': {tid: r for tid, r in zip(trial_ids, shuffled)},
        'Mixed': mixed_rates,
    }

    all_summary = []

    for scen_name, true_rates in scenarios.items():
        print(f"\n[{scen_name}] Running B={B} replications...")
        t0 = time.time()

        # Collect per-rep, per-method MAEs
        method_rep_mae = {m: [] for m in ['KG-DAP', 'rMAP', 'Uniform', 'Pooled', 'EqualWeight']}
        method_rep_cov = {m: [] for m in ['KG-DAP', 'rMAP', 'Uniform', 'Pooled', 'EqualWeight']}
        method_rep_width = {m: [] for m in ['KG-DAP', 'rMAP', 'Uniform', 'Pooled', 'EqualWeight']}
        method_rep_ess = {m: [] for m in ['KG-DAP', 'rMAP', 'Uniform', 'Pooled', 'EqualWeight']}

        for b in range(B):
            rng = np.random.default_rng(42 + b * 1000 + hash(scen_name) % 10000)

            # Generate synthetic data for ALL arms
            syn_y = {}
            for tid in trial_ids:
                syn_y[tid] = rng.binomial(real_n[tid], true_rates[tid])

            # Per-arm metrics
            rep_mae = {m: [] for m in method_rep_mae}
            rep_cov = {m: [] for m in method_rep_mae}
            rep_width = {m: [] for m in method_rep_mae}
            rep_ess = {m: [] for m in method_rep_mae}

            for h in range(K):
                held_id = trial_ids[h]
                n_h = real_n[held_id]
                y_h = syn_y[held_id]
                true_r = true_rates[held_id]

                hist_idx = [i for i in range(K) if i != h]
                n_hist = np.array([real_n[trial_ids[i]] for i in hist_idx], dtype=float)
                y_hist = np.array([syn_y[trial_ids[i]] for i in hist_idx], dtype=float)

                # KG-DAP (precomputed weights)
                prior_dap = construct_prior(precomputed_omega[h], n_hist, y_hist,
                                            w0=0.20, n_cap=200)
                mae, cov, w, e = compute_metrics_fast(prior_dap, true_r, n_h, y_h)
                rep_mae['KG-DAP'].append(mae)
                rep_cov['KG-DAP'].append(cov)
                rep_width['KG-DAP'].append(w)
                rep_ess['KG-DAP'].append(e)

                # rMAP
                prior_rmap = rmap_prior(n_hist, y_hist, w_rob=0.20)
                mae, cov, w, e = compute_metrics_fast(prior_rmap, true_r, n_h, y_h)
                rep_mae['rMAP'].append(mae)
                rep_cov['rMAP'].append(cov)
                rep_width['rMAP'].append(w)
                rep_ess['rMAP'].append(e)

                # Uniform
                prior_unif = uniform_prior()
                mae, cov, w, e = compute_metrics_fast(prior_unif, true_r, n_h, y_h)
                rep_mae['Uniform'].append(mae)
                rep_cov['Uniform'].append(cov)
                rep_width['Uniform'].append(w)
                rep_ess['Uniform'].append(e)

                # Pooled
                prior_pool = pooled_prior(n_hist, y_hist)
                mae, cov, w, e = compute_metrics_fast(prior_pool, true_r, n_h, y_h)
                rep_mae['Pooled'].append(mae)
                rep_cov['Pooled'].append(cov)
                rep_width['Pooled'].append(w)
                rep_ess['Pooled'].append(e)

                # EqualWeight
                prior_ew = equal_weight_prior(n_hist, y_hist, w0=0.20, n_cap=200)
                mae, cov, w, e = compute_metrics_fast(prior_ew, true_r, n_h, y_h)
                rep_mae['EqualWeight'].append(mae)
                rep_cov['EqualWeight'].append(cov)
                rep_width['EqualWeight'].append(w)
                rep_ess['EqualWeight'].append(e)

            # Store per-rep averages
            for m in method_rep_mae:
                method_rep_mae[m].append(np.mean(rep_mae[m]))
                method_rep_cov[m].append(np.mean(rep_cov[m]))
                method_rep_width[m].append(np.mean(rep_width[m]))
                method_rep_ess[m].append(np.mean(rep_ess[m]))

            if (b + 1) % max(1, B // 5) == 0:
                print(f"    Rep {b+1}/{B}")

        elapsed = time.time() - t0
        print(f"    Done in {elapsed:.1f}s")

        # Compute summary
        for m in method_rep_mae:
            arr_mae = np.array(method_rep_mae[m])
            arr_cov = np.array(method_rep_cov[m])
            arr_width = np.array(method_rep_width[m])
            arr_ess = np.array(method_rep_ess[m])

            all_summary.append({
                'scenario': scen_name,
                'method': m,
                'MAE': arr_mae.mean(),
                'MAE_SE': arr_mae.std() / np.sqrt(B),
                'RMSE': np.sqrt((arr_mae**2).mean()),
                'Coverage_95': arr_cov.mean(),
                'Coverage_SE': arr_cov.std() / np.sqrt(B),
                'Mean_Width': arr_width.mean(),
                'Mean_ESS': arr_ess.mean(),
                'N_reps': B,
            })

    summary_df = pd.DataFrame(all_summary)
    summary_df.to_csv(os.path.join(DATA_DIR, "simulation_summary.csv"), index=False)

    print("\n" + "=" * 60)
    print("SIMULATION STUDY SUMMARY")
    print("=" * 60)
    for scenario in ['Favorable', 'Adverse', 'Mixed']:
        print(f"\n--- {scenario} ---")
        sub = summary_df[summary_df['scenario'] == scenario]
        for _, row in sub.iterrows():
            print(f"  {row['method']:15s}  MAE={row['MAE']:.4f} (SE={row['MAE_SE']:.4f})  "
                  f"Cov={row['Coverage_95']:.3f}  Width={row['Mean_Width']:.3f}  "
                  f"ESS={row['Mean_ESS']:.1f}")

    # MC-SE check
    print("\n--- MC-SE / Effect Size Check ---")
    for scenario in ['Favorable', 'Adverse', 'Mixed']:
        sub = summary_df[summary_df['scenario'] == scenario]
        dap = sub[sub['method'] == 'KG-DAP'].iloc[0]
        rmap = sub[sub['method'] == 'rMAP'].iloc[0]
        effect = abs(dap['MAE'] - rmap['MAE'])
        ratio = dap['MAE_SE'] / effect if effect > 0 else float('inf')
        flag = "OK" if ratio < 0.50 else "WARNING"
        print(f"  {scenario}: effect={effect:.4f}, SE={dap['MAE_SE']:.4f}, "
              f"SE/effect={ratio:.2f} [{flag}]")

    print(f"\nSaved to: {os.path.join(DATA_DIR, 'simulation_summary.csv')}")


if __name__ == "__main__":
    main()
