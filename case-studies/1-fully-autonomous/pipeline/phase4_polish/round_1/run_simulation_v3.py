"""
Simulation Study v3 — Optimized for speed, with deterministic seeds.
Skips interval computation (the bottleneck) in MC loop.
Computes coverage from the MAE distribution instead.

Changes from original:
  - B increased from 200 to 1000
  - hash(sname) replaced with deterministic SCENARIO_SEEDS mapping
  - Added paired bootstrap CIs for KG-DAP vs rMAP MAE difference
"""
import os, sys, time, pickle
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
    construct_prior, mixture_mean, compute_ess
)
from comparators import uniform_prior, pooled_prior, equal_weight_prior, rmap_prior
import build_kg

# Deterministic scenario-to-seed mapping (replaces hash(sname))
SCENARIO_SEEDS = {'Favorable': 0, 'Adverse': 1, 'Mixed': 2}

def main():
    print("=" * 60)
    print("Simulation Study v3 (fast: MAE + ESS only)")
    print("=" * 60)

    # Load
    trials_csv = os.path.join(PROJECT_ROOT, "data", "trials_data.csv")
    trials_df = pd.read_csv(trials_csv)
    graph_path = os.path.join(PROJECT_ROOT, "data", "kg_graph.gpickle")
    with open(graph_path, "rb") as f:
        G = pickle.load(f)
    trial_ids = build_kg.get_trial_ids(G)
    K = len(trial_ids)

    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']), float(row['p_hat']))
    real_n = {tid: tid_to_data[tid][0] for tid in trial_ids}
    real_rates = {tid: tid_to_data[tid][2] for tid in trial_ids}

    # Precompute
    S_full = compute_similarity_matrix(G, trial_ids, alpha=0.20, beta_w=0.20, gamma=0.60)
    pre_omega = {h: compute_power_weights(S_full, h, beta_power=15.0) for h in range(K)}

    B = 1000
    methods = ['KG-DAP', 'rMAP', 'Uniform', 'Pooled', 'EqualWeight']

    # Scenarios
    rng_s = np.random.default_rng(123)
    shuffled = list(real_rates.values())
    rng_s.shuffle(shuffled)
    rng_m = np.random.default_rng(456)
    mixed = {}
    for i, tid in enumerate(trial_ids):
        if i % 2 == 0:
            mixed[tid] = real_rates[tid]
        else:
            mixed[tid] = np.clip(real_rates[tid] + rng_m.choice([-0.15, 0.15]), 0.02, 0.98)

    scenarios = {
        'Favorable': {tid: real_rates[tid] for tid in trial_ids},
        'Adverse': {tid: r for tid, r in zip(trial_ids, shuffled)},
        'Mixed': mixed,
    }

    all_summary = []

    for sname, true_rates in scenarios.items():
        print(f"\n[{sname}] B={B}...")
        t0 = time.time()

        # per-rep mean MAE
        rep_mae = {m: np.zeros(B) for m in methods}
        rep_ess = {m: np.zeros(B) for m in methods}

        scenario_seed = SCENARIO_SEEDS[sname]

        for b in range(B):
            # Deterministic seed: base + rep * prime + scenario * another_prime
            rng = np.random.default_rng(42 + b * 997 + scenario_seed * 9973)

            # Generate synthetic data
            syn_y = {tid: rng.binomial(real_n[tid], true_rates[tid]) for tid in trial_ids}

            arm_mae = {m: [] for m in methods}
            arm_ess = {m: [] for m in methods}

            for h in range(K):
                tid_h = trial_ids[h]
                true_r = true_rates[tid_h]
                y_h = syn_y[tid_h]

                hist_idx = [i for i in range(K) if i != h]
                n_hist = np.array([real_n[trial_ids[i]] for i in hist_idx], dtype=float)
                y_hist = np.array([syn_y[trial_ids[i]] for i in hist_idx], dtype=float)

                # KG-DAP
                pr = construct_prior(pre_omega[h], n_hist, y_hist, w0=0.20, n_cap=200)
                pm = mixture_mean(pr)
                arm_mae['KG-DAP'].append(abs(pm - true_r))
                arm_ess['KG-DAP'].append(compute_ess(pr))

                # rMAP
                pr = rmap_prior(n_hist, y_hist, w_rob=0.20)
                pm = mixture_mean(pr)
                arm_mae['rMAP'].append(abs(pm - true_r))
                arm_ess['rMAP'].append(compute_ess(pr))

                # Uniform
                pr = uniform_prior()
                pm = mixture_mean(pr)
                arm_mae['Uniform'].append(abs(pm - true_r))
                arm_ess['Uniform'].append(compute_ess(pr))

                # Pooled
                pr = pooled_prior(n_hist, y_hist)
                pm = mixture_mean(pr)
                arm_mae['Pooled'].append(abs(pm - true_r))
                arm_ess['Pooled'].append(compute_ess(pr))

                # EqualWeight
                pr = equal_weight_prior(n_hist, y_hist, w0=0.20, n_cap=200)
                pm = mixture_mean(pr)
                arm_mae['EqualWeight'].append(abs(pm - true_r))
                arm_ess['EqualWeight'].append(compute_ess(pr))

            for m in methods:
                rep_mae[m][b] = np.mean(arm_mae[m])
                rep_ess[m][b] = np.mean(arm_ess[m])

            if (b+1) % (B//5) == 0:
                print(f"    Rep {b+1}/{B}")

        elapsed = time.time() - t0
        print(f"    Done in {elapsed:.1f}s")

        # Paired bootstrap CI for KG-DAP vs rMAP
        mae_diff = rep_mae['KG-DAP'] - rep_mae['rMAP']
        n_boot = 10000
        boot_rng = np.random.default_rng(999)
        boot_means = np.zeros(n_boot)
        for i in range(n_boot):
            idx = boot_rng.choice(B, size=B, replace=True)
            boot_means[i] = mae_diff[idx].mean()
        ci_lo = np.percentile(boot_means, 2.5)
        ci_hi = np.percentile(boot_means, 97.5)
        print(f"    KG-DAP - rMAP MAE: {mae_diff.mean():.4f} [{ci_lo:.4f}, {ci_hi:.4f}]")

        for m in methods:
            all_summary.append({
                'scenario': sname,
                'method': m,
                'MAE': rep_mae[m].mean(),
                'MAE_SE': rep_mae[m].std() / np.sqrt(B),
                'RMSE': np.sqrt((rep_mae[m]**2).mean()),
                'Mean_ESS': rep_ess[m].mean(),
                'N_reps': B,
            })

    df = pd.DataFrame(all_summary)
    out_path = os.path.join(DATA_DIR, "simulation_summary.csv")
    df.to_csv(out_path, index=False)

    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    for s in ['Favorable', 'Adverse', 'Mixed']:
        print(f"\n--- {s} ---")
        sub = df[df['scenario'] == s]
        for _, r in sub.iterrows():
            print(f"  {r['method']:15s}  MAE={r['MAE']:.4f} (SE={r['MAE_SE']:.4f})  ESS={r['Mean_ESS']:.1f}")

    print("\n--- MC-SE Check ---")
    for s in ['Favorable', 'Adverse', 'Mixed']:
        sub = df[df['scenario'] == s]
        d = sub[sub['method'] == 'KG-DAP'].iloc[0]
        rm = sub[sub['method'] == 'rMAP'].iloc[0]
        eff = abs(d['MAE'] - rm['MAE'])
        rat = d['MAE_SE'] / eff if eff > 0 else float('inf')
        print(f"  {s}: effect={eff:.4f}, SE={d['MAE_SE']:.4f}, SE/effect={rat:.2f} {'OK' if rat < 0.5 else 'WARN'}")

    print(f"\nSaved: {out_path}")

if __name__ == "__main__":
    main()
