"""
KG-DAP Validation Runner
=========================

Runs the full validation protocol:
  1. LOO Cross-Validation (primary evaluation)
  2. Stress Tests: Sensitivity, Contamination, Ablation
  3. Saves all results to CSV and generates figures.

Usage:
    cd pipeline/phase2_validate/validated_code
    python run_validation.py [--quick] [--stress]
"""

import os
import sys
import time
import pickle
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Ensure we can import from this directory and project root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

from kg_dap import (
    build_kg_dap_prior_loocv,
    build_kg_dap_prior_with_external,
    mixture_mean, mixture_variance, mixture_interval,
    compute_ess, log_predictive_pmf, posterior_update,
    compute_diffusion_weights, compute_similarity_matrix,
    construct_prior
)
from comparators import (
    uniform_prior, pooled_prior, equal_weight_prior, rmap_prior
)
import build_kg

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "..", "validated_results")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Default hyperparameters (from methodology specification)
# ---------------------------------------------------------------------------
DEFAULT_PARAMS = {
    'alpha': 0.20,
    'beta_w': 0.20,
    'gamma': 0.60,
    'beta_diff': 15.0,       # Power sharpness (was 2.0 for diffusion)
    'w0': 0.20,
    'alpha0': 1.0,
    'beta0': 1.0,
    'n_cap': 200,
    'weight_method': 'power',  # 'power' or 'diffusion'
}


# ===================================================================
# Data loading
# ===================================================================
def load_data():
    """Load trials data and build knowledge graph."""
    trials_csv = os.path.join(DATA_DIR, "trials_data.csv")
    trials_df = pd.read_csv(trials_csv)

    # Load KG
    graph_path = os.path.join(DATA_DIR, "kg_graph.gpickle")
    if os.path.exists(graph_path):
        with open(graph_path, "rb") as f:
            G = pickle.load(f)
    else:
        # Rebuild
        drug_target = pd.read_csv(os.path.join(DATA_DIR, "drug_target_map.csv"))
        G = build_kg.build_graph(trials_df, drug_target)

    trial_ids = build_kg.get_trial_ids(G)
    return trials_df, G, trial_ids


def load_external_data():
    """Load external adjacent and foreign trial data."""
    adj_path = os.path.join(DATA_DIR, "external_adjacent_trials.csv")
    for_path = os.path.join(DATA_DIR, "external_foreign_trials.csv")
    adj_df = pd.read_csv(adj_path)
    for_df = pd.read_csv(for_path)
    return adj_df, for_df


# ===================================================================
# LOO-CV Evaluation
# ===================================================================
def run_loocv(trials_df, G, trial_ids, params=None, method_name="KG-DAP",
              B=200, seed=42):
    """Run Leave-One-Out Cross-Validation for KG-DAP.

    Parameters
    ----------
    trials_df : pd.DataFrame
    G : nx.MultiDiGraph
    trial_ids : list
    params : dict
        Hyperparameters (defaults used if None).
    method_name : str
    B : int
        Not used for KG-DAP (closed-form), kept for API consistency.
    seed : int
        Random seed.

    Returns
    -------
    results : pd.DataFrame
        One row per held-out arm with columns:
        trial_id, observed_rate, prior_mean, mae_i, coverage_95,
        interval_width, log_pred_score, ess
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()

    K = len(trial_ids)
    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']),
                                         float(row['p_hat']))

    results = []
    for h in range(K):
        held_out_id = trial_ids[h]
        n_h, y_h, p_h = tid_to_data[held_out_id]

        # Build KG-DAP prior holding out arm h
        prior = build_kg_dap_prior_loocv(
            G, trial_ids, trials_df, h,
            alpha=params['alpha'], beta_w=params['beta_w'],
            gamma=params['gamma'], beta_diff=params['beta_diff'],
            w0=params['w0'], alpha0=params['alpha0'],
            beta0=params['beta0'], n_cap=params['n_cap'],
            weight_method=params.get('weight_method', 'power')
        )

        # Metrics
        prior_mean = mixture_mean(prior)
        mae_i = abs(prior_mean - p_h)
        lower, upper = mixture_interval(prior, level=0.95)
        coverage_i = 1 if (lower <= p_h <= upper) else 0
        width_i = upper - lower
        log_pred = log_predictive_pmf(prior, y_h, n_h)
        ess_i = compute_ess(prior)

        results.append({
            'trial_id': held_out_id,
            'observed_rate': p_h,
            'n': n_h,
            'y': y_h,
            'prior_mean': prior_mean,
            'abs_error': mae_i,
            'coverage_95': coverage_i,
            'interval_lower': lower,
            'interval_upper': upper,
            'interval_width': width_i,
            'log_pred_score': log_pred,
            'ess': ess_i,
            'method': method_name,
        })

    return pd.DataFrame(results)


def run_loocv_comparators(trials_df, G, trial_ids, params=None, seed=42):
    """Run LOO-CV for all comparator methods.

    Returns
    -------
    all_results : pd.DataFrame
        Combined results for Uniform, Pooled, EqualWeight, rMAP.
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()

    K = len(trial_ids)
    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']),
                                         float(row['p_hat']))

    all_results = []

    for h in range(K):
        held_out_id = trial_ids[h]
        n_h, y_h, p_h = tid_to_data[held_out_id]

        # Historical data (exclude held-out)
        hist_ids = [trial_ids[i] for i in range(K) if i != h]
        n_hist = np.array([tid_to_data[t][0] for t in hist_ids], dtype=float)
        y_hist = np.array([tid_to_data[t][1] for t in hist_ids], dtype=float)

        # --- Uniform ---
        prior_unif = uniform_prior()
        _append_metrics(all_results, prior_unif, held_out_id, n_h, y_h, p_h,
                       "Uniform")

        # --- Pooled ---
        prior_pool = pooled_prior(n_hist, y_hist)
        _append_metrics(all_results, prior_pool, held_out_id, n_h, y_h, p_h,
                       "Pooled")

        # --- Equal-Weight Mixture ---
        prior_ew = equal_weight_prior(n_hist, y_hist, w0=params['w0'],
                                      n_cap=params['n_cap'])
        _append_metrics(all_results, prior_ew, held_out_id, n_h, y_h, p_h,
                       "EqualWeight")

        # --- rMAP ---
        prior_rmap = rmap_prior(n_hist, y_hist, w_rob=params['w0'])
        _append_metrics(all_results, prior_rmap, held_out_id, n_h, y_h, p_h,
                       "rMAP")

    return pd.DataFrame(all_results)


def _append_metrics(results_list, prior, trial_id, n_h, y_h, p_h, method_name):
    """Compute metrics for a prior and append to results list."""
    prior_mean = mixture_mean(prior)
    mae_i = abs(prior_mean - p_h)
    lower, upper = mixture_interval(prior, level=0.95)
    coverage_i = 1 if (lower <= p_h <= upper) else 0
    width_i = upper - lower
    log_pred = log_predictive_pmf(prior, y_h, n_h)
    ess_i = compute_ess(prior)

    results_list.append({
        'trial_id': trial_id,
        'observed_rate': p_h,
        'n': n_h,
        'y': y_h,
        'prior_mean': prior_mean,
        'abs_error': mae_i,
        'coverage_95': coverage_i,
        'interval_lower': lower,
        'interval_upper': upper,
        'interval_width': width_i,
        'log_pred_score': log_pred,
        'ess': ess_i,
        'method': method_name,
    })


# ===================================================================
# Summary statistics
# ===================================================================
def compute_summary(results_df):
    """Compute summary statistics per method.

    Returns
    -------
    summary : pd.DataFrame
        One row per method with MAE, RMSE, Coverage, MeanWidth, MeanLogPred, MeanESS.
    """
    summary = results_df.groupby('method').agg(
        MAE=('abs_error', 'mean'),
        RMSE=('abs_error', lambda x: np.sqrt(np.mean(x**2))),
        Coverage_95=('coverage_95', 'mean'),
        Mean_Width=('interval_width', 'mean'),
        Mean_Log_Pred=('log_pred_score', 'mean'),
        Mean_ESS=('ess', 'mean'),
        N_arms=('trial_id', 'count'),
    ).reset_index()
    return summary


# ===================================================================
# Stress Test: Sensitivity to beta_diff
# ===================================================================
def sensitivity_beta(trials_df, G, trial_ids, beta_values=None, params=None):
    """Run LOO-CV for a range of diffusion parameter values.

    Returns
    -------
    sens_df : pd.DataFrame
        Summary metrics for each beta value.
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()
    if beta_values is None:
        beta_values = [1.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0]

    sens_results = []
    for bv in beta_values:
        p = params.copy()
        p['beta_diff'] = bv
        res = run_loocv(trials_df, G, trial_ids, params=p,
                       method_name=f"KG-DAP(beta={bv})")
        summary = compute_summary(res)
        summary['beta_diff'] = bv
        sens_results.append(summary)

    return pd.concat(sens_results, ignore_index=True)


# ===================================================================
# Stress Test: Sensitivity to w0
# ===================================================================
def sensitivity_w0(trials_df, G, trial_ids, w0_values=None, params=None):
    """Run LOO-CV for a range of robustness weight values."""
    if params is None:
        params = DEFAULT_PARAMS.copy()
    if w0_values is None:
        w0_values = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]

    sens_results = []
    for wv in w0_values:
        p = params.copy()
        p['w0'] = wv
        res = run_loocv(trials_df, G, trial_ids, params=p,
                       method_name=f"KG-DAP(w0={wv})")
        summary = compute_summary(res)
        summary['w0'] = wv
        sens_results.append(summary)

    return pd.concat(sens_results, ignore_index=True)


# ===================================================================
# Stress Test: Contamination (add external trials)
# ===================================================================
def contamination_test(trials_df, G, trial_ids, adj_df, for_df, params=None):
    """Test robustness by adding adjacent and foreign trials.

    Scenarios:
    1. Clean (baseline, internal only)
    2. +10 adjacent trials
    3. +10 foreign trials
    4. +10 adjacent + 10 foreign
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()

    K = len(trial_ids)
    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']),
                                         float(row['p_hat']))

    scenarios = {
        'Clean': pd.DataFrame(),
        '+10 Adjacent': adj_df.head(10),
        '+10 Foreign': for_df.head(10),
        '+10 Adj +10 For': pd.concat([adj_df.head(10), for_df.head(10)],
                                       ignore_index=True),
    }

    all_results = []

    for scenario_name, ext_df in scenarios.items():
        for h in range(K):
            held_out_id = trial_ids[h]
            n_h, y_h, p_h = tid_to_data[held_out_id]

            if len(ext_df) == 0:
                # Clean — use standard LOO
                prior = build_kg_dap_prior_loocv(
                    G, trial_ids, trials_df, h,
                    alpha=params['alpha'], beta_w=params['beta_w'],
                    gamma=params['gamma'], beta_diff=params['beta_diff'],
                    w0=params['w0'], alpha0=params['alpha0'],
                    beta0=params['beta0'], n_cap=params['n_cap'],
                    weight_method=params.get('weight_method', 'power')
                )
            else:
                # With external trials
                prior = build_kg_dap_prior_with_external(
                    G, trial_ids, trials_df, h, ext_df,
                    alpha=params['alpha'], beta_w=params['beta_w'],
                    gamma=params['gamma'], beta_diff=params['beta_diff'],
                    w0=params['w0'], alpha0=params['alpha0'],
                    beta0=params['beta0'], n_cap=params['n_cap'],
                    weight_method=params.get('weight_method', 'power')
                )

            prior_mean = mixture_mean(prior)
            mae_i = abs(prior_mean - p_h)
            lower, upper = mixture_interval(prior, level=0.95)
            coverage_i = 1 if (lower <= p_h <= upper) else 0

            all_results.append({
                'scenario': scenario_name,
                'trial_id': held_out_id,
                'observed_rate': p_h,
                'prior_mean': prior_mean,
                'abs_error': mae_i,
                'coverage_95': coverage_i,
                'interval_width': upper - lower,
            })

    return pd.DataFrame(all_results)


# ===================================================================
# Stress Test: Ablation
# ===================================================================
def ablation_test(trials_df, G, trial_ids, params=None):
    """Ablation study: remove components one at a time.

    Variants:
    1. Full KG-DAP (baseline)
    2. No graph (equal weights, omega_k=1/K)
    3. No robustness (w0=0)
    4. No sample-size cap (n_cap=inf)
    5. Drug-only similarity (alpha=1, beta_w=0, gamma=0)
    6. Population-only similarity (alpha=0, beta_w=0, gamma=1)
    """
    if params is None:
        params = DEFAULT_PARAMS.copy()

    K = len(trial_ids)
    tid_to_data = {}
    for _, row in trials_df.iterrows():
        tid_to_data[row['trial_id']] = (int(row['n']), int(row['y']),
                                         float(row['p_hat']))

    variants = {
        'Full KG-DAP': params.copy(),
        'No Graph (EqualWeight)': None,  # handled separately
        'No Robustness (w0=0)': {**params, 'w0': 0.001},  # near-zero, avoid division
        'No Cap (n_cap=inf)': {**params, 'n_cap': 1e9},
        'Drug-Only Sim': {**params, 'alpha': 1.0, 'beta_w': 0.0, 'gamma': 0.0},
        'Pop-Only Sim': {**params, 'alpha': 0.0, 'beta_w': 0.0, 'gamma': 1.0},
        'Diffusion Kernel': {**params, 'weight_method': 'diffusion', 'beta_diff': 2.0},
    }

    all_results = []

    for variant_name, variant_params in variants.items():
        for h in range(K):
            held_out_id = trial_ids[h]
            n_h, y_h, p_h = tid_to_data[held_out_id]

            if variant_name == 'No Graph (EqualWeight)':
                # Equal weights, same mixture structure
                hist_ids = [trial_ids[i] for i in range(K) if i != h]
                n_hist = np.array([tid_to_data[t][0] for t in hist_ids], dtype=float)
                y_hist = np.array([tid_to_data[t][1] for t in hist_ids], dtype=float)
                prior = equal_weight_prior(n_hist, y_hist, w0=params['w0'],
                                          n_cap=params['n_cap'])
            else:
                prior = build_kg_dap_prior_loocv(
                    G, trial_ids, trials_df, h,
                    alpha=variant_params['alpha'], beta_w=variant_params['beta_w'],
                    gamma=variant_params['gamma'], beta_diff=variant_params['beta_diff'],
                    w0=variant_params['w0'], alpha0=variant_params['alpha0'],
                    beta0=variant_params['beta0'], n_cap=variant_params['n_cap'],
                    weight_method=variant_params.get('weight_method', 'power')
                )

            prior_mean = mixture_mean(prior)
            mae_i = abs(prior_mean - p_h)
            lower, upper = mixture_interval(prior, level=0.95)
            coverage_i = 1 if (lower <= p_h <= upper) else 0

            all_results.append({
                'variant': variant_name,
                'trial_id': held_out_id,
                'observed_rate': p_h,
                'prior_mean': prior_mean,
                'abs_error': mae_i,
                'coverage_95': coverage_i,
                'interval_width': upper - lower,
            })

    return pd.DataFrame(all_results)


# ===================================================================
# Figures
# ===================================================================
def plot_loocv_comparison(results_df, output_path):
    """Bar chart of MAE by method."""
    summary = compute_summary(results_df)
    methods = summary['method'].values
    mae_vals = summary['MAE'].values

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f']
    bars = ax.bar(range(len(methods)), mae_vals, color=colors[:len(methods)])
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, rotation=30, ha='right')
    ax.set_ylabel('MAE (prior predictive mean vs observed)')
    ax.set_title('LOO-CV: MAE by Method')
    ax.set_ylim(0, max(mae_vals) * 1.3)
    for bar, val in zip(bars, mae_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
               f'{val:.4f}', ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {output_path}")


def plot_sensitivity_beta(sens_df, output_path):
    """Line plot of MAE vs beta_diff."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(sens_df['beta_diff'], sens_df['MAE'], 'o-', color='#4e79a7',
           linewidth=2, markersize=7)
    ax.set_xlabel('Diffusion parameter (beta)')
    ax.set_ylabel('MAE')
    ax.set_title('Sensitivity of KG-DAP MAE to Diffusion Parameter')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {output_path}")


def plot_sensitivity_w0(sens_df, output_path):
    """Line plot of MAE and Coverage vs w0."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(sens_df['w0'], sens_df['MAE'], 'o-', color='#4e79a7',
            linewidth=2, markersize=7)
    ax1.set_xlabel('Robustness weight (w0)')
    ax1.set_ylabel('MAE')
    ax1.set_title('MAE vs w0')
    ax1.grid(True, alpha=0.3)

    ax2.plot(sens_df['w0'], sens_df['Coverage_95'], 's-', color='#e15759',
            linewidth=2, markersize=7)
    ax2.axhline(y=0.95, color='gray', linestyle='--', alpha=0.5, label='Nominal 95%')
    ax2.axhline(y=0.88, color='gray', linestyle=':', alpha=0.5, label='Lower bound')
    ax2.axhline(y=0.98, color='gray', linestyle=':', alpha=0.5, label='Upper bound')
    ax2.set_xlabel('Robustness weight (w0)')
    ax2.set_ylabel('95% Coverage')
    ax2.set_title('Coverage vs w0')
    ax2.set_ylim(0.5, 1.05)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {output_path}")


def plot_contamination(contam_df, output_path):
    """Grouped bar chart of MAE by contamination scenario."""
    summary = contam_df.groupby('scenario').agg(
        MAE=('abs_error', 'mean'),
        Coverage=('coverage_95', 'mean'),
    ).reset_index()

    # Preserve order
    order = ['Clean', '+10 Adjacent', '+10 Foreign', '+10 Adj +10 For']
    summary = summary.set_index('scenario').loc[order].reset_index()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    x = range(len(order))

    ax1.bar(x, summary['MAE'], color='#4e79a7')
    ax1.set_xticks(x)
    ax1.set_xticklabels(order, rotation=30, ha='right')
    ax1.set_ylabel('MAE')
    ax1.set_title('KG-DAP MAE Under Contamination')

    ax2.bar(x, summary['Coverage'], color='#e15759')
    ax2.axhline(y=0.85, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(order, rotation=30, ha='right')
    ax2.set_ylabel('95% Coverage')
    ax2.set_title('KG-DAP Coverage Under Contamination')
    ax2.set_ylim(0.5, 1.05)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {output_path}")


def plot_ablation(ablation_df, output_path):
    """Horizontal bar chart of MAE by ablation variant."""
    summary = ablation_df.groupby('variant').agg(
        MAE=('abs_error', 'mean'),
        Coverage=('coverage_95', 'mean'),
    ).reset_index()

    # Sort by MAE
    summary = summary.sort_values('MAE', ascending=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = range(len(summary))
    colors = ['#59a14f' if 'Full' in v else '#4e79a7' for v in summary['variant']]
    ax.barh(y_pos, summary['MAE'], color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(summary['variant'])
    ax.set_xlabel('MAE')
    ax.set_title('Ablation Study: MAE by Variant')
    for i, (mae, cov) in enumerate(zip(summary['MAE'], summary['Coverage'])):
        ax.text(mae + 0.002, i, f'MAE={mae:.4f}, Cov={cov:.2f}',
               va='center', fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {output_path}")


def plot_arm_level(results_df, output_path):
    """Scatter plot of prior mean vs observed rate, colored by method."""
    fig, ax = plt.subplots(figsize=(8, 8))
    methods = results_df['method'].unique()
    colors = {'KG-DAP': '#4e79a7', 'rMAP': '#f28e2b', 'EqualWeight': '#76b7b2',
              'Uniform': '#e15759', 'Pooled': '#59a14f'}
    markers = {'KG-DAP': 'o', 'rMAP': 's', 'EqualWeight': '^',
               'Uniform': 'v', 'Pooled': 'D'}

    for method in methods:
        subset = results_df[results_df['method'] == method]
        ax.scatter(subset['observed_rate'], subset['prior_mean'],
                  c=colors.get(method, 'gray'), marker=markers.get(method, 'o'),
                  label=method, alpha=0.7, s=40)

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Perfect prediction')
    ax.set_xlabel('Observed rate')
    ax.set_ylabel('Prior predictive mean')
    ax.set_title('LOO-CV: Prior Mean vs Observed Rate')
    ax.legend()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {output_path}")


# ===================================================================
# Main
# ===================================================================
def main():
    parser = argparse.ArgumentParser(description="KG-DAP Validation")
    parser.add_argument('--quick', action='store_true',
                       help='Run only LOO-CV (skip stress tests)')
    parser.add_argument('--stress', action='store_true',
                       help='Run only stress tests')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    np.random.seed(args.seed)
    rng = np.random.default_rng(seed=args.seed)

    print("=" * 60)
    print("KG-DAP Validation Runner")
    print("=" * 60)

    # Load data
    print("\n[1] Loading data...")
    t0 = time.time()
    trials_df, G, trial_ids = load_data()
    print(f"    {len(trial_ids)} trial arms, {G.number_of_nodes()} graph nodes")
    print(f"    Loaded in {time.time()-t0:.2f}s")

    # ---------------------------------------------------------------
    # QUICK VALIDATION: LOO-CV
    # ---------------------------------------------------------------
    if not args.stress:
        print("\n[2] Running LOO-CV: KG-DAP...")
        t0 = time.time()
        dap_results = run_loocv(trials_df, G, trial_ids,
                                method_name="KG-DAP", seed=args.seed)
        t_dap = time.time() - t0
        print(f"    KG-DAP LOO-CV completed in {t_dap:.2f}s")

        print("\n[3] Running LOO-CV: Comparators...")
        t0 = time.time()
        comp_results = run_loocv_comparators(trials_df, G, trial_ids,
                                             seed=args.seed)
        t_comp = time.time() - t0
        print(f"    Comparators LOO-CV completed in {t_comp:.2f}s")

        # Combine
        all_loocv = pd.concat([dap_results, comp_results], ignore_index=True)
        all_loocv.to_csv(os.path.join(RESULTS_DIR, "loocv_results.csv"),
                         index=False)

        # Summary
        summary = compute_summary(all_loocv)
        summary.to_csv(os.path.join(RESULTS_DIR, "loocv_summary.csv"),
                       index=False)

        print("\n" + "=" * 60)
        print("LOO-CV SUMMARY")
        print("=" * 60)
        print(summary.to_string(index=False))

        # Seed robustness check
        print("\n[4] Seed robustness check (seed=123)...")
        t0 = time.time()
        dap_check = run_loocv(trials_df, G, trial_ids,
                              method_name="KG-DAP_seed123", seed=123)
        mae_42 = dap_results['abs_error'].mean()
        mae_123 = dap_check['abs_error'].mean()
        print(f"    MAE(seed=42)={mae_42:.6f}, MAE(seed=123)={mae_123:.6f}")
        print(f"    Difference: {abs(mae_42-mae_123):.6f} "
              f"({'STABLE' if abs(mae_42-mae_123) < 0.001 else 'UNSTABLE'})")
        print(f"    Completed in {time.time()-t0:.2f}s")

        # Figures
        print("\n[5] Generating figures...")
        plot_loocv_comparison(all_loocv,
                             os.path.join(FIGURES_DIR, "loocv_mae_comparison.png"))
        plot_arm_level(all_loocv,
                      os.path.join(FIGURES_DIR, "loocv_arm_level.png"))

    # ---------------------------------------------------------------
    # STRESS TESTS
    # ---------------------------------------------------------------
    if not args.quick:
        print("\n[6] Stress Test: Sensitivity to beta_diff...")
        t0 = time.time()
        sens_beta_df = sensitivity_beta(trials_df, G, trial_ids)
        sens_beta_df.to_csv(os.path.join(RESULTS_DIR, "sensitivity_beta.csv"),
                           index=False)
        print(f"    Completed in {time.time()-t0:.2f}s")
        plot_sensitivity_beta(sens_beta_df,
                             os.path.join(FIGURES_DIR, "sensitivity_beta.png"))

        print("\n[7] Stress Test: Sensitivity to w0...")
        t0 = time.time()
        sens_w0_df = sensitivity_w0(trials_df, G, trial_ids)
        sens_w0_df.to_csv(os.path.join(RESULTS_DIR, "sensitivity_w0.csv"),
                         index=False)
        print(f"    Completed in {time.time()-t0:.2f}s")
        plot_sensitivity_w0(sens_w0_df,
                           os.path.join(FIGURES_DIR, "sensitivity_w0.png"))

        print("\n[8] Stress Test: Contamination...")
        t0 = time.time()
        adj_df, for_df = load_external_data()
        contam_df = contamination_test(trials_df, G, trial_ids, adj_df, for_df)
        contam_df.to_csv(os.path.join(RESULTS_DIR, "contamination_results.csv"),
                        index=False)
        print(f"    Completed in {time.time()-t0:.2f}s")
        plot_contamination(contam_df,
                          os.path.join(FIGURES_DIR, "contamination.png"))

        print("\n[9] Stress Test: Ablation...")
        t0 = time.time()
        ablation_df = ablation_test(trials_df, G, trial_ids)
        ablation_df.to_csv(os.path.join(RESULTS_DIR, "ablation_results.csv"),
                          index=False)
        print(f"    Completed in {time.time()-t0:.2f}s")
        plot_ablation(ablation_df,
                     os.path.join(FIGURES_DIR, "ablation.png"))

    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)
    print(f"Results saved to: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
