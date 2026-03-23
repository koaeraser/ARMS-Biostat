"""
Production re-run of KG-CAR validation with B=2000 and ESS analysis.

Reuses functions from run_validation.py but with publication-quality settings.
Also computes effective sample size (ESS) and variance decomposition.
"""

import sys
import os
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.special import expit, logit as sp_logit

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent.parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR.parent / "validated_results"

sys.path.insert(0, str(SCRIPT_DIR))

from method import KGCAR, KGCARNoRob, KGCARNoBYM, build_adjacency_and_laplacian
from comparator import NoBorrowing, StandardBHM, RobustMAP, KGPP, FullPooling
from run_validation import (load_data, crps_sample, run_loocv,
                            run_contamination_test, run_sensitivity_analysis,
                            run_ablation, load_external_data, make_figures)


B_PROD = 2000
SEED = 42


def compute_ess(model, h_target, y, n, W):
    """
    Compute effective sample size borrowed for arm h_target.

    ESS_h = rho * sum_j(W_hj * n_j) / (rho * sum_j(W_hj) + (1 - rho))
    modulated by sigma_phi^2 / (sigma_phi^2 + sigma_eps^2).

    Also returns variance decomposition:
    - frac_structured = sigma_phi^2 / (sigma_phi^2 + sigma_eps^2)
    - frac_unstructured = 1 - frac_structured
    """
    H = len(y)
    mask = np.ones(H, dtype=bool)
    mask[h_target] = False

    # Fit on remaining arms
    model_copy = KGCAR(W, model.L, w_rob=model.w_rob)
    model_copy.fit(y, n, mask=mask)

    rho = model_copy.rho_hat
    sigma_phi = model_copy.sigma_phi_hat
    sigma_eps = model_copy.sigma_eps_hat

    W_row = W[h_target, mask]
    n_hist = n[mask]

    w_sum = W_row.sum()
    denom = rho * w_sum + (1 - rho)

    # Raw ESS from graph structure (in units of "equivalent patients")
    ess_raw = rho * np.sum(W_row * n_hist) / max(denom, 1e-10)

    # Modulate by variance fraction: how much of the total variance is structured
    var_total = sigma_phi**2 + sigma_eps**2
    frac_structured = sigma_phi**2 / max(var_total, 1e-10)

    # Effective ESS = raw ESS * fraction explained by structure
    ess_effective = ess_raw * frac_structured

    # Per-neighbor ESS contributions (top contributors)
    neighbor_ess = rho * W_row * n_hist / max(denom, 1e-10) * frac_structured

    return {
        'ess_raw': ess_raw,
        'ess_effective': ess_effective,
        'rho': rho,
        'sigma_phi': sigma_phi,
        'sigma_eps': sigma_eps,
        'frac_structured': frac_structured,
        'w_sum': w_sum,
        'neighbor_ess': neighbor_ess,
    }


def run_ess_analysis(y, n, p_hat, trial_ids, W, L):
    """Compute ESS for all 40 arms via LOO."""
    H = len(y)
    model = KGCAR(W, L, w_rob=0.10)
    records = []

    print("\n  Computing ESS for each arm...")
    t0 = time.time()

    for h in range(H):
        ess_result = compute_ess(model, h, y, n, W)

        # Top 5 neighbors by ESS contribution
        neighbor_ess = ess_result['neighbor_ess']
        mask = np.ones(H, dtype=bool)
        mask[h] = False
        hist_ids = [trial_ids[j] for j in range(H) if j != h]
        top_idx = np.argsort(neighbor_ess)[::-1][:5]
        top_neighbors = [(hist_ids[i], neighbor_ess[i]) for i in top_idx]

        records.append({
            'trial_id': trial_ids[h],
            'arm_idx': h,
            'n_h': n[h],
            'p_hat': p_hat[h],
            'ess_raw': ess_result['ess_raw'],
            'ess_effective': ess_result['ess_effective'],
            'rho': ess_result['rho'],
            'sigma_phi': ess_result['sigma_phi'],
            'sigma_eps': ess_result['sigma_eps'],
            'frac_structured': ess_result['frac_structured'],
            'w_sum': ess_result['w_sum'],
            'top1_neighbor': top_neighbors[0][0],
            'top1_ess': top_neighbors[0][1],
            'top2_neighbor': top_neighbors[1][0],
            'top2_ess': top_neighbors[1][1],
            'top3_neighbor': top_neighbors[2][0],
            'top3_ess': top_neighbors[2][1],
        })

    print(f"    Done in {time.time()-t0:.1f}s")

    ess_df = pd.DataFrame(records)
    return ess_df


def add_mc_standard_errors(results_df):
    """Add MC standard errors to a summary dataframe."""
    summary = results_df.groupby('method').agg(
        MSE=('se', 'mean'),
        MSE_SE=('se', lambda x: x.std() / np.sqrt(len(x))),
        MAE=('ae', 'mean'),
        MAE_SE=('ae', lambda x: x.std() / np.sqrt(len(x))),
        RMSE=('se', lambda x: np.sqrt(x.mean())),
        Coverage=('covered', 'mean'),
        Coverage_SE=('covered', lambda x: x.std() / np.sqrt(len(x))),
        CRPS=('crps', 'mean'),
        CRPS_SE=('crps', lambda x: x.std() / np.sqrt(len(x))),
        CI_Width=('ci_width', 'mean'),
        CI_Width_SE=('ci_width', lambda x: x.std() / np.sqrt(len(x))),
        N_valid=('se', 'count'),
    ).reset_index()
    return summary


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 70)
    print("KG-CAR PRODUCTION VALIDATION (B=2000)")
    print("=" * 70)

    # Load data
    print("\n[Step 0] Loading data...")
    trials, y, n, p_hat, trial_ids, S_composite, J_drug, J_target, S_pop = load_data()
    W, L = build_adjacency_and_laplacian(S_composite)
    H = len(y)
    print(f"  {H} arms, n=[{n.min():.0f}, {n.max():.0f}], p=[{p_hat.min():.3f}, {p_hat.max():.3f}]")

    # =============================================
    # Step 1: LOO-CV (B=2000)
    # =============================================
    print(f"\n[Step 1] LOO-CV (B={B_PROD})...")
    methods = {
        'No Borrowing': NoBorrowing(),
        'Full Pooling': FullPooling(),
        'BHM': StandardBHM(),
        'rMAP': RobustMAP(w_rob=0.10),
        'KG-PP': KGPP(W),
        'KG-CAR': KGCAR(W, L, w_rob=0.10),
    }

    loocv_results, _ = run_loocv(methods, y, n, p_hat, trial_ids, W, B=B_PROD, seed=SEED)
    loocv_summary = add_mc_standard_errors(loocv_results)

    loocv_results.to_csv(RESULTS_DIR / 'loocv_results.csv', index=False)
    loocv_summary.to_csv(RESULTS_DIR / 'loocv_summary.csv', index=False)

    print("\n  === LOO-CV Summary (B=2000) ===")
    for _, row in loocv_summary.iterrows():
        print(f"  {row['method']:>14s}  "
              f"MSE={row['MSE']:.5f}({row['MSE_SE']:.5f})  "
              f"Cov={row['Coverage']:.3f}({row['Coverage_SE']:.3f})  "
              f"CRPS={row['CRPS']:.5f}({row['CRPS_SE']:.5f})  "
              f"CIW={row['CI_Width']:.3f}")

    # =============================================
    # Step 2: ESS Analysis
    # =============================================
    print(f"\n[Step 2] ESS Analysis...")
    ess_df = run_ess_analysis(y, n, p_hat, trial_ids, W, L)
    ess_df.to_csv(RESULTS_DIR / 'ess_analysis.csv', index=False)

    print(f"\n  === ESS Summary ===")
    print(f"  ESS (effective): mean={ess_df['ess_effective'].mean():.1f}, "
          f"median={ess_df['ess_effective'].median():.1f}, "
          f"range=[{ess_df['ess_effective'].min():.1f}, {ess_df['ess_effective'].max():.1f}]")
    print(f"  Frac structured: mean={ess_df['frac_structured'].mean():.3f}, "
          f"sd={ess_df['frac_structured'].std():.3f}")
    print(f"  rho: mean={ess_df['rho'].mean():.3f}, sd={ess_df['rho'].std():.3f}")
    print(f"  sigma_phi: mean={ess_df['sigma_phi'].mean():.3f}")
    print(f"  sigma_eps: mean={ess_df['sigma_eps'].mean():.3f}")

    # =============================================
    # Step 3: Stress tests (B=2000)
    # =============================================
    print(f"\n[Step 3] Stress tests (B={B_PROD})...")

    # 3a: Contamination
    print("\n  [3a] Contamination test...")
    trials_adj, trials_foreign = load_external_data()
    contamination_results, contamination_summary = run_contamination_test(
        methods, y, n, p_hat, trial_ids, S_composite,
        trials, trials_adj, trials_foreign, B=B_PROD, seed=SEED
    )
    contamination_results.to_csv(RESULTS_DIR / 'stress_test_contamination.csv', index=False)

    # Robustness ratio
    clean_mse = loocv_summary.set_index('method')['MSE']
    dirty_mse = contamination_summary.set_index('method')['MSE']
    robustness = pd.DataFrame()
    for method in clean_mse.index:
        if method in dirty_mse.index:
            robustness = pd.concat([robustness, pd.DataFrame([{
                'method': method,
                'MSE_clean': clean_mse[method],
                'MSE_contaminated': dirty_mse[method],
                'robustness_ratio': clean_mse[method] / max(dirty_mse[method], 1e-10),
                'MSE_change_pct': (dirty_mse[method] - clean_mse[method]) / max(clean_mse[method], 1e-10) * 100,
            }])], ignore_index=True)
    robustness.to_csv(RESULTS_DIR / 'stress_test_robustness.csv', index=False)
    print("\n  Robustness ratios:")
    print(robustness.to_string(index=False))

    # 3b: Sensitivity
    print("\n  [3b] Sensitivity analysis...")
    sensitivity_results = run_sensitivity_analysis(
        y, n, p_hat, trial_ids, S_composite, J_drug, J_target, S_pop,
        B=B_PROD, seed=SEED
    )
    sensitivity_results.to_csv(RESULTS_DIR / 'stress_test_sensitivity.csv', index=False)

    # 3c: Ablation
    print("\n  [3c] Ablation study...")
    ablation_results, ablation_summary = run_ablation(
        y, n, p_hat, trial_ids, W, L, B=B_PROD, seed=SEED
    )
    ablation_results.to_csv(RESULTS_DIR / 'stress_test_ablation.csv', index=False)
    print("\n  Ablation summary:")
    print(ablation_summary.to_string(index=False))

    # =============================================
    # Step 4: Figures
    # =============================================
    print("\n[Step 4] Generating figures...")
    make_figures(loocv_results, sensitivity_results, ablation_summary,
                 loocv_summary, contamination_summary)

    # =============================================
    # Final
    # =============================================
    print("\n" + "=" * 70)
    print("PRODUCTION VALIDATION COMPLETE (B=2000)")
    print("=" * 70)

    # Final comparison
    kg = loocv_summary[loocv_summary['method'] == 'KG-CAR'].iloc[0]
    rm = loocv_summary[loocv_summary['method'] == 'rMAP'].iloc[0]
    mse_red = (1 - kg['MSE'] / rm['MSE']) * 100
    print(f"\n  KG-CAR vs rMAP: MSE reduction = {mse_red:+.1f}%")
    print(f"  KG-CAR Coverage: {kg['Coverage']:.3f} ({kg['Coverage_SE']:.3f})")
    print(f"  ESS effective: mean={ess_df['ess_effective'].mean():.1f}")
    print(f"\n  Results: {RESULTS_DIR}/")


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings('ignore')
    main()
