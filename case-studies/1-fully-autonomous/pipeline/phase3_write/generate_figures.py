"""
Production Figure Generation for KG-DAP Manuscript
====================================================
Generates publication-quality PDF figures from validated results.
Uses the validated code as foundation - does NOT rewrite anything.
"""

import os
import sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# Publication style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'lines.linewidth': 1.5,
    'lines.markersize': 6,
})

# Color palette (colorblind-friendly)
COLORS = {
    'KG-DAP': '#0072B2',      # blue
    'rMAP': '#D55E00',         # vermillion
    'EqualWeight': '#009E73',  # green
    'Uniform': '#CC79A7',      # pink
    'Pooled': '#F0E442',       # yellow
    'accent': '#56B4E9',       # light blue
}

METHOD_ORDER = ['KG-DAP', 'rMAP', 'EqualWeight', 'Uniform', 'Pooled']


def fig1_loocv_mae_bar():
    """Figure 1: LOO-CV MAE comparison bar chart."""
    summary = pd.read_csv(os.path.join(DATA_DIR, "loocv_summary.csv"))

    # Reorder
    summary = summary.set_index('method').loc[METHOD_ORDER].reset_index()

    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    colors = [COLORS[m] for m in summary['method']]
    bars = ax.bar(range(len(summary)), summary['MAE'], color=colors,
                  edgecolor='black', linewidth=0.5, width=0.7)

    ax.set_xticks(range(len(summary)))
    labels = ['KG-DAP', 'rMAP', 'Equal\nWeight', 'Uniform', 'Pooled']
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Mean Absolute Error (MAE)')
    ax.set_ylim(0, 0.22)
    ax.set_title('Leave-One-Out Cross-Validation: MAE by Method', fontsize=11)

    # Add value labels
    for bar, val in zip(bars, summary['MAE']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.003,
                f'{val:.4f}', ha='center', va='bottom', fontsize=8)

    # Add horizontal reference line at KG-DAP
    ax.axhline(y=summary['MAE'].iloc[0], color='gray', linestyle='--',
               alpha=0.4, linewidth=0.8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig1_loocv_mae.pdf"))
    plt.close()
    print("  fig1_loocv_mae.pdf")


def fig2_arm_level_scatter():
    """Figure 2: Arm-level scatter plot of prior mean vs observed rate."""
    results = pd.read_csv(os.path.join(DATA_DIR, "loocv_results.csv"))

    fig, ax = plt.subplots(figsize=(5.0, 5.0))

    markers = {'KG-DAP': 'o', 'rMAP': 's', 'EqualWeight': '^',
               'Uniform': 'v', 'Pooled': 'D'}

    for method in METHOD_ORDER:
        subset = results[results['method'] == method]
        ax.scatter(subset['observed_rate'], subset['prior_mean'],
                   c=COLORS[method], marker=markers[method],
                   label=method, alpha=0.65, s=30, edgecolors='black',
                   linewidth=0.3, zorder=3 if method == 'KG-DAP' else 2)

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=0.8,
            label='Perfect prediction', zorder=1)

    ax.set_xlabel('Observed MRD Negativity Rate')
    ax.set_ylabel('Prior Predictive Mean')
    ax.set_title('LOO-CV: Prior Mean vs. Observed Rate', fontsize=11)
    ax.legend(loc='upper left', frameon=True, framealpha=0.9, fontsize=8)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig2_arm_level.pdf"))
    plt.close()
    print("  fig2_arm_level.pdf")


def fig3_sensitivity_beta():
    """Figure 3: MAE and Coverage vs beta (dual axis)."""
    sens = pd.read_csv(os.path.join(DATA_DIR, "sensitivity_beta.csv"))

    fig, ax1 = plt.subplots(figsize=(5.0, 3.5))

    # MAE on left axis
    line1, = ax1.plot(sens['beta_diff'], sens['MAE'], 'o-',
                      color=COLORS['KG-DAP'], linewidth=2, markersize=6,
                      label='MAE')
    ax1.set_xlabel(r'Power parameter $\beta$')
    ax1.set_ylabel('MAE', color=COLORS['KG-DAP'])
    ax1.tick_params(axis='y', labelcolor=COLORS['KG-DAP'])
    ax1.set_xscale('log')
    ax1.set_ylim(0.11, 0.17)

    # Coverage on right axis
    ax2 = ax1.twinx()
    line2, = ax2.plot(sens['beta_diff'], sens['Coverage_95'], 's--',
                      color=COLORS['rMAP'], linewidth=1.5, markersize=5,
                      label='95% Coverage')
    ax2.set_ylabel('95% Coverage', color=COLORS['rMAP'])
    ax2.tick_params(axis='y', labelcolor=COLORS['rMAP'])
    ax2.set_ylim(0.85, 1.05)

    # Reference lines
    ax2.axhline(y=0.95, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)
    ax2.axhline(y=0.88, color='gray', linestyle=':', alpha=0.3, linewidth=0.6)
    ax2.axhline(y=0.98, color='gray', linestyle=':', alpha=0.3, linewidth=0.6)

    # Shade safe operating range
    ax1.axvspan(8, 30, alpha=0.08, color='green', label=r'Safe range $\beta \in [8, 30]$')

    # Combined legend
    lines = [line1, line2]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper right', fontsize=8, framealpha=0.9)

    ax1.set_title(r'Sensitivity to Power Parameter $\beta$', fontsize=11)
    ax1.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig3_sensitivity_beta.pdf"))
    plt.close()
    print("  fig3_sensitivity_beta.pdf")


def fig4_ablation():
    """Figure 4: Ablation study horizontal bar chart."""
    ablation = pd.read_csv(os.path.join(DATA_DIR, "ablation_results.csv"))

    summary = ablation.groupby('variant').agg(
        MAE=('abs_error', 'mean'),
        Coverage=('coverage_95', 'mean'),
    ).reset_index()
    summary = summary.sort_values('MAE', ascending=True)

    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    y_pos = range(len(summary))

    colors = []
    for v in summary['variant']:
        if 'Full' in v:
            colors.append(COLORS['KG-DAP'])
        elif 'No Graph' in v or 'Diffusion' in v:
            colors.append(COLORS['rMAP'])
        else:
            colors.append(COLORS['accent'])

    bars = ax.barh(y_pos, summary['MAE'], color=colors,
                   edgecolor='black', linewidth=0.5, height=0.65)
    ax.set_yticks(y_pos)

    # Clean variant labels
    clean_labels = {
        'Full KG-DAP': 'Full KG-DAP',
        'No Graph (EqualWeight)': 'No Graph (Equal Wt)',
        'No Robustness (w0=0)': r'No Robustness ($w_0 \approx 0$)',
        'No Cap (n_cap=inf)': r'No Cap ($n_{\mathrm{cap}} = \infty$)',
        'Drug-Only Sim': 'Drug-Only Similarity',
        'Pop-Only Sim': 'Population-Only Similarity',
        'Diffusion Kernel': 'Diffusion Kernel (original)',
    }
    ax.set_yticklabels([clean_labels.get(v, v) for v in summary['variant']],
                       fontsize=8)
    ax.set_xlabel('MAE')
    ax.set_title('Component Contribution Analysis', fontsize=11)

    # Add MAE values
    for i, (mae, cov) in enumerate(zip(summary['MAE'], summary['Coverage'])):
        ax.text(mae + 0.002, i, f'{mae:.4f}', va='center', fontsize=7)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig4_ablation.pdf"))
    plt.close()
    print("  fig4_ablation.pdf")


def fig5_contamination():
    """Figure 5: Contamination robustness."""
    contam = pd.read_csv(os.path.join(DATA_DIR, "contamination_stress_results.csv"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.5))

    # Filter to show sim=0.8 and sim=0.9 separately
    sim08 = contam[(contam['fake_sim'] == 0.8) | (contam['scenario'] == 'Clean')]
    sim09 = contam[(contam['fake_sim'] == 0.9) | (contam['scenario'] == 'Clean')]

    # Plot sim=0.8
    x_08 = [0, 5, 5, 10, 10]
    mae_08 = [contam[contam['scenario'] == 'Clean']['MAE'].values[0]]
    labels_08 = ['Clean']
    for _, row in contam[contam['fake_sim'] == 0.8].iterrows():
        mae_08.append(row['MAE'])
        labels_08.append(f"+{int(row['n_fake'])}\n(r={row['fake_rate']})")

    ax1.bar(range(len(mae_08)), mae_08, color=COLORS['KG-DAP'],
            edgecolor='black', linewidth=0.5, width=0.6)
    ax1.set_xticks(range(len(mae_08)))
    ax1.set_xticklabels(labels_08, fontsize=7, rotation=0)
    ax1.set_ylabel('MAE')
    ax1.set_title('Contamination (sim = 0.80)', fontsize=10)
    ax1.axhline(y=0.1257 + 0.05, color='red', linestyle='--', alpha=0.5,
                linewidth=0.8, label='5% degradation threshold')
    ax1.set_ylim(0.10, 0.20)
    ax1.legend(fontsize=7)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Plot sim=0.9
    mae_09 = [contam[contam['scenario'] == 'Clean']['MAE'].values[0]]
    labels_09 = ['Clean']
    for _, row in contam[contam['fake_sim'] == 0.9].iterrows():
        mae_09.append(row['MAE'])
        labels_09.append(f"+{int(row['n_fake'])}\n(r={row['fake_rate']})")

    ax2.bar(range(len(mae_09)), mae_09, color=COLORS['rMAP'],
            edgecolor='black', linewidth=0.5, width=0.6)
    ax2.set_xticks(range(len(mae_09)))
    ax2.set_xticklabels(labels_09, fontsize=7, rotation=0)
    ax2.set_ylabel('MAE')
    ax2.set_title('Contamination (sim = 0.90)', fontsize=10)
    ax2.axhline(y=0.1257 + 0.05, color='red', linestyle='--', alpha=0.5,
                linewidth=0.8, label='5% degradation threshold')
    ax2.set_ylim(0.10, 0.25)
    ax2.legend(fontsize=7)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig5_contamination.pdf"))
    plt.close()
    print("  fig5_contamination.pdf")


def fig6_sensitivity_w0():
    """Figure 6: Sensitivity to w0 (supplementary)."""
    sens = pd.read_csv(os.path.join(DATA_DIR, "sensitivity_w0.csv"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.5))

    ax1.plot(sens['w0'], sens['MAE'], 'o-', color=COLORS['KG-DAP'],
             linewidth=2, markersize=6)
    ax1.set_xlabel(r'Robustness weight $w_0$')
    ax1.set_ylabel('MAE')
    ax1.set_title(r'MAE vs. $w_0$', fontsize=10)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(True, alpha=0.2)

    ax2.plot(sens['w0'], sens['Coverage_95'], 's-', color=COLORS['rMAP'],
             linewidth=2, markersize=6)
    ax2.axhline(y=0.95, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax2.axhspan(0.88, 0.98, alpha=0.08, color='green', label='Target range')
    ax2.set_xlabel(r'Robustness weight $w_0$')
    ax2.set_ylabel('95% Coverage')
    ax2.set_title(r'Coverage vs. $w_0$', fontsize=10)
    ax2.set_ylim(0.80, 1.05)
    ax2.legend(fontsize=8)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig6_sensitivity_w0.pdf"))
    plt.close()
    print("  fig6_sensitivity_w0.pdf")


if __name__ == "__main__":
    print("Generating publication-quality figures (PDF)...")
    fig1_loocv_mae_bar()
    fig2_arm_level_scatter()
    fig3_sensitivity_beta()
    fig4_ablation()
    fig5_contamination()
    fig6_sensitivity_w0()
    print("Done. All figures saved to:", FIGURES_DIR)
