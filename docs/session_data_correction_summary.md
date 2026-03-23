# Session: Data Correction & Manuscript Update

**Date**: 2026-03-23
**Branch**: `auto-research-test`

## Completed

### Step 1: Data Corrections Applied to `data_curation.py`
- **DROPPED** (4 arms): ELOQUENT-1 (2 arms, MRD never reported), SWOG-S1211 (2 arms, MRD never reported)
- **Source-only fixes** (5 trials): Elo-KRd, CEPHEUS, GRIFFIN, MANHATTAN, IsKia
- **Data fixes** (6 trials):
  - ADVANCE: n=176/174 → 148/139, p=0.80/0.65 → 0.59/0.33 (Landgren JCO 2025)
  - DETERMINATION: arms swapped, n=365/357, p=0.545/0.398
  - ENDURANCE: rates 10x off, y=179/188 → 54/38, threshold 1e-5 → 1e-4
  - IFM2009: ASCT y=228 → 220
  - IFM2018-04: source + y=40 → 31 (Touzeau Blood 2024)
  - GMMG-CONCEPT: n=50 → 99, source → Leypoldt JCO 2024
- **Restructured** (2 trials):
  - MIDAS: 2 arms → 1 (induction only, n=791, p=0.63, Perrot Blood 2025)
  - COBRA → Derman-DKRd: fixed identity, NCT, n=42, p=0.619
- **Final dataset**: 35 arms from 21 trials (was 40/23)

### Step 2: Data Regeneration
- `trials_data.csv`: 35 rows
- `similarity_matrices.pkl`: 35×35
- KG: 56 nodes, 201 edges

### Step 3: All Analyses Re-run
- LOO-CV (6 methods, B=2000)
- ESS via ELIR (M=200,000)
- Contamination robustness (20 external trials)
- Sensitivity (rho, w_rob, KG weights)
- Component contribution analysis (ablation)
- Simulation study (5 scenarios, R=1000, B=2000)
- Design OC (R=10,000, nC=25 and nC=50)

### Step 4: Both Manuscripts Updated
- `manuscript-manual-edit.tex` — all numbers, tables, text updated
- `manuscript.tex` — all numbers, tables, text updated
- Red-highlighted versions reviewed, then colors cleared
- Both compile cleanly

## Key Result Changes (Old → New)

| Metric | Old (40 arms) | New (35 arms) |
|--------|--------------|---------------|
| KG-CAR MSE | 0.0433 | 0.0382 |
| rMAP MSE | 0.0490 | 0.0406 |
| MSE reduction | 11.6% | 5.9% |
| CRPS reduction | 7.2% | 3.7% |
| Coverage | 0.950 | 0.886 |
| Wilcoxon p | 0.028 | 0.001 |
| Wins | 24/40 | 26/35 |
| Pairwise corr(W,|Δp|) | -0.50 | +0.11 |
| Sim Sc.2 MSE reduction | 56.7% | 56.5% |
| Sim Sc.3 MSE change | -2.0% | -1.8% |
| Sim Sc.4 MSE reduction | 37.2% | 37.4% |
| ESS_ELIR mean | 5.3 | 5.7 |
| Frac structured | 0.544 | 0.499 |

### Coverage Analysis
4 arms uncovered by ALL borrowing methods (KG-CAR, rMAP, BHM):
- ALCYONE VMP (p=0.070)
- ENDURANCE VRd (p=0.072)
- ENDURANCE KRd (p=0.103)
- MAIA Rd (p=0.111)

All 4 are extreme-low-rate arms below the predictive distribution's lower bound.

### Correlation Change
Pairwise r(W_ij, |p_i-p_j|) went from -0.50 to +0.11. Driven by ENDURANCE correction (34%→10%/7%): these VRd/KRd arms are structurally similar to high-rate arms but have very different outcomes. Reframed in manuscript as motivation for multivariate CAR (captures conditional dependence beyond pairwise correlation).

## Not Yet Done
- [ ] Re-run full research-pipeline skill on corrected data (user will run in separate terminal)
- [ ] Assemble ARMS GitHub repo (blocked on pipeline re-run)

## Pipeline Command
```bash
cd "/Users/yuanj/Dropbox (Personal)/YuanJi/research/Yunxuan-Zhang/KG"
claude -p "/research-pipeline research_brief_v2.md" --dangerouslySkipPermissions
```

## Files Modified
- `data_curation.py` — trial data corrections
- `data/trials_data.csv` — regenerated (35 arms)
- `data/similarity_matrices.pkl` — regenerated (35×35)
- `pipeline/phase2_validate/validated_results/*.csv` — all re-run
- `pipeline/phase3_write/data/*.csv` — simulation results
- `pipeline/phase4_polish/round_1/manuscript-manual-edit.tex` — all numbers updated
- `pipeline/phase4_polish/round_1/manuscript.tex` — all numbers updated
- `pipeline/phase4_polish/round_1/*.pdf` — recompiled
