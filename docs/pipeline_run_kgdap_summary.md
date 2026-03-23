# KG-DAP Pipeline Run Summary

**Date**: 2026-03-23 (20:14 - 23:28)
**Duration**: ~3h 14min
**Branch**: `auto-research-test`
**Brief**: `research_brief_v2.md` (method-agnostic, no reference summaries)
**Outcome**: SUCCESS (42/50 self-graded, 38/50 independent grade)

## Method Discovered: KG-DAP
Knowledge-Graph-Driven Adaptive Prior — a finite Beta mixture prior with
power-sharpened composite similarity weights from a biomedical KG.

Key design decisions made by the pipeline:
- Started with diffusion kernel (from methodology spec)
- Phase 2 discovered spectral compression issue on dense graphs
- Pivoted to power-based weights (omega_k = S^beta / sum S^beta)
- Single hyperparameter beta controls borrowing concentration
- Fully pre-specifiable, closed-form, conjugate to binomial

## Phase Timing
| Phase | Duration | Key Output |
|-------|----------|-----------|
| 1. THINK | ~68 min | methodology_specification.md (481 lines) |
| 2. VALIDATE | ~15 min | GO verdict, 25.9% MAE improvement |
| 3. WRITE | ~46 min | manuscript.tex (645 lines), 6 figs, 20 refs |
| 4. POLISH | ~62 min | 3 rounds: 40→41→42/50 |
| **Total** | **~3h 14min** | |

## Key Results
- MAE: 0.126 (KG-DAP) vs 0.170 (rMAP) — 25.9% improvement
- Coverage: 94.3% (well-calibrated)
- ESS: ~5.6 (modest, intentional)
- Simulation: favorable +25.9%, adverse -0.8%, mixed +16.2%
- Design OC: +7-12pp go probability at clinically relevant effects

## Score Breakdown (Pipeline Self-Grade vs Independent)
| Dimension | Self | Independent |
|-----------|------|------------|
| Correctness | 5 | 4 |
| Completeness | 5 | 4 |
| Rigor | 5 | 4 |
| Clarity | 5 | 4 |
| Novelty | 4 | 3 |
| Impact | 4 | 4 |
| Performance | 4 | 4 |
| **Weighted** | **42/50** | **38/50** |

Self-grading inflation: +4 points (consistent with known 1-4 point pattern).

## Estimated Token Consumption
- ~15-18 agent conversations (orchestrator + phase agents + sub-agents)
- Model: Claude Opus 4.6 (1M context)
- No exact token counts available; check claude.ai/settings/usage

## Files Produced
- pipeline/phase1_think/: methodology spec, rationale, literature briefing
- pipeline/phase2_validate/: validation report, kg_dap.py, comparators.py, 8 CSVs, 6 PNGs
- pipeline/phase3_write/: manuscript.tex, 6 PDF figs, 8 CSVs, references.bib
- pipeline/phase4_polish/round_2/: final manuscript.tex + PDF (30 pages)
- pipeline/final_report.md, pipeline_log.md, pipeline_state.md
