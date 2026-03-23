# Case Studies

Two papers were produced from the same research brief and corrected dataset (35 NDMM trial arms, 21 studies), using different workflows. Both target *Biometrics*.

## Comparison

|  | 1. Fully Autonomous | 2. Human-Revised |
|--|---------------------|------------------|
| **Method** | KG-DAP (Beta mixture + power-sharpened KG weights) | KG-CAR (Leroux CAR prior + BYM decomposition) |
| **Workflow** | `/research-pipeline` end-to-end, zero human input | ~13 AI sessions + ~6h human revision |
| **Duration** | ~3h | ~3h + ~6h human revision |
| **Self-graded score** | 42/50 (Grade A) | 42/50 (Grade A) |
| **Independent score** | 38/50 (Grade B) | 39/50 (Grade B) |
| **Coverage** | 94.3% | 88.6% |
| **Theory depth** | 4 propositions, simple proofs | 2 propositions + 2 remarks, deeper proofs |
| **Comparators** | 4 methods | 6 methods |
| **Simulation** | 3 scenarios, B=1000 | 5 scenarios, R=1000, B=2000 |
| **Design OC** | Go/no-go power curves | Permutation-based centering study |
| **Key strength** | Larger empirical advantage, fully closed-form | Deeper novelty (spatial stats bridge), more comprehensive evaluation |
| **Key weakness** | Methodologically simpler (weighted Beta mixture) | Smaller real-data advantage |

## What to look at

### 1-fully-autonomous/
The complete `pipeline/` directory as produced by the research-pipeline skill. Every file was written by Claude Code with no human intervention. Start with:
- `pipeline/final_report.md` -- outcome summary and lessons learned
- `pipeline/pipeline_log.md` -- phase-by-phase decision log
- `pipeline/phase4_polish/round_2/manuscript.pdf` -- final paper

### 2-human-revised/
The KG-CAR paper developed over ~13 sessions of human-AI collaboration, then revised manually. Contains:
- `manuscript.tex` / `.pdf` -- human-revised final version
- `manuscript-auto.tex` / `.pdf` -- LLM-generated version before human editing
- `code/` -- validated analysis scripts (method.py, comparator.py, etc.)
- `results/` -- all analysis output CSVs

## Key Observations

1. **The pipeline independently discovered a different (and arguably better-performing) method.** Given the same problem brief and data, the methodology architect chose KG-DAP (Beta mixture with power-sharpened weights) rather than KG-CAR (spatial CAR prior). Neither method was prescribed in the brief.

2. **Self-grading inflates by 1-4 points.** Both papers self-graded at 42/50; independent assessment gives 38-39/50. This is consistent across all pipeline runs.

3. **Autonomous vs human-revised quality gap is small (1 point).** The pipeline produces a competent paper (38/50) that is close to but slightly below human-revised quality (39/50). The gap is in novelty (simpler method) and completeness (fewer evaluation scenarios), not in correctness or clarity.

4. **Speed vs depth tradeoff.** The autonomous pipeline produces a complete paper in 3 hours. The human-revised version took weeks but has deeper theoretical contribution and more comprehensive evaluation.
