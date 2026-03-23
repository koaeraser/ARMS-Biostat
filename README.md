# Research Pipeline: Autonomous Academic Paper Production with Claude Code

A system of 10 coordinated Claude Code skills that automate the full lifecycle of a statistical methodology paper — from idea to polished manuscript. Built and tested on a Bayesian clinical trial methods paper targeting *Biometrics*.

## Why This Exists

### The Problem

We attempted to automate paper writing with a simple loop: write a full paper, grade it, fix the weak spots, repeat. Across three test runs, this approach consistently **plateaued at ~65–75% of the target quality** and could not push further:

| Run | Cycles | Score Trajectory | Outcome |
|-----|--------|-----------------|---------|
| 1   | 3      | 16 → 17 → 19    | Plateau |
| 2   | 5      | 19 → 20 → 21 → 21 → 21 | Plateau |
| 3   | 3      | 22 → 21 → 23    | Plateau |

**Root cause:** the system wrote a full paper *before* validating whether the method actually worked. Once the paper existed, the polish loop could only fix surface issues (tables, prose, citations) — it couldn't fix fundamental methodology problems. An RMSE-vs-MAE bug survived three polish cycles undetected because the grader reviews LaTeX, not Python.

### The Insight

> Auto-research is a **gap-closing** system, not a **frontier-pushing** system. It can close the gap between a validated method and a polished paper, but it cannot close the gap between an unvalidated idea and a valid contribution.

This led to a complete architectural redesign: **validate the method before writing anything**.

## Architecture

```
research-pipeline (outer orchestrator)
├── Phase 1: THINK     → methodology-architect
├── Phase 2: VALIDATE  → validate-method
├── Phase 3: WRITE     → write-manuscript
│   ├── literature-lead
│   ├── paper-modeler
│   ├── paper-writer
│   └── paper-critic
└── Phase 4: POLISH    → paper-grader + paper-fixer
```

Each phase is a **separate agent invocation** with its own context window. Communication between phases happens exclusively through **files on disk** (the "anti-telephone-game" pattern — no information is passed through agent summaries that could degrade).

### Phase Gates and Feedback Loops

- **Phase 2 → Phase 1 (RETHINK):** If validation fails, the failure diagnosis feeds back to the methodology architect. Max 2 rethink cycles.
- **Phase 4 (POLISH loop):** Grade → fix → re-grade, up to 3 rounds. Fixes are restricted to **execution quality only** (formulas, tables, citations, clarity) — not methodology.
- **Kill criteria:** The system can honestly report failure (NO-GO, KILL) rather than producing a paper about a method that doesn't work.

## The 10 Skills

### Phase 1: THINK

| Skill | Role | Lines |
|-------|------|-------|
| **`research-pipeline`** | Outer orchestrator. Manages the 4-phase flow, checks phase gates, handles RETHINK loops, writes pipeline logs. Does no research itself. | 564 |
| **`methodology-architect`** | Senior researcher agent. Reads literature (via parallel subagent readers), reasons about method combinations using a Provides/Needs matrix, stress-tests candidates, assesses societal impact, and produces a formal methodology specification. Includes a "wildcard search" phase that looks in adjacent fields (ML, spatial stats, decision theory) for importable ideas. | 659 |

### Phase 2: VALIDATE

| Skill | Role | Lines |
|-------|------|-------|
| **`validate-method`** | Validation scientist. Implements the proposed method and comparators, runs quick experiments (B=200), judges results against success criteria, then stress-tests (B=500) for robustness. Max 3 iterations with structured diagnosis. Includes a 9-point code audit checklist and the "Beauvais Rule" — don't iterate past fundamental limits. | 432 |

### Phase 3: WRITE

| Skill | Role | Lines |
|-------|------|-------|
| **`write-manuscript`** | Manuscript orchestrator. Plans production evaluations, dispatches sub-agents in sequence (literature → modeling → writing → critique), enforces scope constraints. The method is validated; its job is exposition, not discovery. | 553 |
| **`literature-lead`** | Coordinates parallel paper readers for a **writing-focused** review (positioning and citations, not method design). Produces a synthesized briefing with comparative tables, gap analysis, and a draft Related Work section. | 330 |
| **`paper-modeler`** | Scales validated code to publication quality (B=5000). Runs formula-code consistency checks, generates PDF figures, and produces a modeling briefing. Reuses Phase 2 code — extends, doesn't rewrite. | 211 |
| **`paper-writer`** | Writes LaTeX manuscript sections one at a time. Includes mandatory protocols: anomaly detection (scan all results for unexplained patterns), claim-vs-data verification (every interpretive claim checked against CSVs), and ablation interpretation. | 358 |
| **`paper-critic`** | Adversarial reviewer (up to 3 rounds). Severity-based escalation: Critical issues must be fixed, Major issues should be fixed, Minor issues are logged. Scope-limited to exposition quality — cannot re-litigate validated methodology. | 148 |

### Phase 4: POLISH

| Skill | Role | Lines |
|-------|------|-------|
| **`paper-grader`** | Associate editor agent. Scores the manuscript on 7 dimensions (Correctness, Completeness, Rigor, Clarity, Novelty, Impact, Performance) using a calibrated rubric. Dispatches a code auditor subagent for a 15-item computational correctness check. Research dimensions (Novelty, Impact, Performance) carry 2× weight. | 419 |
| **`paper-fixer`** | Copy-editor agent. Applies targeted fixes from the grade report. Strict whitelist/blacklist: may fix formula transcription errors, table mismatches, broken references, and clarity issues. May NOT change methodology, re-run experiments, delete unfavorable results, or add new evaluations. | 298 |

## Key Design Decisions

### 1. Validate Before Writing
The single most important decision. Phase 2 must return GO before any LaTeX is produced. This prevents the v1 failure mode of polishing a paper about a broken method.

### 2. Files on Disk, Not Agent Memory
All inter-phase communication uses files. The methodology spec, validation report, and result CSVs are written to disk and read from disk. No agent summary is passed through an intermediate agent — every consumer reads the primary source. This eliminates the "telephone game" degradation we observed in v1.

### 3. Scope Constraints at Every Level
Each skill has an explicit whitelist (what it MAY do) and blacklist (what it MAY NOT do). The paper-fixer cannot change methodology. The paper-writer cannot re-run experiments. The grader cannot modify files. These constraints prevent scope creep that blurs phase boundaries.

### 4. Honest Failure as a First-Class Outcome
The pipeline has 6 possible outcomes: SUCCESS, VALIDATED_BELOW_TARGET, MARGINAL (surface to user), NO_GO, KILL, and PHASE_1_INCOMPLETE. Four of these are graceful failures. The v1 system always produced a paper, even about methods that didn't work.

### 5. Reuse Before Rewrite
A hard lesson from Run 3: the modeler rewrote an ESS formula from the spec instead of importing the validated code, introducing a transcription bug that survived to the final paper. The v2 system enforces code reuse at every phase.

### 6. The Beauvais Rule
Named after Beauvais Cathedral (which collapsed from iterating past structural limits): if the method fundamentally doesn't work, stop. Three runs of the same experiment with minor parameter tweaks is one iteration with noise, not three genuine attempts.

## Lessons Learned from Testing

### What Works
- **Phase separation** eliminates the "writing about a broken method" failure mode
- **Structured diagnosis** in Phase 2 produces actionable failure reports
- **File-based communication** prevents information degradation across agents
- **Scope constraints** prevent agents from "helping" in ways that blur boundaries

### What Doesn't Work (Yet)
- **Self-grading inflates by 1–4 points** compared to independent human grading
- **The system is gap-closing, not frontier-pushing** — it can produce a competent paper about a validated method, but it cannot push the methodological frontier beyond what the architect conceives
- **Context pressure** in Phase 3 (write-manuscript) is the tightest bottleneck — the orchestrator manages 4 sub-agents while tracking a complex decision log

### Cross-Domain Iteration Insights
These lessons from other fields informed the v2 design:

1. **Wiles heuristic** (mathematics): When stuck, analyze *why* it fails instead of trying harder
2. **AlphaFold 1→2** (ML): Proof of concept ≠ production system — sometimes throw away the architecture
3. **Kuhn phases** (philosophy of science): Normal science → anomaly accumulation → crisis. Different phases need different strategies
4. **Beauvais Cathedral** (engineering): Iterating past fundamental limits = catastrophic failure. Need kill criteria
5. **Double Diamond** (design): Separate problem-framing from solution-finding

## Usage

### As Claude Code Skills

1. Clone or copy the `.claude/skills/` directory into your project
2. Ensure each skill's `SKILL.md` is in its own subdirectory under `.claude/skills/`
3. Write a `research_brief.md` describing your problem (what to solve, not how)
4. Invoke: `/research-pipeline research_brief.md`

### Requirements

- Claude Code with access to the Agent/Task tools
- Python environment with numpy, scipy, pandas, matplotlib (for Phase 2–3 code execution)
- LaTeX installation (for Phase 3–4 compilation checks)

### Adapting to Your Domain

The skills are written for Bayesian clinical trial methods targeting *Biometrics*, but most of the architecture is domain-agnostic. To adapt:

- **`methodology-architect`**: Update the "adjacent fields" list in Phase 2.5 for your domain
- **`paper-writer`**: Change the journal conventions (document class, citation style, section structure)
- **`paper-grader`**: Adjust the rubric calibration notes for your target journal
- **`paper-critic`**: Update the review dimensions for your field's standards

The pipeline orchestrator, validate-method, paper-modeler, and paper-fixer are largely domain-independent.

## Repository Structure

```
ARMS/
├── README.md                              # This file
├── LICENSE                                # CC BY-NC 4.0
├── DISCLAIMER.md                          # LLM verification disclaimer
├── research_brief.md                      # Problem statement (method-agnostic)
├── skills/                                # 10 pipeline skills
│   ├── research-pipeline/SKILL.md         #   Outer orchestrator
│   ├── methodology-architect/SKILL.md     #   Phase 1: THINK
│   ├── validate-method/SKILL.md           #   Phase 2: VALIDATE
│   ├── write-manuscript/SKILL.md          #   Phase 3: WRITE
│   ├── literature-lead/SKILL.md           #   Phase 3 sub-agent
│   ├── paper-modeler/SKILL.md             #   Phase 3 sub-agent
│   ├── paper-writer/SKILL.md              #   Phase 3 sub-agent
│   ├── paper-critic/SKILL.md              #   Phase 3 sub-agent
│   ├── paper-grader/SKILL.md              #   Phase 4: POLISH
│   └── paper-fixer/SKILL.md               #   Phase 4: POLISH
├── data/                                  # Corrected NDMM dataset (35 arms, 21 trials)
├── code/                                  # Data curation + KG builder
├── case-studies/
│   ├── README.md                          # Comparison table + what to look at
│   ├── 1-fully-autonomous/                # KG-DAP: 3h, zero human input, 38/50
│   │   └── pipeline/                      #   Complete Phase 1-4 artifacts
│   └── 2-human-revised/                   # KG-CAR: ~3 weeks, 6h human, 39/50
│       ├── manuscript.tex                 #   Human-edited final
│       ├── manuscript-auto.tex            #   LLM baseline before editing
│       ├── code/                          #   Analysis scripts
│       └── results/                       #   All output CSVs
└── docs/                                  # Session summaries
```

## Case Studies

Two papers were produced from the **same research brief and corrected dataset**, using different workflows:

|  | Fully Autonomous | Human-Revised |
|--|-----------------|---------------|
| Method discovered | KG-DAP (Beta mixture) | KG-CAR (spatial CAR) |
| Time | 3h 14min | ~3 weeks |
| Human effort | 0 | ~6 hours |
| Independent score | 38/50 (B) | 39/50 (B) |
| Real-data advantage | 25.9% MAE over rMAP | 5.9% MSE over rMAP |

Neither method was prescribed in the brief. The pipeline independently discovered a different approach each time it ran. See [`case-studies/README.md`](case-studies/README.md) for detailed comparison.

## Provenance

Developed over sessions D-K (March 2026) as part of a project to automate production of Bayesian clinical trial methods papers. The v1 system (sessions D-E, simple write-grade-fix loop) identified the plateau problem; the v2 system (sessions J-K) was designed from scratch to address it. The autonomous case study (KG-DAP) was produced on 2026-03-23 in a single pipeline run.

## License

CC BY-NC 4.0 (Creative Commons Attribution-NonCommercial 4.0 International)
