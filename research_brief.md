# Research Brief: Leveraging Structural Knowledge for Clinical Trial Inference

## Problem

In clinical trials — especially for rare diseases or early-phase studies — sample
sizes are often too small for reliable inference about treatment efficacy. A pool of
completed historical trials in the same therapeutic area exists and contains relevant
information, but it is unclear how to use it responsibly.

The core challenge: **how much information to borrow from which historical trials.**

Borrow too much from an irrelevant trial → biased estimates, inflated confidence.
Borrow too little → waste available information, no improvement over ignoring history.
Borrow uniformly → treat all historical trials as equally informative, ignoring known
structural differences between them.

We have a **biomedical knowledge graph** that encodes structural relationships between
trials — which drugs they use, which molecular targets those drugs hit, which patient
populations they enroll, which endpoints they measure. This structural information
could inform borrowing decisions, but no existing method fully exploits it.

## Research Question

How can structural knowledge encoded in a biomedical knowledge graph be used to
improve inference for a new clinical trial by selectively borrowing from historical
trials — in a way that is principled, robust, computationally tractable, and
compatible with regulatory requirements?

## Abstract Problem (for cross-disciplinary thinking)

Strip away the clinical trial specifics and the core problem is:

> **Given a target task and a pool of related source tasks with known structural
> relationships, how do you optimally combine information from sources to improve
> estimation on the target — while being robust to sources that turn out to be
> irrelevant?**

This abstract problem appears across many fields:
- Multi-source transfer learning and domain adaptation
- Sensor fusion with heterogeneous reliability
- Expert opinion aggregation and forecast combination
- Meta-learning and few-shot learning
- Causal transportability across populations
- Bayesian model averaging and stacking

The methodology architect should actively search OUTSIDE the clinical trial
statistics literature for ideas. The best solution may combine a borrowing
mechanism from one field with an uncertainty framework from another.

## Available Data

### Clinical Trial Data
- **35 trial arms from 21 newly diagnosed multiple myeloma (NDMM) studies** reporting
  MRD (minimal residual disease) negativity rates, in `data/trials_data.csv`
  - Columns: trial name, arm, sample size (n), responders (y), observed rate (phat),
    regimen details, transplant status
  - Response rates range from ~0.07 to ~0.79 across arms
- **Patient-level data** for one trial (Elo-KRd, 30 patients, 28 evaluable)
- **Simulated patient-level data** for MASTER (n=123) and MANHATTAN (n=41) trials
- **External trials for robustness testing**:
  - 10 adjacent MM trials (different endpoints: CR, ORR) in `data/external_adjacent_trials.csv`
  - 10 foreign indication trials (breast cancer, CLL, AML) in `data/external_foreign_trials.csv`

### Knowledge Graph
- **56-node NetworkX MultiDiGraph** in `data/kg_graph.gpickle`
  - Node types: drug, molecular target, trial, endpoint, population
  - Edge types: drug-targets-protein, trial-uses-drug, trial-measures-endpoint,
    trial-enrolls-population, etc.
  - 11 drugs, 6 molecular targets
- **Drug-target map** in `data/drug_target_map.csv`
- Trial-to-trial structural similarity can be computed from the KG based on shared
  drugs, shared molecular targets, and population characteristics

### Data Infrastructure
- `data_curation.py` — data loading, cleaning, patient-level simulation
- `build_kg.py` — knowledge graph construction, similarity matrix computation

## Target Journal

**Biometrics**. Expects:
- Rigorous theory (propositions/theorems with proofs or proof sketches)
- Comprehensive evaluation (simulation + real data + sensitivity/robustness analysis)
- Comparison to established methods in the field
- Clear practical value for the target audience (biostatisticians, clinical trialists)
- Typical length: 20-30 pages double-spaced

## Reference Papers

Papers are available in `reference/`. The methodology architect should also conduct
its own literature search — the papers on disk are starting points, not the complete
landscape.

## Success Criteria

- **Estimation accuracy**: Measurable improvement over a non-informative baseline,
  especially at small sample sizes where borrowing matters most
- **Calibration**: Credible intervals should have appropriate coverage (neither
  over-confident nor overly conservative)
- **Robustness**: Adding irrelevant historical trials to the pool should cause
  minimal degradation — ideally less than existing borrowing methods
- **Computational tractability**: Fast enough for interactive use and hyperparameter
  exploration (seconds, not minutes)
- **Pre-specifiability**: Ideally, the borrowing mechanism can be fully determined
  before any new trial data is observed — this matters for regulatory statistical
  analysis plans
- **Interpretability**: Practitioners should understand why the method borrows more
  from some trials than others

## What We Are NOT Prescribing

This brief intentionally does not specify:
- Which statistical framework to use (Bayesian, frequentist, or hybrid)
- Which specific borrowing mechanism to employ
- Which existing methods to build upon
- What the mathematical form of the solution should look like

The methodology architect should reason about the problem from first principles,
informed by a broad literature search across multiple disciplines. The best ideas
often come from importing a concept from a different field rather than incrementally
refining within the home field.
