# Data Correction Plan

**Date**: 2026-03-23
**Source**: Verification of all 23 NDMM trial sources by web search (3 parallel agents)
**Status**: PLAN — not yet executed

---

## Summary

- **8 trials fully verified** (no changes needed): PERSEUS, CASSIOPEIA, MAIA, ALCYONE, GMMG-HD7, BENEFIT, FORTE, MASTER
- **5 trials with wrong journal citation only** (data correct): CEPHEUS, GRIFFIN, MANHATTAN, IsKia, Elo-KRd
- **3 trials with wrong MRD data**: ADVANCE, ENDURANCE, DETERMINATION
- **2 trials with structural issues**: MIDAS (conditional rates), COBRA (wrong trial entirely)
- **2 trials with wrong rates + source**: IFM2018-04, GMMG-CONCEPT
- **2 trials to DROP**: ELOQUENT-1, SWOG-S1211 (MRD never reported)

**Final dataset**: 35 arms from 21 trials (was 40 arms from 23 trials)

---

## Step 1: DROP These Trials

### ELOQUENT-1 (NCT01335399) — 2 arms
**Reason**: MRD negativity was never a reported endpoint. The p=0.150/0.099 values cannot be verified.
**Correct publication**: Dimopoulos et al. Lancet Haematol 2022; 9: e403-e414.
**Action**: Remove `eloquent1_elord` and `eloquent1_rd` from `data_curation.py`.

### SWOG-S1211 (NCT01668719) — 2 arms
**Reason**: MRD negativity was never reported. The p=0.558/0.529 values are fabricated.
**Correct publication**: Usmani et al. Lancet Haematol 2021; 8(1): e45-e54 (PMID: 33357482).
**Action**: Remove `swog_s1211_elovrd` and `swog_s1211_vrd` from `data_curation.py`.

---

## Step 2: FIX Source Citations Only (data correct)

### Elo-KRd
- **Current**: "Usmani et al. Blood Adv 2024"
- **Correct**: "Derman et al. JAMA Oncol 2022" (PMID: 35862034)

### CEPHEUS
- **Current**: "Facon et al. NEJM 2024"
- **Correct**: "Facon et al. Nat Med 2025" (doi: 10.1038/s41591-024-03485-7)

### GRIFFIN
- **Current**: "Voorhees et al. Blood 2023"
- **Correct**: "Voorhees et al. Lancet Haematol 2023" (PMID: 37708911)

### MANHATTAN
- **Current**: "Landgren et al. Lancet Haematol 2021"
- **Correct**: "Landgren et al. JAMA Oncol 2021" (PMID: not available; doi: 10.1001/jamaoncol.2021.0611)

### IsKia
- **Current**: "Gay et al. NEJM 2024"
- **Correct**: "Gay et al. Blood 2023; 142(Suppl 1): 4" (ASH 2023 abstract)

---

## Step 3: FIX Data Errors

### ADVANCE (NCT04268498) — 2 arms
**Current**: source "Kumar et al. NEJM 2024", n=176/174, p=0.801/0.655
**Correct**:
- Source: "Landgren et al. JCO 2025; 43(suppl 16): 7503" (ASCO 2025 abstract)
- DKRd arm: n=148, y=87 (0.59 × 148 ≈ 87), p=0.588
- KRd arm: n=139, y=46 (0.33 × 139 ≈ 46), p=0.331
- Note: Exact y values should be verified from abstract; 59% and 33% are ITT rates

### DETERMINATION (NCT01208662) — 2 arms
**Current**: n=357 assigned to ASCT, n=365 assigned to VRd-alone
**Correct**: **Arms are swapped**. n=365 is ASCT, n=357 is VRd-alone.
- Source: Richardson et al. NEJM 2022 (correct, keep)
- VRd+ASCT arm: n=365, y=199 (0.544 × 365 ≈ 199), p=0.545
- VRd-alone arm: n=357, y=142 (0.398 × 357 ≈ 142), p=0.398
- Note: MRD rates (54.4%, 39.8%) from substudy of 198 patients, extrapolated to ITT

### ENDURANCE (NCT01863550) — 2 arms
**Current**: n=526/535, p=0.340/0.351 (y=179/188)
**Correct**: MRD rates are ~10%/7%, NOT 34%/35%. The current values are catastrophically wrong (likely confused with VGPR rates).
- Source: Kumar et al. Lancet Oncol 2020 (correct, keep)
- KRd arm: n=526, y=54 (0.103 × 526 ≈ 54), p=0.103
- VRd arm: n=527 (not 535), y=38 (0.072 × 527 ≈ 38), p=0.072
- Note: MRD threshold is unspecified in publication (flow cytometry, NOT NGS). Consider setting threshold to "unspecified" or "1e-4".

### IFM2009 (NCT01191060) — 2 arms
**Current**: ASCT arm y=228 (p=0.651), no-ASCT arm y=171 (p=0.489)
**Correct**: Minor correction to ASCT arm.
- ASCT arm: y=220 (not 228), p=0.629. Source: Attal et al. NEJM 2017, Table 2 (220/278 among CR/VGPR, divided by ITT n=350)
- no-ASCT arm: y=171 (p=0.489) — correct as-is
- Note: These are flow cytometry at 10^-4, not NGS at 10^-5. Keep mrd_threshold as "1e-4".

### IFM2018-04 (NCT03606577) — 1 arm
**Current**: source "Corre et al. Blood 2023", p=0.800
**Correct**:
- Source: "Touzeau et al. Blood 2024; 143(20): 2029-2036" (PMID: 38394666)
- n=50 (correct)
- Post-induction MRD-neg at 10^-5: 62% (23/37 evaluable). Using n=50 ITT: y=31, p=0.62
- Note: The 94% rate is at 10^-6 premaintenance (after tandem ASCT), which is a different threshold and timepoint.

### GMMG-CONCEPT (NCT03104842) — 1 arm
**Current**: source "Leypoldt et al. Lancet Haematol 2024", n=50, p=0.680
**Correct** (use the full JCO 2024 analysis):
- Source: "Leypoldt et al. JCO 2024; 42(3): 283-293" (PMID: 37753960)
- n=93 (evaluable TE patients, post-consolidation), y=63, p=0.677
- Alternative: n=99 (ITT TE), y=67 (0.677 × 99 ≈ 67), p=0.677
- Note: The old n=50 matched the Leukemia 2022 interim. The JCO 2024 full analysis is more reliable.

---

## Step 4: RESTRUCTURE These Trials

### MIDAS (NCT04934475) — merge 2 arms → 1 arm
**Current**: Two arms (Isa-KRd n=242 p=0.860, ASCT n=243 p=0.840)
**Problem**: The 86%/84% rates are CONDITIONAL (10^-6 given already 10^-5 negative) from Perrot et al. NEJM 2025. The randomization occurs AFTER MRD assessment — all patients received the same Isa-KRd induction.
**Correct**: Single-arm induction entry:
- trial_id: "midas_isakrd"
- arm_label: "Isa-KRd (induction)"
- n=791 (ITT) or n=757 (completed 6 cycles)
- p=0.63 (unconditional MRD-neg at 10^-5 after induction)
- y=498 (if n=791) or y=477 (if n=757)
- Source: "Perrot et al. Blood 2025; 146(1): 52-61"
- Remove the second arm (midas_asct)

### COBRA → Derman-DKRd (NCT03500445) — 1 arm, fix everything
**Current**: trial_id "cobra_dakrd", NCT04400500, n=100, p=0.800
**Problem**: "COBRA" is actually a different trial (KRd vs VRd, NCT03729804, Dytfeld et al.). The Derman Dara-KRd trial has a different NCT number and different data.
**Correct**:
- trial_id: "derman_dakrd"
- trial_name: "Derman-DKRd"
- nct_number: "NCT03500445"
- arm_label: "Dara-KRd"
- n=42 (enrolled) or n=40 (evaluable)
- MRD-neg at 10^-5: ~63% (20/32 with evaluable MRD samples). Using n=42: y=26, p=0.619
- Source: "Derman et al. Blood Cancer J 2024; 14(1): 87" (PMID: 38811560)
- population: mixed (transplant-ineligible intent, UChicago)

---

## Step 5: Execution Checklist

After applying all corrections:

- [ ] Update `data_curation.py` with corrected entries
- [ ] Re-run `python data_curation.py` to regenerate `trials_data.csv`
- [ ] Re-run `python build_kg.py` to regenerate similarity matrices (new arm set)
- [ ] Verify: `trials_data.csv` has 35 rows (header + 35 data rows)
- [ ] Verify: similarity matrix is 35×35
- [ ] Re-run LOO-CV: `run_production.py` or equivalent
- [ ] Re-run simulation study: `simulation_study.py`
- [ ] Re-run ESS analysis: `compute_elir_ess.py`
- [ ] Re-run design OC: `design_oc_structured.py`
- [ ] Re-run sensitivity analysis
- [ ] Re-run contamination robustness
- [ ] Re-run component contribution analysis
- [ ] Regenerate all figures
- [ ] Update all tables in `manuscript-manual-edit.tex`
- [ ] Update abstract (numbers will change)
- [ ] Update Section 1.2 (40 arms → 35 arms, 23 trials → 21 trials)
- [ ] Update Supplementary Table S5 (trial arms table)
- [ ] Regenerate similarity heatmap (35×35)
- [ ] Update all text references to "40 arms" / "23 studies"
- [ ] Verify LaTeX compiles cleanly

---

## Expected Impact on Results

The corrected dataset has:
- **Lower average MRD rate** (ADVANCE drops from 80%/65% to 59%/33%; ENDURANCE drops from 34%/35% to 10%/7%)
- **More heterogeneity** (range will be wider: ~7% to ~86%)
- **Fewer arms** (35 vs 40) but all verifiable
- The KG similarity structure is preserved — drug/target overlaps don't change for retained trials
- KG-CAR's advantage should persist or increase (more heterogeneity → more value from differential borrowing)
