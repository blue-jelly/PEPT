# PEPT — Performance Pharmacology & Endophenotype Translator

**Final version scope (locked 3 Nov 2026, from the Sep–13 Oct 2026 ChatGPT export)**

A reproducible, research-style computational pharmacology project for a
targeted research internship — **not** a dosing, ranking, or
"best compound" tool. Compounds such as yohimbine, clenbuterol,
trenbolone and GW501516 appear only as reference perturbations in public
data, where the dataset contains them.

## Research question

Can molecular structure predict activity at performance-relevant targets,
and can those predictions be translated — honestly, with assay-aware
filtering — into a mechanistic picture of tissue signalling?

## What the final version actually is

1. **Target-panel feasibility audit (the gate, before any modelling)**
   - Targets: AR (`CHEMBL1871`), hERG / KCNH2 (`CHEMBL240`),
     ADRB2 (`CHEMBL210`)
   - Data: ChEMBL, human only, exact IC50 only (`standard_relation = '='`,
     `standard_type = 'IC50'`, `standard_units = 'nM'`)
   - A target passes only if **one scientifically coherent IC50 assay
     family** yields **≥300 unique compounds and ≥50 Murcko scaffolds**.
     If ADRB2 fails, the pre-registered fallback is DRD2 → 5-HT2A —
     the gate is not massaged to force a pass.
   - **Family definition, stated plainly:** `chembl_audit.py` first reports
     the strict lower bound — a single `assay_chembl_id`. The export's
     intended family is broader and pooled, e.g. the AR HEK293 WT
     luciferase cohort (4 assays, 163 measurements, 156 compounds).
     Both levels are reported in the final report; a target is only
     claimed to pass on a pooled family if that family's assays share
     cell system, endpoint and readout, with the assay IDs listed —
     never by pooling unrelated assays to reach 300.
   - This matters because the numbers on file are misleading in isolation:
     AR has 3,599 exact human IC50 records overall, but the preliminary
     modelling cohort in the export was only 156 compounds / 55 scaffolds
     (below the gate). hERG has 12,010 exact human IC50 records overall.
     Strict AR–hERG IC50×IC50 overlap was only 50 shared / 14 exact paired
     compounds — so no selectivity claim is built on overlap alone.

2. **Baseline QSAR — only for targets that pass the gate**
   - Morgan fingerprints (2048-bit, radius 2) + RDKit descriptors
   - Model ladder: majority → logistic regression → random forest →
     XGBoost, compared on scaffold-grouped splits, not random splits
   - Metrics: RMSE / MAE / R² for pIC50 regression; ROC-AUC / PR-AUC /
     F1 / MCC if a classification cut is ever used
   - Interpretability: feature importance + SHAP, mapped back to
     molecular substructures where RDKit allows it

3. **What is deliberately NOT in this repo**
   - No optimal dose, no "best steroid / fat-burner" ranking, no
     human-performance prediction from a molecule, no therapeutic-index
     claim. Those need exposure, tissue, dose and physiology data this
     project does not have.

## Results so far (3 Nov 2026, ChEMBL live audit)

Strict single-assay gate (≥300 compounds / ≥50 scaffolds): AR, hERG,
ADRB2, DRD2 and 5-HT2A all FAIL — see
`reports/target_panel_audit_summary.md` for the full table.

Exploratory pooled hERG model (`src/pooled_model.py`, NOT gate-passing,
assay-heterogeneous, scaffold-grouped test, 7,925 train / 1,877 test):
Random Forest RMSE 0.619, MAE 0.449, R² 0.492, vs dummy mean R² ~0.
Present this only alongside the audit and its limitation.

## Repo layout

```
pept/
├── README.md
├── requirements.txt
├── environment.yml          # conda/mamba (recommended for RDKit)
├── .devcontainer/           # GitHub Codespaces (Student Pack)
├── .github/workflows/ci.yml # re-runs the audit on push
├── src/
│   ├── chembl_audit.py      # Step 1: target-panel feasibility audit
│   ├── features.py          # RDKit fingerprints + descriptors
│   └── models.py            # baseline QSAR ladder + scaffold split
├── data/
│   ├── raw/                 # ChEMBL downloads (not committed if large)
│   ├── interim/
│   └── processed/
├── reports/                 # audit tables (CSV + summary)
├── figures/
├── notebooks/               # exploratory only; scripts are canonical
└── paper/                   # manuscript draft, when results warrant it
```

## Quick start

```bash
# Recommended: conda/mamba, because RDKit via plain pip is platform-sensitive
mamba env create -f environment.yml
mamba activate pept

# Or, if you have the GitHub Student Pack: open this repo in Codespaces
# — the .devcontainer builds the same environment, no local RDKit pain.

# Step 1 — the audit (no modelling until a target passes the gate)
python src/chembl_audit.py --out reports/target_panel_audit.csv

# Step 2 — only for passing targets, e.g. hERG if it passes
python src/models.py --target CHEMBL240 --scaffold-split
```

## GitHub Student Pack — what is used here, and what is not needed

- **Used:** a GitHub repo (private while drafting, public for the
  internship application), GitHub Actions for the audit CI, and
  Codespaces / devcontainer as the no-install fallback for RDKit.
- **Not required:** Copilot, cloud credits, or any paid compute.
  The ChEMBL audit is small enough to run locally or in Actions.
  Nothing in this repo requires a Student Pack subscription to reproduce
  — it just removes the Windows RDKit-install friction.

## Provenance

- ChEMBL (EBI) — compound / target / bioactivity data, via its public API.
- Tox21 (NCATS) — independent assay validation in later versions only,
  not mixed into the ChEMBL IC50 gate.
- LINCS/L1000, Open Targets — V2+ transcriptomics / phenotype layers,
  deliberately **not** in this final V1 audit. They are roadmap, not results.

## Reproducibility note

Every number in `reports/` must be regenerable by re-running
`src/chembl_audit.py`. If a number in this README cannot be regenerated
that way, treat it as a draft claim from the earlier export, not a result
of this repo.
