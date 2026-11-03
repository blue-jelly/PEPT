# PEPT target-panel audit — strict single-assay results

Run: `python src/chembl_audit.py --out reports/target_panel_audit.csv`
Source: ChEMBL API, 3 Nov 2026. Filter: exact IC50 only
(`standard_type=IC50`, `standard_relation='='`, `standard_units='nM'`),
human targets. Scaffolds: Murcko, RDKit (acyclic compounds excluded
from scaffold counts, not counted as one '' scaffold).

| Target | ChEMBL | Records | Unique compounds | Assay families (single assay) | Best single assay | Best: compounds / scaffolds | Gate 300/50 |
|---|---|---|---|---|---|---|---|
| AR | CHEMBL1871 | 3,603 | 2,370 | 377 | CHEMBL5733828 | 53 / 28 | FAIL |
| hERG (KCNH2) | CHEMBL240 | 12,010 | 9,804 | 1,952 | CHEMBL1827362 | 217 / 118 | FAIL |
| ADRB2 | CHEMBL210 | 685 | 623 | 64 | CHEMBL1947791 | 50 / 26 | FAIL |

## What this does and does not mean

- **It does NOT mean there is no hERG or AR data.** Overall counts are
  large (hERG 9,804 compounds). The failure is specific: no *single*
  assay reaches 300 compounds. hERG is closest (217, with 118 scaffolds —
  scaffold diversity is not the blocker, compound count in one assay is).
- **It does NOT justify pooling unrelated assays to reach 300.**
  Pooling all patch-clamp hERG assays gives 3,366 compounds, but that
  mixes binding, displacement, manual and automated patch protocols —
  that is not a coherent family, and it is not claimed as a pass here.
- **The export's AR HEK293 pooled cohort (4 assays, 156 compounds /
  55 scaffolds) is consistent with this audit:** it is below 300 too,
  under either definition. The earlier pilot result (Random Forest not
  beating the mean on that cohort) therefore stands as a pilot, not a model.
- **Consequence for V1:** `src/models.py` correctly refuses to run for
  AR, hERG and ADRB2, because their gate status is FAIL. No QSAR
  metrics are claimed in this version.

## Fallback panel — same gate, also FAIL

| Target | ChEMBL | Records | Unique compounds | Best single assay | Best: compounds / scaffolds | Gate |
|---|---|---|---|---|---|---|
| DRD2 | CHEMBL217 | 1,100 | 891 | CHEMBL5732166 | 94 / 53 | FAIL |
| 5-HT2A | CHEMBL224 | 1,432 | 1,215 | CHEMBL1827733 | 112 / 8 | FAIL |

DRD2 passes the scaffold part (53) but not compounds (94).
5-HT2A's best assay is a narrow analogue series (112 compounds, 8 scaffolds).

## Next step, pre-registered in the export

Test the fallback panel DRD2 → 5-HT2A with the *same* gate, rather than
lowering the gate for AR/hERG/ADRB2. Done above: neither passes either.
So V1 is an audit + repo, and the model layer waits for a target that
actually passes — that is a publishable negative feasibility result,
not a failure to try hard enough.

Separately, and clearly labelled as NOT gate-passing, an exploratory
pooled hERG model (`src/pooled_model.py`) is provided: all 9,804 hERG
compounds pooled across assays, scaffold-grouped split. Present it only
with its assay-heterogeneity limitation, alongside this audit.

## Pooled hERG exploratory model — results (3 Nov 2026)

9,802 compounds modelled (2 failed featurisation), 5,132 scaffolds,
scaffold-grouped split: train 7,925 / test 1,877, no scaffold in both.

| Model | RMSE (pIC50) | MAE | R² |
|---|---|---|---|
| Dummy mean | 0.870 | 0.640 | 0.00 |
| Ridge | 0.757 | 0.583 | 0.243 |
| Random Forest | 0.619 | 0.449 | 0.492 |

Interpretation: RF beats the mean baseline substantially on unseen
scaffolds, but R² 0.49 is not a production / safety model, and the
pooled assay heterogeneity means these numbers are an upper-bound
pipeline demo, not a gate-passing result.
