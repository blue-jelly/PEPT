#!/usr/bin/env python3
"""
Exploratory pooled target-level QSAR — NOT a gate-passing model.

This is deliberately separate from src/models.py, which refuses to run
for gate-FAIL targets. This script trains on ALL exact human IC50 records
for one target, pooled across assays, with a scaffold-grouped test split
(so no Murcko scaffold appears in both train and test).

Limitations, stated in the output and in the report:
  - Assay heterogeneity is NOT removed. Different protocols, cell systems
    and readouts are pooled, so pIC50 values are not strictly comparable.
  - This must NOT be presented as "passing the 300/50 single-family gate".
  - It is useful as: (a) a working end-to-end pipeline demo, (b) a baseline
    for a future coherent-family dataset, (c) a check that the pipeline
    runs at target scale (e.g. hERG, 9,804 compounds).
  - If you present this for an internship, present it alongside
    reports/target_panel_audit_summary.md, with this limitation on the slide.

Usage:
  python src/pooled_model.py --raw data/raw/hERG_CHEMBL240_ic50.json --out reports/hERG_pooled_qsar.csv
"""
import argparse, csv, json
from collections import defaultdict
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from sklearn.model_selection import GroupShuffleSplit
    from sklearn.dummy import DummyRegressor
    from sklearn.linear_model import Ridge
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
    from features import featurize

    with open(args.raw) as f:
        acts = json.load(f)
    by_compound = defaultdict(list); smiles = {}
    for a in acts:
        cid = a.get("molecule_chembl_id")
        try: v = float(a.get("standard_value"))
        except (TypeError, ValueError): continue
        if not v or v <= 0: continue
        by_compound[cid].append(-np.log10(v * 1e-9))
        smiles.setdefault(cid, a.get("canonical_smiles"))

    X, y, groups = [], [], []
    for cid, vals in by_compound.items():
        feat = featurize(smiles.get(cid))
        if feat is None: continue
        mol = Chem.MolFromSmiles(smiles[cid])
        try:
            scaf_mol = MurckoScaffold.GetScaffoldForMol(mol)
            scaf = Chem.MolToSmiles(scaf_mol) if (scaf_mol is not None and scaf_mol.GetNumAtoms() > 0) else "NO_SCAFFOLD"
        except Exception:
            scaf = "NO_SCAFFOLD"
        X.append(feat); y.append(float(np.median(vals))); groups.append(scaf)
    X = np.array(X); y = np.array(y); groups = np.array(groups)
    print(f"Compounds: {len(y)}  scaffolds: {len(set(groups))}  (pooled, assay-heterogeneous)", flush=True)

    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    tr, te = next(gss.split(X, y, groups))
    Xtr, Xte, ytr, yte = X[tr], X[te], y[tr], y[te]
    assert not set(groups[tr]) & set(groups[te])

    models = [
        ("dummy_mean", DummyRegressor(strategy="mean")),
        ("ridge", Ridge(alpha=1.0)),
        ("random_forest", RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=42)),
    ]
    rows = []
    for name, m in models:
        m.fit(Xtr, ytr); pred = m.predict(Xte)
        row = {"model": name, "n_train": len(ytr), "n_test": len(yte),
               "rmse": float(np.sqrt(mean_squared_error(yte, pred))),
               "mae": float(mean_absolute_error(yte, pred)),
               "r2": float(r2_score(yte, pred))}
        rows.append(row); print(row, flush=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"Wrote {args.out} — POOLED, NOT GATE-PASSING. See module docstring.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
