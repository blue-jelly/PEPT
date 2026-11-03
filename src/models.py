#!/usr/bin/env python3
"""
Baseline QSAR ladder — ONLY for targets that PASS reports/target_panel_audit.csv.

Refuses to run for a FAIL/UNKNOWN/NO_AUDIT target. That refusal is intentional:
the export's AR pilot (156 compounds / 55 scaffolds) is below the gate,
so it must not be presented as a trained model.

Pipeline (for a passing target):
  raw ChEMBL IC50 JSON -> one row per compound (median pIC50, pIC50 = -log10(IC50_M))
  -> Morgan fingerprint (2048, r=2) + 7 RDKit descriptors (src/features.py)
  -> scaffold-grouped test split (Murcko scaffolds never shared train/test)
  -> ladder: dummy mean, Ridge, RandomForest, XGBoost (if installed)
  -> metrics: RMSE / MAE / R^2 on the scaffold-held-out test set

Usage:
  python src/models.py --target CHEMBL240 --raw data/raw/hERG_CHEMBL240_ic50.json \
      --out reports/hERG_qsar.csv
"""
import argparse, csv, json, sys
import numpy as np

def gate_status(target_id):
    try:
        with open("reports/target_panel_audit.csv") as f:
            for row in csv.DictReader(f):
                if row["target_chembl_id"] == target_id:
                    return row["gate_300_compounds_50_scaffolds"]
    except FileNotFoundError:
        return "NO_AUDIT"
    return "NOT_FOUND"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    status = gate_status(args.target)
    if status != "PASS":
        print(f"REFUSED: target {args.target} gate status is {status}, not PASS. "
              "Run src/chembl_audit.py first and only model passing targets.", file=sys.stderr)
        return 2

    from collections import defaultdict
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

    # One row per compound, median pIC50. Skip non-numeric / non-positive values.
    by_compound = defaultdict(list)
    smiles = {}
    for a in acts:
        cid = a.get("molecule_chembl_id")
        try:
            v = float(a.get("standard_value"))
        except (TypeError, ValueError):
            continue
        if not v or v <= 0:
            continue
        # standard_units is nM by the audit filter; pIC50 = -log10(v * 1e-9)
        by_compound[cid].append(-np.log10(v * 1e-9))
        smiles.setdefault(cid, a.get("canonical_smiles"))

    X, y, groups = [], [], []
    for cid, vals in by_compound.items():
        feat = featurize(smiles.get(cid))
        if feat is None:
            continue
        mol = Chem.MolFromSmiles(smiles[cid])
        try:
            scaf_mol = MurckoScaffold.GetScaffoldForMol(mol)
            scaf = Chem.MolToSmiles(scaf_mol) if (scaf_mol is not None and scaf_mol.GetNumAtoms() > 0) else "NO_SCAFFOLD"
        except Exception:
            scaf = "NO_SCAFFOLD"
        X.append(feat); y.append(float(np.median(vals))); groups.append(scaf)
    X = np.array(X); y = np.array(y); groups = np.array(groups)
    print(f"Compounds modelled: {len(y)}  scaffolds: {len(set(groups))}")

    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    tr, te = next(gss.split(X, y, groups))
    Xtr, Xte, ytr, yte = X[tr], X[te], y[tr], y[te]
    assert not set(groups[tr]) & set(groups[te]), "scaffold leak in split"

    models = [
        ("dummy_mean", DummyRegressor(strategy="mean")),
        ("ridge", Ridge(alpha=1.0)),
        ("random_forest", RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=42)),
    ]
    try:
        from xgboost import XGBRegressor
        models.append(("xgboost", XGBRegressor(n_estimators=300, max_depth=6, learning_rate=0.05, n_jobs=-1, random_state=42)))
    except Exception as e:
        print(f"xgboost unavailable, skipping: {e}")

    rows = []
    for name, m in models:
        m.fit(Xtr, ytr)
        pred = m.predict(Xte)
        row = {
            "target": args.target, "model": name, "n_train": len(ytr), "n_test": len(yte),
            "rmse": float(np.sqrt(mean_squared_error(yte, pred))),
            "mae": float(mean_absolute_error(yte, pred)),
            "r2": float(r2_score(yte, pred)),
        }
        rows.append(row)
        print(row)

    out = args.out or f"reports/{args.target}_qsar.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"Wrote {out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
