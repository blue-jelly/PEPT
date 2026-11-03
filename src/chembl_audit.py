#!/usr/bin/env python3
"""
PEPT target-panel feasibility audit.

Queries the public ChEMBL API for human targets, exact IC50 only,
and tests the pre-registered gate:
  one coherent IC50 assay family -> >=300 unique compounds AND >=50 scaffolds

Assay family here = (assay_chembl_id, assay_type, assay_description, cell_type, document)
grouped, so a count pooled across unrelated assays cannot pass the gate.

Scaffolds: Murcko scaffolds via RDKit if available; otherwise the script
still reports compound/assay counts and marks scaffold counts as NA
instead of inventing them. Do NOT present an NA scaffold count as a pass.

Targets (verified 3 Nov 2026 against ChEMBL API):
  AR    CHEMBL1871  Androgen receptor, Homo sapiens
  hERG  CHEMBL240   Voltage-gated inwardly rectifying potassium channel KCNH2
  ADRB2 CHEMBL210   Beta-2 adrenergic receptor, Homo sapiens

Usage:
  python src/chembl_audit.py --out reports/target_panel_audit.csv [--limit-per-target 0]
  --limit-per-target 0 = no limit (default). Use a small limit only for smoke tests.
"""
import argparse
import csv
import json
import sys
import time
import urllib.request
import urllib.parse
from collections import Counter, defaultdict

BASE = "https://www.ebi.ac.uk/chembl/api/data"

TARGETS = {
    "AR": "CHEMBL1871",
    "hERG": "CHEMBL240",
    "ADRB2": "CHEMBL210",
}

try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False


def get_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def fetch_activities(target_id, limit=0):
    """Fetch exact human IC50 activities for a target, paginated."""
    rows = []
    offset = 0
    page = 1000
    while True:
        q = {
            "target_chembl_id": target_id,
            "standard_type": "IC50",
            "standard_relation": "=",
            "standard_units": "nM",
            "limit": page,
            "offset": offset,
        }
        url = f"{BASE}/activity.json?{urllib.parse.urlencode(q)}"
        data = get_json(url)
        acts = data.get("activities", [])
        if not acts:
            break
        rows.extend(acts)
        if limit and len(rows) >= limit:
            rows = rows[:limit]
            break
        meta = data.get("page_meta", {})
        if meta.get("next") is None:
            break
        offset += page
        time.sleep(0.1)
    return rows


def murcko(smiles):
    if not HAS_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        # MurckoScaffoldSmiles() takes a SMILES string in this RDKit build;
        # the Mol-based API is GetScaffoldForMol().
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        if scaf is None or scaf.GetNumAtoms() == 0:
            return None  # acyclic: no Murcko scaffold, do not count '' as one
        return Chem.MolToSmiles(scaf)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="reports/target_panel_audit.csv")
    ap.add_argument("--limit-per-target", type=int, default=0)
    ap.add_argument("--raw-dir", default="data/raw")
    args = ap.parse_args()

    import os
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    os.makedirs(args.raw_dir, exist_ok=True)

    summary = []
    for name, tid in TARGETS.items():
        print(f"[{name} {tid}] querying ChEMBL exact human IC50 ...", flush=True)
        acts = fetch_activities(tid, limit=args.limit_per_target)
        # Keep only human assays if the API returned any non-human rows.
        # ChEMBL activity rows carry assay organism via the assay endpoint;
        # the activity filter above is target-based, so record the filter
        # actually applied in the report rather than assuming organism.
        with open(os.path.join(args.raw_dir, f"{name}_{tid}_ic50.json"), "w") as f:
            json.dump(acts, f)

        compounds = set()
        smiles_by_compound = {}
        for a in acts:
            cid = a.get("molecule_chembl_id")
            if cid:
                compounds.add(cid)
                if cid not in smiles_by_compound:
                    smiles_by_compound[cid] = a.get("canonical_smiles")

        # Assay-family grouping
        families = defaultdict(list)
        for a in acts:
            key = (
                a.get("assay_chembl_id"),
                a.get("assay_type"),
            )
            families[key].append(a)

        family_stats = []
        for (assay_id, assay_type), fam in families.items():
            fam_compounds = {a.get("molecule_chembl_id") for a in fam if a.get("molecule_chembl_id")}
            fam_scaffolds = set()
            if HAS_RDKIT:
                for cid in fam_compounds:
                    s = murcko(smiles_by_compound.get(cid))
                    if s:
                        fam_scaffolds.add(s)
            family_stats.append({
                "assay_chembl_id": assay_id,
                "assay_type": assay_type,
                "n_records": len(fam),
                "n_compounds": len(fam_compounds),
                "n_scaffolds": len(fam_scaffolds) if HAS_RDKIT else None,
            })
        family_stats.sort(key=lambda x: x["n_compounds"], reverse=True)
        best = family_stats[0] if family_stats else {"n_compounds": 0, "n_scaffolds": None}
        passes = (
            HAS_RDKIT
            and best["n_compounds"] >= 300
            and (best["n_scaffolds"] or 0) >= 50
        )
        row = {
            "target": name,
            "target_chembl_id": tid,
            "n_exact_ic50_records": len(acts),
            "n_unique_compounds_overall": len(compounds),
            "n_assay_families": len(family_stats),
            "best_family_assay_chembl_id": best.get("assay_chembl_id"),
            "best_family_assay_type": best.get("assay_type"),
            "best_family_n_records": best.get("n_records", 0),
            "best_family_n_compounds": best["n_compounds"],
            "best_family_n_scaffolds": best["n_scaffolds"],
            "gate_300_compounds_50_scaffolds": "PASS" if passes else ("UNKNOWN_RDKIT_MISSING" if not HAS_RDKIT else "FAIL"),
            "rdkit_available": HAS_RDKIT,
        }
        summary.append(row)
        print(f"  records={len(acts)} compounds={len(compounds)} "
              f"best_family={best['n_compounds']} compounds, "
              f"{best['n_scaffolds']} scaffolds -> {row['gate_300_compounds_50_scaffolds']}", flush=True)
        with open(os.path.join("reports", f"{name}_assay_families.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(family_stats[0].keys()) if family_stats else ["assay_chembl_id"])
            w.writeheader()
            w.writerows(family_stats)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    print(f"Wrote {args.out}")
    if not HAS_RDKIT:
        print("WARNING: RDKit not available — scaffold gate is UNKNOWN, not PASS. "
              "Install via environment.yml (conda/mamba) or Codespaces devcontainer.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
