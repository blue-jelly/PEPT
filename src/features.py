"""RDKit features for PEPT: Morgan fingerprints + descriptors."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def featurize(smiles, n_bits=2048, radius=2):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    fp_arr = np.zeros((n_bits,), dtype=np.float32)
    from rdkit import DataStructs
    DataStructs.ConvertToNumpyArray(fp, fp_arr)
    desc = np.array([
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol),
        rdMolDescriptors.CalcNumRings(mol),
    ], dtype=np.float32)
    return np.concatenate([fp_arr, desc])
