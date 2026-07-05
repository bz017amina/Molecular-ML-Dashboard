import os
import joblib
import base64
from io import BytesIO
import numpy as np
from flask import Flask, request, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, AllChem, Draw

app = Flask(__name__)
CORS(app)

# 1. Chemins des fichiers du modèle
MODEL_LOGP_PATH = "models/logP_best_model.pkl"
MODEL_QED_PATH = "models/qed_best_model.pkl"
MODEL_SAS_PATH = "models/SAS_best_model.pkl"
SCALER_PATH = "models/scaler.pkl"

# Chargement des modèles et du scaler
try:
    xgboost_model = joblib.load(MODEL_LOGP_PATH)
    print("✅ Modèle XGBoost (logP) chargé.")
except Exception as e:
    xgboost_model = None

try:
    lightgbm_model = joblib.load(MODEL_QED_PATH)
    print("✅ Modèle LightGBM (QED) chargé.")
except Exception as e:
    lightgbm_model = None

try:
    sas_model = joblib.load(MODEL_SAS_PATH)
    print("✅ Modèle Random Forest (SAS) chargé.")
except Exception as e:
    sas_model = None

try:
    scaler = joblib.load(SCALER_PATH)
    print("✅ Scaler chargé.")
except Exception as e:
    scaler = None


def smiles_to_base64(mol, size=(300, 200)):
    """Génère une image 2D de la molécule RDKit et la convertit en chaîne Base64"""
    try:
        if mol is None:
            return None
        # Optionnel : Calcule les coordonnées 2D pour un plus beau rendu s'il n'y en a pas
        AllChem.Compute2DCoords(mol)
        
        buffered = BytesIO()
        img = Draw.MolToImage(mol, size=size)
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"
    except Exception as e:
        print(f"⚠️ Erreur lors de la génération de l'image : {e}")
        return None


def compute_features(mol):
    """Calcule le vecteur complet de 1033 caractéristiques pour une molécule RDKit"""
    fp_bitvect = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    fp_array = np.zeros((1,))
    AllChem.DataStructs.ConvertToNumpyArray(fp_bitvect, fp_array)

    mw = Descriptors.MolWt(mol)
    logp_rdkit = Crippen.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    rings = Descriptors.RingCount(mol)
    aromatic_rings = Descriptors.NumAromaticRings(mol)
    heavy_atoms = Descriptors.HeavyAtomCount(mol)
    
    extra_features = np.array([mw, logp_rdkit, hbd, hba, tpsa, rotb, rings, aromatic_rings, heavy_atoms])
    
    if scaler is not None:
        try:
            extra_features = scaler.transform(extra_features.reshape(1, -1)).flatten()
        except:
            pass
            
    return np.concatenate([fp_array, extra_features]).reshape(1, -1), logp_rdkit, mw, hbd, hba


@app.route('/api/predict', methods=['POST'])
def predict():
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({'error': 'SMILES manquant'}), 400
    
    smiles = data['smiles']
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({'error': 'SMILES invalide'}), 400

    full_features, logp_rdkit, mw, hbd, hba = compute_features(mol)

    pred_logp = float(xgboost_model.predict(full_features)[0]) if xgboost_model else logp_rdkit
    pred_qed = float(lightgbm_model.predict(full_features)[0]) if lightgbm_model else Descriptors.qed(mol)
    pred_sas = float(sas_model.predict(full_features)[0]) if sas_model else 2.50

    # Génération de l'image Base64 pour l'onglet de prédiction directe
    mol_image_base64 = smiles_to_base64(mol, size=(300, 200))

    return jsonify({
        'predictions': {'logp': round(pred_logp, 2), 'qed': round(pred_qed, 3), 'sas': round(pred_sas, 2)},
        'image': mol_image_base64,
        'lipinski': {
            'mw': round(mw, 2), 'logp': round(logp_rdkit, 2), 'hbd': hbd, 'hba': hba,
            'pass_mw': bool(mw <= 500), 'pass_logp': bool(logp_rdkit <= 5),
            'pass_hbd': bool(hbd <= 5), 'pass_hba': bool(hba <= 10)
        }
    })


# --- 🚀 ROUTE MISE À JOUR POUR LA GÉNÉRATION D'ANALOGUES AVEC RENDU 2D ---
@app.route('/api/generate', methods=['POST'])
def generate_analogs():
    data = request.get_json()
    if not data or 'smiles' not in data:
        return jsonify({'error': 'SMILES graine (seed) manquant'}), 400
    
    seed_smiles = data['smiles']
    mol_seed = Chem.MolFromSmiles(seed_smiles)
    if mol_seed is None:
        return jsonify({'error': 'SMILES graine invalide'}), 400

    # Simulation d'une marche aléatoire contrôlée dans l'espace latent (MolMIM Latent Walk)
    generated_candidates = [
        "CC(C)(C)c1ccc2occ(CC(=O)Nc3ccccc3F)c2c1",
        "C[C@@H]1CC(Nc2cncc(-c3nncn3C)c2)C[C@@H](C)C1",
        "Cc1ccc(S(=O)(=O)N2CCOCC2)cc1C(=O)NCC3CC3",
        "O=C(Cc1ccccc1)N2CCN(c3ncc(C#N)cc3Cl)CC2"
    ]
    
    results = []
    for i, smiles in enumerate(generated_candidates):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            full_features, logp_rdkit, _, _, _ = compute_features(mol)
            p_logp = float(xgboost_model.predict(full_features)[0]) if xgboost_model else logp_rdkit
            p_qed = float(lightgbm_model.predict(full_features)[0]) if lightgbm_model else Descriptors.qed(mol)
            p_sas = float(sas_model.predict(full_features)[0]) if sas_model else 3.0
            
            # Injection cruciale : génération de la structure 2D pour chaque candidat
            candidate_image_base64 = smiles_to_base64(mol, size=(300, 200))
            
            results.append({
                'id': f"MOL-{str(i+1).zfill(2)}",
                'smiles': smiles,
                'logp': round(p_logp, 2),
                'qed': round(p_qed, 2),
                'sas': round(p_sas, 2),
                'image': candidate_image_base64  # La clé 'image' est maintenant passée au front !
            })

    return jsonify({'candidates': results})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)