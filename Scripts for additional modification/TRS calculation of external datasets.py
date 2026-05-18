"""
TRS prediction using cross-validated CNN ensemble

Purpose:
1. Load the ensemble of 5 CV-trained models
2. Load a new expression matrix (rows=samples, cols=genes)
3. Predict TRS scores only
4. Save TRS results
"""

import pandas as pd
import numpy as np
import tensorflow as tf
import random
import os
import joblib
import json
import warnings

warnings.filterwarnings('ignore')

print("--- TRS prediction using cross-validated CNN ensemble ---")

# =========================================================
# Step 0: Global settings
# =========================================================
SEED = 42
os.environ['PYTHONHASHSEED'] = str(SEED)
random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)

print("--- Random seed set to 42 ---")

# =========================================================
# Step 1: File paths
# =========================================================
print("\n--- Step 1: Configure file paths ---")

# Note: Local paths have been translated to English for SCI submission.
BASE_DIR = r"D:/CRC_Liver_Metastasis_Biomarker_Diagnosis/New_Strategy/Autoencoder"

# Cross-validation model output paths
CV_OUTPUT_DIR = r"D:/temp_output_cv"
ENSEMBLE_MODEL_DIR = os.path.join(CV_OUTPUT_DIR, "ensemble_models")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Modify here: Your expression matrix file to be predicted
# Requirements: Rows = Samples, Columns = Genes
# The first column / index column must be the sample ID
NEW_DATA_FILE = r"D:/CRC_Liver_Metastasis_Biomarker_Diagnosis/New_Strategy/Autoencoder/validation_datasets/TCGA_exp_matrix.csv"
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Results output path
OUTPUT_DIR = os.path.join(BASE_DIR, "trs_prediction_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Required files
CV_LABEL_ENCODER_FILE = os.path.join(ENSEMBLE_MODEL_DIR, "label_encoder.pkl")
CV_GENES_FILE = os.path.join(CV_OUTPUT_DIR, "used_functional_genes_cv.txt")

print(f"Ensemble model directory: {ENSEMBLE_MODEL_DIR}")
print(f"Input expression matrix: {NEW_DATA_FILE}")
print(f"Results output directory: {OUTPUT_DIR}")

required_files = [CV_LABEL_ENCODER_FILE, CV_GENES_FILE, NEW_DATA_FILE]

if not os.path.exists(ENSEMBLE_MODEL_DIR):
    print(f"ERROR: Ensemble model directory not found: {ENSEMBLE_MODEL_DIR}")
    print("Please make sure the trained models exist.")
    raise SystemExit(1)

missing_files = [f for f in required_files if not os.path.exists(f)]
if missing_files:
    print("ERROR: The following required files are missing:")
    for f in missing_files:
        print(f"  - {f}")
    raise SystemExit(1)
else:
    print("All required files are present.")

# =========================================================
# Step 2: Load ensemble models and preprocessors
# =========================================================
print("\n--- Step 2: Load ensemble models and associated preprocessors ---")

try:
    models = []
    scalers = []

    for i in range(1, 6):
        model_path = os.path.join(ENSEMBLE_MODEL_DIR, f"model_fold{i}.keras")
        scaler_path = os.path.join(ENSEMBLE_MODEL_DIR, f"scaler_fold{i}.pkl")

        if not os.path.exists(model_path) or not os.path.exists(scaler_path):
            raise FileNotFoundError(f"Model or scaler for fold {i} not found.")

        models.append(tf.keras.models.load_model(model_path))
        scalers.append(joblib.load(scaler_path))
        print(f"  - Loaded model and scaler for fold {i}")

    print(f"Loaded {len(models)} cross-validation models successfully.")

    # Keep consistent with original script; although not used for scoring only
    cv_label_encoder = joblib.load(CV_LABEL_ENCODER_FILE)

    # Load gene list used for modeling
    with open(CV_GENES_FILE, 'r') as f:
        cv_model_genes = [line.strip() for line in f.readlines() if line.strip()]

    print(f"Functional gene list loaded: {len(cv_model_genes)} genes")
    print(f"Model input shape: {models[0].input_shape}")

except Exception as e:
    print(f"ERROR: Failed to load models or resources: {e}")
    raise

# =========================================================
# Step 3: Load new dataset
# =========================================================
print("\n--- Step 3: Load new expression matrix ---")

# Assume the first column is the sample ID by default
new_data = pd.read_csv(NEW_DATA_FILE, index_col=0)

print(f"Input data shape: {new_data.shape}")
print(f"Number of samples: {new_data.shape[0]}")
print(f"Number of genes in input: {new_data.shape[1]}")

# Ensure data type is numeric
new_data = new_data.apply(pd.to_numeric, errors='coerce')

# Warn if there are missing values
if new_data.isnull().any().any():
    print("WARNING: Missing values detected in input matrix. They will be filled with 0.")
    new_data = new_data.fillna(0)

# =========================================================
# Step 4: Build feature matrix in the exact training-gene order
# =========================================================
print("\n--- Step 4: Prepare feature matrix in training gene order ---")

available_genes = [gene for gene in cv_model_genes if gene in new_data.columns]
missing_genes = [gene for gene in cv_model_genes if gene not in new_data.columns]

print(f"Of {len(cv_model_genes)} required model genes:")
print(f"  - Found in input data: {len(available_genes)} ({len(available_genes)/len(cv_model_genes)*100:.1f}%)")
print(f"  - Missing: {len(missing_genes)} ({len(missing_genes)/len(cv_model_genes)*100:.1f}%)")

if len(missing_genes) > 0:
    print("WARNING: Missing genes will be filled with zeros.")
    print(f"First 20 missing genes: {missing_genes[:20]}")

# Build feature matrix in same column order as training
X_new = pd.DataFrame(index=new_data.index, columns=cv_model_genes)

# Fill available genes
X_new[available_genes] = new_data[available_genes]

# Fill missing genes with 0
if missing_genes:
    X_new[missing_genes] = 0

# Ensure numeric and no NA
X_new = X_new.apply(pd.to_numeric, errors='coerce').fillna(0)

print(f"Final feature matrix shape: {X_new.shape}")
print(f"Contains any NA: {X_new.isnull().any().any()}")

# =========================================================
# Step 5: Ensemble TRS prediction
# =========================================================
print("\n--- Step 5: Predict TRS scores with 5-model ensemble ---")

all_predictions = []

for i, (model, scaler) in enumerate(zip(models, scalers), start=1):
    # Apply fold-specific scaler (same as original script)
    X_scaled = scaler.transform(X_new)

    # CNN expects shape: (samples, features, 1)
    X_cnn = np.expand_dims(X_scaled, axis=-1)

    pred_proba = model.predict(X_cnn, verbose=0).flatten()
    all_predictions.append(pred_proba)

    print(f"  - Prediction completed for fold {i}")

# Average predicted probabilities across 5 models
trs_scores_raw = np.mean(all_predictions, axis=0)
print(f"\nRaw ensemble score range: [{trs_scores_raw.min():.4f}, {trs_scores_raw.max():.4f}]")

# Keep exactly the same semantic correction as original script
trs_scores_final = 1 - trs_scores_raw
print(f"Final TRS score range: [{trs_scores_final.min():.4f}, {trs_scores_final.max():.4f}]")

print("TRS prediction completed.")

# =========================================================
# Step 6: Save TRS results
# =========================================================
print("\n--- Step 6: Save TRS results ---")

trs_results = pd.DataFrame({
    "Sample_ID": X_new.index,
    "TRS_Score": trs_scores_final
})

# Optional: Keep the scores for each fold to facilitate subsequent checks
for i, fold_pred in enumerate(all_predictions, start=1):
    trs_results[f"Fold{i}_RawScore"] = fold_pred
    trs_results[f"Fold{i}_TRS"] = 1 - fold_pred

trs_results_file = os.path.join(OUTPUT_DIR, "predicted_TRS_scores.csv")
trs_results.to_csv(trs_results_file, index=False)

summary_info = {
    "Model_Type": "Cross_Validation_Ensemble_CNN",
    "Total_Samples": int(X_new.shape[0]),
    "Required_Features": int(len(cv_model_genes)),
    "Available_Features": int(len(available_genes)),
    "Missing_Features": int(len(missing_genes)),
    "Missing_Feature_Ratio": float(len(missing_genes) / len(cv_model_genes)),
    "TRS_Range": {
        "Min": float(trs_scores_final.min()),
        "Max": float(trs_scores_final.max()),
        "Mean": float(np.mean(trs_scores_final)),
        "Median": float(np.median(trs_scores_final))
    }
}

summary_file = os.path.join(OUTPUT_DIR, "trs_prediction_summary.json")
with open(summary_file, "w", encoding="utf-8") as f:
    json.dump(summary_info, f, indent=4, ensure_ascii=False)

print(f"TRS result file saved to: {trs_results_file}")
print(f"Summary file saved to: {summary_file}")

print("\nPreview of predicted TRS:")
print(trs_results[["Sample_ID", "TRS_Score"]].head(10))

print("\n" + "="*80)
print("🎉 TRS prediction finished successfully!")
print("="*80)
print(f"All outputs are saved to: {OUTPUT_DIR}")