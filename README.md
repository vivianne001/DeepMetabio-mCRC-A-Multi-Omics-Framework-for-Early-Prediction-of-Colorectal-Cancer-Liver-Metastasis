
# DeepMetabio-mCRC Screener: A Multi-Omics Framework for Early Prediction of Colorectal Liver Metastasis

This repository contains the code and analysis pipeline for the "DeepMetabio-mCRC Screener" project. This study develops and validates an integrated multi-omics framework to predict colorectal liver metastasis (CRLM) and identifies novel, multifunctional biomarker for CRLM diagnosis, prognosis, and therapeutic strategy.

---

## Abstract

Colorectal liver metastasis (CRLM) is the leading cause of mortality in colorectal cancer, yet early predictive tools or biomarkers remain scarce. We present DeepMetabio‑mCRC Screener, an integrated multi‑omics framework coupling large‑scale transcriptomic profiles with serum metabolomics. Using 1,077 CRC samples, 620 metabolism‑related genes trained a convolutional neural network, achieving AUROC = 0.92 in validation and 0.97 in independent testing, outperforming 10 classical machine‑learning models. Model‑derived transcriptomic risk scores (TRS) identified 22 core metabolic features contributing to mCRC and CRLM occurrence, enriched in retinol and tryptophan metabolism. Cross‑omics analysis identified ACMSD as a novel biomarker, mechanistically linked to impaired NAD⁺ biosynthesis. Clinical validation in 100 CRC patients confirmed elevated ACMSD in CRLM, correlating with advanced stage, recurrence risk, immune‑inflamed tumor microenvironment, and increased sensitivity to EGFR/VEGFR‑targeted therapy. Our study establishes DeepMetabio‑mCRC Screener as a robust early‑detection tool and positions ACMSD as a multifunctional biomarker for CRLM diagnosis, prognosis, molecular analysis and therapeutic strategy.

---

## Framework

The project follows a systematic, multi-stage analysis pipeline, from data acquisition to multi-omics integration and clinical application.

```
+--------------------------------+
|   Data Acquisition & Curation    |
| (Multiple GEO & Nat.Med Datasets)|
+--------------------------------+
                 |
                 v
+--------------------------------+
|   1. Functional Gene Pooling   |
|  (DEG Analysis & Metabolic Gene|
|          Intersection)         |
+--------------------------------+
                 |
                 v
+--------------------------------+
|   2. Model Development          |
| (Data Integration, ComBat,      |
|  1D-CNN Training with 5-fold CV)|
+--------------------------------+
                 |
                 v
+--------------------------------+
|   3. Model Validation          |
| (Internal & External Cohorts,  |
|   Benchmark vs. 10 ML Models)  |
+--------------------------------+
                 |
                 v
+--------------------------------+
|   4. Biomarker Discovery       |
| (TRS Calculation, Correlation, |
|   Core Feature Selection)      |
+--------------------------------+
                 |
                 v
+--------------------------------+
|   5. Multi-Omics Integration   |
| (Transcriptome-Metabolome      |
|      Pathway Intersection)     |
+--------------------------------+
                 |
                 v
+--------------------------------+
|   6. Clinical & Mechanistic    |
|          Validation            |
| (Survival, Stage, TME, GSEA)   |
+--------------------------------+
                 |
                 v
+--------------------------------+
|   7. Precision Medicine App    |
| (OncoPredict on Patients,      |
|      GDSC on Cell Lines)       |
+--------------------------------+
                 |
                 v
+--------------------------------+
|  Final Biomarker: ACMSD        |
+--------------------------------+
```

---

## Workflow and Code Structure

This repository is organized into directories corresponding to the major stages of the analysis pipeline.

### 1. Functional Gene Pooling
*   **Description:** This stage identifies differentially expressed genes (DEGs) between metastatic and primary CRC samples from the GSE131418 dataset and intersects them with a known list of metabolic genes to create a "functional gene pool" for model training.
*   **Scripts:**
    *   `Functional Gene Pooling/script-data preprocessing.R`: Downloads and preprocesses the GSE131418 dataset.
    *   `Functional Gene Pooling/Differentially expression analysis and Pathway Enrichment.R`: Performs differential expression analysis (mCRC vs. Primary CRC) and GO enrichment.
    *   `Functional Gene Pooling/Identification and Visualization of Functional DEGs.R`: Intersects DEGs with metabolic genes to define the functional gene pool.

### 2. Model Development and Validation
*   **Description:** This core stage involves integrating multiple transcriptomic datasets, correcting for batch effects, and then building, training, and validating the 1D-Convolutional Neural Network (CNN) model.
*   **Scripts:**
    *   `Samples merge and split/data integration and combat-new.R`: Merges multiple GEO datasets and applies the ComBat algorithm to remove batch effects. PCA and silhouette scores are used to evaluate the correction.
    *   `Model/Model Construct-Training-Cross Validation.ipynb`: Uses the batch-corrected data and the functional gene pool to train a 1D-CNN classifier with a 5-fold stratified cross-validation strategy.
    *   `Model/External Validation.ipynb`: Applies the trained ensemble of 5 CNN models to a completely independent external validation cohort (GSE50760) to assess its generalization performance.
    *   `Benchmark Testing/script-10ml_benchmark_test.ipynb`: Compares the performance of the 1D-CNN model against 10 classical machine learning models (e.g., SVM, Random Forest, XGBoost) on the same datasets to demonstrate its superiority.

### 3. Biomarker Discovery and Validation
*   **Description:** This stage uses the trained model to calculate a Transcriptomic Risk Score (TRS) for each sample and identifies the core genes most strongly associated with this risk score.
*   **Scripts:**
    *   `Model/script-extract_trs.ipynb` & `Model/script-correlate_trs_with_genes.ipynb`: Calculate the TRS for each sample and compute the Pearson correlation between the TRS and the expression of each input gene.
    *   `Biomarker Selection/Core TRS Features/script-selection_core_trs_features.R`: Identifies the top 50 genes most correlated with TRS in both the training and external validation sets, and finds the intersection to define robust "Core TRS Features".
    *   `Biomarker Selection/script-ACMSD_batch_single_gene_AUROC.R`: Evaluates the diagnostic performance (AUROC) of individual core features, comparing them with known CRLM biomarkers to highlight the superior performance of `ACMSD`.

### 4. Multi-Omics Integration
*   **Description:** This stage integrates the transcriptomic pathway analysis with results from an independent serum metabolomics study to identify core dysregulated metabolic networks in CRLM.
*   **Scripts:**
    *   `Metabolic Network/script-core_trs_features_enrichment.R`: Performs KEGG pathway enrichment on the "Core TRS Features".
    *   `Metabolic Network/script-transcriptome_metabolome_pathway_intersection.R` & `Metabolic Network/script-final_integration_analysis_ORA_version.ipynb`: Intersects the KEGG pathways enriched in the transcriptome with those enriched from serum metabolomics data, visualizing the overlap to identify the most critical metabolic pathways.

### 5. Clinical and Immune Microenvironment (TME) Analysis
*   **Description:** This stage explores the clinical relevance and potential biological mechanisms of the key biomarker, `ACMSD`.
*   **Scripts:**
    *   `Clinical and Immune TME Analysis/Clinical Analysis/Kaplan-Meier Survival Analysis.R`: Validates the prognostic value of `ACMSD` by performing Kaplan-Meier survival analysis in an independent CRC cohort (GSE17536).
    *   `Clinical and Immune TME Analysis/Clinical Analysis/Gene Expression Analysis Across Pathological Stages in TCGA.R`: Analyzes the correlation between `ACMSD` expression and tumor progression using pathological stage data from the TCGA-CRC dataset.
    *   `Clinical and Immune TME Analysis/Immune TME Analysis/`: A series of scripts to dissect the tumor immune microenvironment.
        *   `ssGSEA Analysis...R`: Uses ssGSEA to quantify immune cell infiltration and correlates it with `ACMSD` expression.
        *   `multibox plot analysis...R`: Compares the expression of key immune and inflammatory factors between `ACMSD`-high and `ACMSD`-low patient groups.
        *   `GSEA analysis of pathway.R`: Performs GSEA to identify immune-related pathways activated in `ACMSD`-high patients.
        *   `Calculation of ESTIMATE Scores...R`: Uses the ESTIMATE algorithm to assess overall immune and stromal infiltration and its relationship with `ACMSD`.



### 6. Precision Medicine Application
*   **Description:** This final stage explores the translational potential of `ACMSD` as a predictive biomarker for therapy selection.
*   **Scripts:**
    *   `Precision Medicine/Oncopredict/`: Uses the `OncoPredict` R package to predict patient-specific drug sensitivity based on tumor gene expression.
        *   `Drug Sensitivity Prediction...R`: Predicts IC50 values for various drugs in CRLM patients.
        *   `Visualization of Drug Sensitivity...R`: Compares the predicted sensitivity to targeted therapies (e.g., VEGFR/EGFR inhibitors) between `ACMSD`-high and `ACMSD`-low patient groups.
    *   `Precision Medicine/gdsc/`: Validates drug sensitivity findings using the GDSC cancer cell line database.
    *   `Precision Medicine/gdsc/script-correlation_group_comparison_analysis.R`: Correlates `ACMSD` expression with drug IC50 values across dozens of CRC cell lines.
    *   `Precision Medicine/gdsc/script-drug_target_attribute_classification.R` & `Precision Medicine/gdsc/script-GDSC_CTRP_target_network_analysis.R`: Analyzes and visualizes the targets of `ACMSD`-correlated drugs to uncover mechanistic links to specific drug classes.

---

## Requirements

*   **R** (version 4.0 or higher) and the following packages:
    *   `GEOquery`, `limma`, `BiocManager`
    *   `ggplot2`, `ggvenn`, `ggpubr`, `pheatmap`, `corrplot`
    *   `sva`, `GSVA`, `IOBR`, `estimate`
    *   `survival`, `survminer`
    *   `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`
    *   `tidyverse`, `dplyr`, `reshape2`, `readxl`
    *   `oncoPredict`

*   **Python** (version 3.8 or higher) and the following packages:
    *   `pandas`, `numpy`, `scikit-learn`
    *   `tensorflow` (version 2.x)
    *   `matplotlib`, `seaborn`
    *   `joblib`, `openpyxl`
    *   `lightgbm`, `xgboost`

## How to Use

1.  **Clone the repository:**
    ```bash
    git clone [repository-url]
    ```
2.  **Install dependencies:** Ensure all required R and Python packages are installed.
3.  **Prepare data:** Download the required raw data from GEO/TCGA or place the provided intermediate data files in the correct directories.
4.  **Run the scripts:** Follow the numerical order of the directories. The scripts within each directory are designed to be run sequentially to reproduce the analysis. Please adjust file paths within the scripts to match your local environment.

---

## Citation

If you use the code or findings from this project in your research, please cite our publication (details to be added upon acceptance).
