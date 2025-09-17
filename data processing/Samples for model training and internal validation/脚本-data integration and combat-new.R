rm(list=ls())
# data_integration_and_combat.R
# Purpose: integrate GEO and TCGA data, perform batch correction, and save cleaned data for downstream Python modeling
#
# Before running ensure:
# 1. R objects `exp`, `ph`, and `dat_silu` are available in the environment (GEO expression, GEO phenotype, TCGA expression)
# 2. `sva` package is installed (if not: BiocManager::install("sva"))
# 3. `limma` package is installed (if not: BiocManager::install("limma"))

# --- Check and install required packages ---
if (!require("sva", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("sva")
  library(sva)
}

if (!require("limma", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("limma")
  library(limma)
}

# Load visualization packages
required_packages <- c("ggplot2", "RColorBrewer", "pheatmap", "corrplot")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# --- Step 1: Data checks and preprocessing ---
message("--- Step 1: Validate and preprocess data ---")

# Check required objects
if (!exists("exp")) {
  stop("Error: 'exp' object not found. Please load GEO expression matrix.")
}

dat_silu = read.csv("dat_silu.csv")
dat_silu = as.data.frame(dat_silu)
rownames(dat_silu) = dat_silu$X
dat_silu = dat_silu[,-1]
if (!exists("ph")) {
  stop("Error: 'ph' object not found. Please load GEO phenotype/metadata.")
}
if (!exists("dat_silu")) {
  stop("Error: 'dat_silu' object not found. Please ensure TCGA data is available.")
}

message(paste("GEO expression matrix dimensions:", paste(dim(exp), collapse=" x ")))
# note: exp_geo variable not always present; keep message informative
if (exists("exp_geo")) message(paste("GEO transposed expression dimensions (exp_geo):", paste(dim(exp_geo), collapse=" x ")))
message(paste("TCGA expression matrix dimensions (dat_silu):", paste(dim(dat_silu), collapse=" x ")))

# --- Step 2: Harmonize data formats ---
message("--- Step 2: Harmonize data formats ---")

# For GEO: ensure exp column names match sample IDs in ph
if (!all(colnames(exp) %in% ph$geo_accession)) {
  message("Warning: GEO expression sample names do not fully match phenotype data. Attempting to subset to common samples...")
  common_samples_geo <- intersect(colnames(exp), ph$geo_accession)
  exp <- exp[, common_samples_geo]
  ph <- ph[ph$geo_accession %in% common_samples_geo, ]
}

# Transpose GEO so rows = samples, columns = genes
exp_geo_t <- t(exp)
exp_geo_df <- as.data.frame(exp_geo_t)

# For TCGA: separate status column
status_tcga <- dat_silu$status
exp_tcga_df <- dat_silu[, !colnames(dat_silu) %in% "status", drop = FALSE]

# For other GEO-like object (exp_geo) separate status if present
status_geo <- NULL
if (exists("exp_geo")) {
  if ("status" %in% colnames(exp_geo)) {
    status_geo <- exp_geo$status
    exp_gse_df <- exp_geo[, !colnames(exp_geo) %in% "status", drop = FALSE]
  } else {
    exp_gse_df <- exp_geo
  }
} else {
  # If exp_geo not present, create empty data frame
  exp_gse_df <- data.frame()
}

# --- Step 3: Find common genes and align datasets ---
message("--- Step 3: Find common genes and align datasets ---")

# Debugging for ACMSD presence
message("--- Debugging gene presence: ACMSD ---")
geo_genes <- colnames(exp_geo_df)
tcga_genes <- colnames(exp_tcga_df)
gse_genes <- if (ncol(exp_gse_df) > 0) colnames(exp_gse_df) else character(0)

is_acmsd_in_geo <- "ACMSD" %in% geo_genes
is_acmsd_in_tcga <- "ACMSD" %in% tcga_genes
is_acmsd_in_gse <- "ACMSD" %in% gse_genes

message(paste("Is 'ACMSD' in GEO data?", is_acmsd_in_geo))
message(paste("Is 'ACMSD' in TCGA data?", is_acmsd_in_tcga))

if (!is_acmsd_in_geo) {
  message("'ACMSD' not found in GEO. Searching for case-insensitive matches...")
  print(grep("acmsd", geo_genes, value = TRUE, ignore.case = TRUE))
}

if (!is_acmsd_in_tcga) {
  message("'ACMSD' not found in TCGA. Searching for case-insensitive matches...")
  print(grep("acmsd", tcga_genes, value = TRUE, ignore.case = TRUE))
}
message("--- End Debugging ---")

# Find intersecting genes across datasets
common_genes <- intersect(colnames(exp_geo_df), colnames(exp_tcga_df))
if (length(exp_gse_df) > 0) common_genes <- intersect(common_genes, colnames(exp_gse_df))
message(paste("Found", length(common_genes), "common genes"))

if (length(common_genes) < 100) {
  stop("Error: Too few common genes found. Please check gene name conventions.")
}

# Keep only common genes and align columns
exp_geo_aligned <- exp_geo_df[, common_genes, drop = FALSE]
exp_tcga_aligned <- exp_tcga_df[, common_genes, drop = FALSE]
exp_gse_aligned <- if (ncol(exp_gse_df) > 0) exp_gse_df[, common_genes, drop = FALSE] else data.frame()

# --- Step 4: Build combined expression matrix and metadata ---
message("--- Step 4: Combine datasets ---")

# Combine expression matrices (rows = samples)
exp_combined <- rbind(exp_geo_aligned, exp_tcga_aligned)
if (nrow(exp_gse_aligned) > 0) exp_combined <- rbind(exp_combined, exp_gse_aligned)

# Build metadata for each source
meta_geo <- data.frame(
  sample_id = rownames(exp_geo_aligned),
  group = ph$group[match(rownames(exp_geo_aligned), ph$geo_accession)],
  batch = "GeneChip Platform GEO Samples",
  stringsAsFactors = FALSE
)

meta_silu <- data.frame(
  sample_id = rownames(exp_tcga_aligned),
  group = status_tcga[match(rownames(exp_tcga_aligned), names(status_tcga))],
  batch = "RNA-Seq Platform COAD_Silu Samples",
  stringsAsFactors = FALSE
)

meta_gse <- if (nrow(exp_gse_aligned) > 0) data.frame(
  sample_id = rownames(exp_gse_aligned),
  group = status_geo[match(rownames(exp_gse_aligned), names(status_geo))],
  batch = "RNA-Seq Platform GEO Samples",
  stringsAsFactors = FALSE
) else data.frame()

identical(meta_silu$sample_id, rownames(dat_silu))
if (nrow(meta_silu) > 0) meta_silu$group = dat_silu$status
if (nrow(meta_gse) > 0) meta_gse$group = if (!is.null(status_geo)) status_geo else meta_gse$group

meta_combined <- rbind(meta_geo, meta_silu)
if (nrow(meta_gse) > 0) meta_combined <- rbind(meta_combined, meta_gse)
rownames(meta_combined) <- meta_combined$sample_id

# Ensure meta and expression have matching order
meta_combined <- meta_combined[rownames(exp_combined), , drop = FALSE]

message("Combined data dimensions:")
message(paste("Expression matrix:", paste(dim(exp_combined), collapse=" x ")))
message(paste("Metadata:", paste(dim(meta_combined), collapse=" x ")))
message("Batch distribution:")
print(table(meta_combined$batch))
message("Group distribution:")
print(table(meta_combined$group))

# --- Step 5: Run ComBat batch correction ---
message("--- Step 5: Execute ComBat batch correction ---")

# ComBat requires genes as rows and samples as columns
exp_for_combat <- t(exp_combined)

# Build model to protect biological variable 'group'
mod <- model.matrix(~ group, data = meta_combined)
batch_vector <- meta_combined$batch

message("Starting ComBat correction...")
exp_combat_corrected <- ComBat(
  dat = exp_for_combat,
  batch = batch_vector,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

# Transpose back to samples x genes
exp_final <- as.data.frame(t(exp_combat_corrected))

message("Batch correction completed!")

# --- Step 5.5: Visual comparison before/after correction ---
message("--- Step 5.5: Generate comparison plots before/after batch correction ---")

# Create output directory for plots (plots_dir and output_dir should be defined before use)
if (!exists("output_dir")) output_dir <- getwd()
plots_dir <- file.path(output_dir, "batch_correction_plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# === 1. PCA comparison ===
message("Generating PCA comparison plots...")

# Select top variable genes for PCA (e.g., top 2000)
gene_vars_before <- apply(exp_for_combat, 1, var)
top_var_genes <- names(sort(gene_vars_before, decreasing = TRUE))[1:min(2000, length(gene_vars_before))]

exp_before_pca <- t(exp_for_combat[top_var_genes, , drop = FALSE])
exp_after_pca <- t(exp_combat_corrected[top_var_genes, , drop = FALSE])

# Perform PCA
pca_before <- prcomp(exp_before_pca, center = TRUE, scale. = TRUE)
pca_after <- prcomp(exp_after_pca, center = TRUE, scale. = TRUE)

# Prepare PCA plotting data
pca_before_df <- data.frame(
  PC1 = pca_before$x[,1],
  PC2 = pca_before$x[,2],
  Batch = meta_combined$batch,
  Group = meta_combined$group,
  Sample = rownames(meta_combined)
)

pca_after_df <- data.frame(
  PC1 = pca_after$x[,1],
  PC2 = pca_after$x[,2],
  Batch = meta_combined$batch,
  Group = meta_combined$group,
  Sample = rownames(meta_combined)
)

# Variance explained
var_explained_before <- round(100 * pca_before$sdev^2 / sum(pca_before$sdev^2), 2)
var_explained_after <- round(100 * pca_after$sdev^2 / sum(pca_after$sdev^2), 2)

# PCA plots colored by batch
p1 <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("GeneChip Platform GEO Samples" = "#E31A1C", "RNA-Seq Platform COAD_Silu Samples" = "#1F78B4","RNA-Seq Platform GEO Samples" = "yellowgreen")) +
  scale_shape_manual(values = c("primary" = 16, "metastasis" = 17)) +
  labs(title = "PCA - Before Batch Correction",
       x = paste0("PC1 (", var_explained_before[1], "% variance)"),
       y = paste0("PC2 (", var_explained_before[2], "% variance)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

p2 <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("GeneChip Platform GEO Samples" = "#E31A1C", "RNA-Seq Platform COAD_Silu Samples" = "#1F78B4","RNA-Seq Platform GEO Samples" = "yellowgreen")) +
  scale_shape_manual(values = c("primary" = 16, "metastasis" = 17)) +
  labs(title = "PCA - After Batch Correction",
       x = paste0("PC1 (", var_explained_after[1], "% variance)"),
       y = paste0("PC2 (", var_explained_after[2], "% variance)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

# PCA plots colored by group
p3 <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("primary" = "#2E8B57", "metastasis" = "#DC143C")) +
  scale_shape_manual(values = c("GEO" = 16, "TCGA" = 17)) +
  labs(title = "PCA - Before Correction (colored by group)",
       x = paste0("PC1 (", var_explained_before[1], "% variance)"),
       y = paste0("PC2 (", var_explained_before[2], "% variance)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

p4 <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("primary" = "#2E8B57", "metastasis" = "#DC143C")) +
  scale_shape_manual(values = c("GEO" = 16, "TCGA" = 17)) +
  labs(title = "PCA - After Correction (colored by group)",
       x = paste0("PC1 (", var_explained_after[1], "% variance)"),
       y = paste0("PC2 (", var_explained_after[2], "% variance)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

# Save PCA figures (requires gridExtra::arrangeGrob)
if (!require("gridExtra", quietly = TRUE)) install.packages("gridExtra"); library(gridExtra)
ggsave(file.path(plots_dir, "PCA_batch_before_after.png"), 
       arrangeGrob(p1, p2, ncol = 2), width = 16, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "PCA_group_before_after.png"), 
       arrangeGrob(p3, p4, ncol = 2), width = 16, height = 6, dpi = 300)



# === 1.5 Batch effect quantitative metrics ===
message("Computing batch effect quantitative metrics...")

# 1. silhouette score
if (!require("cluster", quietly = TRUE)) install.packages("cluster")
library(cluster)
# Use PCA components for silhouette (choose appropriate dimensions)
sil_before <- silhouette(as.numeric(as.factor(meta_combined$batch)), dist(pca_before$x[,1:2]))
sil_after <- silhouette(as.numeric(as.factor(meta_combined$batch)), dist(pca_after$x[,1:2]))
sil_score_before <- mean(sil_before[,3])
sil_score_after <- mean(sil_after[,3])

sil_before = as.data.frame(sil_before)
sil_after = as.data.frame(sil_after)

save(sil_before, sil_after, file = "silhouette_Score_Matrix.rdata")

