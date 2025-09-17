
###########################################################
# Data Integration and Batch Correction (ComBat) - Publication-ready
# Purpose:
# - Integrate GEO and silu expression data
# - Align genes and samples, perform ComBat batch correction
# - Produce PCA / t-SNE / UMAP visualizations before/after correction
# - Save corrected expression data and plots for downstream analysis
#
# Requirements:
# - R objects or CSV/RData files providing expression matrices and phenotype:
#   * GEO expression matrix 'exp' (genes x samples) and phenotype 'ph' with 'geo_accession'
#   * silu expression matrix 'dat_silu' (samples x genes) with 'status' column (or adjust accordingly)
# - Installed packages: sva, limma, ggplot2, RColorBrewer, pheatmap, corrplot, Rtsne, umap, cluster, gridExtra
###########################################################

rm(list = ls())

# ---------------------------
# Package checks / load
# ---------------------------
required_bioc <- c("sva", "limma")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_bioc) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

required_cran <- c("ggplot2", "RColorBrewer", "pheatmap", "corrplot", "Rtsne", "umap", "cluster", "gridExtra", "dplyr")
for (pkg in required_cran) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ---------------------------
# User settings
# ---------------------------
IO_DIR <- getwd()  # change to your data/working folder
plots_dir <- file.path(IO_DIR, "batch_correction_plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

OUT_COMBAT_RDS <- file.path(IO_DIR, "expression_data_combat_corrected.csv")
OUT_COMBINED_RDS <- file.path(IO_DIR, "metadata_combined.csv")
OUT_SIL_SCORE_RDATA <- file.path(IO_DIR, "silhouette Score Matrix.rdata")

# ---------------------------
# Data loading / checks
# ---------------------------
# Option A: If you have a processed RData that contains 'exp', 'ph', 'dat_silu', load it here:
# load(file.path(IO_DIR, "processed_expression_and_metadata.RData"))

# Option B: If the individual objects are already in the environment, script will use them.
if (!exists("exp")) stop("'exp' (GEO expression matrix, genes x samples) not found in environment.")
if (!exists("dat_silu") && !exists("exp_gse")) stop("'dat_silu' (silu expression) or 'exp_gse' not found.")
if (!exists("ph")) warning("'ph' (GEO phenotype) not found. Some sample filtering may be skipped.")

# If dat_silu not present but dat_silu.csv exists, user can load it:
if (!exists("dat_silu") && file.exists(file.path(IO_DIR, "dat_silu.csv"))) {
  dat_silu <- read.csv(file.path(IO_DIR, "dat_silu.csv"), stringsAsFactors = FALSE, row.names = 1)
  message("Loaded dat_silu from CSV.")
}

# If exp_gse or exp_geo object naming inconsistencies exist, standardize them
# Expectation: 'exp' is GEO expression (genes x samples)
# If you have another GEO expression matrix named 'exp_geo' or 'exp_gse', adapt below
if (exists("exp_gse") && !exists("exp_geo")) exp_geo <- exp_gse
if (!exists("exp_geo") && exists("exp")) {
  # If exp appears to be gene x samples as expected, we use it
  exp_geo <- exp
}

# Ensure dat_silu is samples x genes. If your dat_silu is genes x samples, transpose as needed.
# Here we assume dat_silu rows are samples (as in your original script)
if (is.null(rownames(dat_silu))) stop("dat_silu must have rownames as sample IDs.")

# ---------------------------
# Prepare GEO and silu expression frames: samples x genes
# ---------------------------
# For GEO: transpose 'exp_geo' if needed to samples x genes
if (nrow(exp_geo) > ncol(exp_geo)) {
  # if genes x samples, transpose to samples x genes
  exp_geo_df <- as.data.frame(t(exp_geo))
} else {
  exp_geo_df <- as.data.frame(exp_geo)
}

# For silu: if dat_silu has a 'status' column, separate it
if ("status" %in% colnames(dat_silu)) {
  status_silu <- dat_silu$status
  exp_silu_df <- dat_silu[, setdiff(colnames(dat_silu), "status"), drop = FALSE]
} else {
  exp_silu_df <- dat_silu
  status_silu <- rep(NA_character_, nrow(exp_silu_df))
  names(status_silu) <- rownames(exp_silu_df)
}

# If exp_gse exists separately (additional GEO samples), prepare it too
if (exists("exp_gse") && !is.null(exp_gse)) {
  if (nrow(exp_gse) > ncol(exp_gse)) {
    exp_gse_df <- as.data.frame(t(exp_gse))
  } else {
    exp_gse_df <- as.data.frame(exp_gse)
  }
} else {
  exp_gse_df <- data.frame()
}

# ---------------------------
# Debugging check for a key gene (ACMSD)
# ---------------------------
message("Checking presence of gene 'ACMSD' across datasets...")
geo_genes <- colnames(exp_geo_df)
silu_genes <- colnames(exp_silu_df)
gse_genes <- if (ncol(exp_gse_df) > 0) colnames(exp_gse_df) else character(0)

message("ACMSD in GEO:", "ACMSD" %in% geo_genes)
message("ACMSD in silu:", "ACMSD" %in% silu_genes)
if (length(gse_genes) > 0) message("ACMSD in additional GEO:", "ACMSD" %in% gse_genes)

# ---------------------------
# Find intersecting genes across datasets
# ---------------------------
common_genes <- Reduce(intersect, list(colnames(exp_geo_df), colnames(exp_silu_df), colnames(exp_gse_df)))
# If exp_gse_df is empty, ignore it
if (ncol(exp_gse_df) == 0) {
  common_genes <- intersect(colnames(exp_geo_df), colnames(exp_silu_df))
}

message("Found ", length(common_genes), " common genes across datasets.")
if (length(common_genes) < 100) stop("Too few common genes. Please check gene identifiers and formats (e.g., gene symbols).")

# Subset and align columns
exp_geo_aligned <- exp_geo_df[, common_genes, drop = FALSE]
exp_silu_aligned <- exp_silu_df[, common_genes, drop = FALSE]
if (ncol(exp_gse_df) > 0) exp_gse_aligned <- exp_gse_df[, common_genes, drop = FALSE] else exp_gse_aligned <- data.frame()

# ---------------------------
# Merge expression matrices (samples x genes)
# ---------------------------
exp_combined <- rbind(exp_geo_aligned, exp_silu_aligned)
if (ncol(exp_gse_aligned) > 0) exp_combined <- rbind(exp_combined, exp_gse_aligned)

# Build combined metadata (sample, group, batch)
meta_geo <- data.frame(
  sample_id = rownames(exp_geo_aligned),
  group = if (exists("ph") && "group" %in% colnames(ph)) ph$group[match(rownames(exp_geo_aligned), ph$geo_accession)] else NA_character_,
  batch = "GEO_GeneChip",
  stringsAsFactors = FALSE
)
meta_silu <- data.frame(
  sample_id = rownames(exp_silu_aligned),
  group = ifelse(!is.null(status_silu[rownames(exp_silu_aligned)]), status_silu[rownames(exp_silu_aligned)], NA_character_),
  batch = "silu_RNAseq",
  stringsAsFactors = FALSE
)
if (ncol(exp_gse_aligned) > 0) {
  meta_gse <- data.frame(
    sample_id = rownames(exp_gse_aligned),
    group = if (exists("exp_gse") && "status" %in% colnames(exp_gse)) exp_gse$status[match(rownames(exp_gse_aligned), rownames(exp_gse))] else NA_character_,
    batch = "GEO_RNAseq",
    stringsAsFactors = FALSE
  )
} else {
  meta_gse <- data.frame()
}

meta_combined <- rbind(meta_geo, meta_silu)
if (nrow(meta_gse) > 0) meta_combined <- rbind(meta_combined, meta_gse)
rownames(meta_combined) <- meta_combined$sample_id

# Ensure ordering
meta_combined <- meta_combined[rownames(exp_combined), , drop = FALSE]

# Save combined expression for record
saveRDS(exp_combined, OUT_COMBINED_RDS)
message("Saved combined expression to: ", OUT_COMBINED_RDS)

# ---------------------------
# Prepare for ComBat: ComBat expects genes x samples
# ---------------------------
exp_for_combat <- t(exp_combined)  # genes as rows, samples as columns

# Model matrix to retain biological group differences if available
if ("group" %in% colnames(meta_combined) && any(!is.na(meta_combined$group))) {
  mod <- model.matrix(~ group, data = meta_combined)
} else {
  mod <- NULL
}

batch_vector <- meta_combined$batch

# Run ComBat
message("Running ComBat for batch correction...")
exp_combat_corrected <- ComBat(dat = exp_for_combat, batch = batch_vector, mod = mod, par.prior = TRUE, prior.plots = FALSE)

exp_final <- as.data.frame(t(exp_combat_corrected))  # samples x genes
saveRDS(exp_final, OUT_COMBAT_RDS)
message("Saved ComBat-corrected expression to: ", OUT_COMBAT_RDS)

# ---------------------------
# Visualization before/after: PCA, t-SNE, UMAP
# ---------------------------
message("Generating pre/post batch-correction visualizations...")

# Choose a subset of top variable genes for dimensionality reduction
gene_vars_before <- apply(exp_for_combat, 1, var, na.rm = TRUE)
top_var_genes <- names(sort(gene_vars_before, decreasing = TRUE))[1:min(2000, length(gene_vars_before))]

exp_before_pca <- t(exp_for_combat[top_var_genes, , drop = FALSE])
exp_after_pca  <- t(exp_combat_corrected[top_var_genes, , drop = FALSE])

# PCA
pca_before <- prcomp(exp_before_pca, center = TRUE, scale. = TRUE)
pca_after <- prcomp(exp_after_pca, center = TRUE, scale. = TRUE)
var_explained_before <- round(100 * pca_before$sdev^2 / sum(pca_before$sdev^2), 2)
var_explained_after  <- round(100 * pca_after$sdev^2 / sum(pca_after$sdev^2), 2)

pca_before_df <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2], Batch = meta_combined$batch, Group = meta_combined$group, Sample = rownames(meta_combined))
pca_after_df  <- data.frame(PC1 = pca_after$x[,1],  PC2 = pca_after$x[,2],  Batch = meta_combined$batch, Group = meta_combined$group, Sample = rownames(meta_combined))

library(gridExtra)
p1 <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 2, alpha = 0.8) + labs(title = paste0("PCA before batch correction (PC1: ", var_explained_before[1], "%; PC2: ", var_explained_before[2], "%)")) + theme_bw()
p2 <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Batch, shape = Group)) +
  geom_point(size = 2, alpha = 0.8) + labs(title = paste0("PCA after batch correction (PC1: ", var_explained_after[1], "%; PC2: ", var_explained_after[2], "%)")) + theme_bw()

ggsave(file.path(plots_dir, "PCA_before_after.png"), arrangeGrob(p1, p2, ncol = 2), width = 14, height = 6, dpi = 300)
message("Saved PCA comparison to: ", file.path(plots_dir, "PCA_before_after.png"))

# ---------------------------
# Batch effect quantification: Silhouette score on first few PCs
# ---------------------------
message("Computing silhouette scores to quantify batch effect before/after correction...")
library(cluster)
n_pcs_for_sil <- min(10, ncol(pca_before$x))  # use up to 10 PCs
dist_before <- dist(pca_before$x[, 1:n_pcs_for_sil, drop = FALSE])
dist_after  <- dist(pca_after$x[, 1:n_pcs_for_sil, drop = FALSE])
sil_before <- silhouette(as.numeric(as.factor(meta_combined$batch)), dist_before)
sil_after  <- silhouette(as.numeric(as.factor(meta_combined$batch)), dist_after)
sil_score_before <- mean(sil_before[, 3])
sil_score_after  <- mean(sil_after[, 3])
message("Silhouette score (before): ", round(sil_score_before, 4))
message("Silhouette score (after):  ", round(sil_score_after, 4))

save(sil_before, sil_after, file = OUT_SIL_SCORE_RDATA)
message("Saved silhouette score objects to: ", OUT_SIL_SCORE_RDATA)

# End of script