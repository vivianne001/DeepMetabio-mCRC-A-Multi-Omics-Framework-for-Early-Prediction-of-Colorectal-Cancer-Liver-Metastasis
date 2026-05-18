
##########################################################
# Correlation Analysis: Drug IC50 vs Gene Expression (GDSC)
# Purpose:
# - Compute correlation between a gene (e.g., ACMSD) and drug IC50 across CRC cell lines
# - Identify significant drug associations and visualize example correlations and group comparisons
# - Save results for downstream analyses
#
# Requirements:
# - R objects or files:
#   * expression matrix 'exp' (genes x samples) or RData containing it
#   * drug IC50 matrices (e.g., GDSC1_LNIC50, GDSC2_LNIC50) with same sample column names as 'exp'
#   * phenotype table 'ph' with sample metadata (rownames matching sample IDs, contains columns like 'GDSC Tissue descriptor 2')
# - Packages: tidyverse (or dplyr, ggplot2), ggpubr, corrplot, tinyarray (optional)
##########################################################

rm(list = ls())

# ---------------------------
# Libraries
# ---------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(corrplot)

# Optional: tinyarray for cor.one if available
# library(tinyarray)

# ---------------------------
# User-editable settings
# ---------------------------
IO_DIR <- getwd()  # change to your data folder as needed
OUT_RDATA <- file.path(IO_DIR, "dat_of_GDSC_CRC_cells_IC50.rdata")
OUT_COR_RESULTS <- file.path(IO_DIR, "GDSC_drug_gene_correlations.csv")
OUT_SCATTER_PNG <- file.path(IO_DIR, "PD173074_vs_ACMSD_scatter.png")
OUT_BOX_PNG <- file.path(IO_DIR, "Tamoxifen_by_ACMSD_group.png")

# Gene and datasets to analyze
gene_symbol <- "ACMSD"   # gene to correlate with drug IC50
# Provide names of drug matrices; the script expects objects like GDSC1_LNIC50, GDSC2_LNIC50
# If such objects exist in environment, they'll be used; otherwise load from RData before running
# Example: drug_matrix <- GDSC2_LNIC50
drug_matrix <- NULL  # set to your drug matrix object name (e.g., GDSC2_LNIC50)

# ---------------------------
# Data loading (user can load RData manually before running)
# ---------------------------
# If you have a processed RData file with 'exp', 'GDSC2_LNIC50', 'ph' etc., load it:
# load(file.path(IO_DIR, "GDSC_processed_data.RData"))

# Check required objects
if (!exists("exp")) stop("Expression matrix 'exp' not found. Load it before running this script.")
# choose drug matrix: prefer GDSC2_LNIC50 if present
if (is.null(drug_matrix)) {
  if (exists("GDSC2_LNIC50")) {
    drug <- get("GDSC2_LNIC50")
  } else if (exists("GDSC1_LNIC50")) {
    drug <- get("GDSC1_LNIC50")
  } else {
    stop("No drug IC50 matrix found (e.g., GDSC2_LNIC50 or GDSC1_LNIC50). Please load one into the workspace.")
  }
} else {
  if (exists(drug_matrix)) {
    drug <- get(drug_matrix)
  } else {
    stop("Specified drug_matrix object not found in environment: ", drug_matrix)
  }
}

if (!exists("ph")) {
  warning("Phenotype table 'ph' not found. Some downstream filtering (by tissue) will be skipped.")
  ph <- NULL
}

# ---------------------------
# Synchronize samples between expression and drug data
# ---------------------------
if (!identical(colnames(exp), colnames(drug))) {
  common_cols <- intersect(colnames(exp), colnames(drug))
  if (length(common_cols) == 0) stop("No overlapping sample IDs between 'exp' and drug matrix.")
  exp <- exp[, common_cols, drop = FALSE]
  drug <- drug[, common_cols, drop = FALSE]
}

# Convert gene expression and drug matrices to a combined data.frame
if (!(gene_symbol %in% rownames(exp))) stop("Gene symbol not found in expression matrix: ", gene_symbol)
dat_combined <- rbind(exp[gene_symbol, , drop = FALSE], drug)
dat_combined <- as.data.frame(t(dat_combined))  # rows = samples, cols = gene + drugs

# ---------------------------
# Filter to CRC cell lines if phenotype available
# ---------------------------
if (!is.null(ph)) {
  # Ensure sample identifiers match
  if (!all(rownames(dat_combined) %in% rownames(ph))) {
    # try to match using common sample naming
    common_ids <- intersect(rownames(dat_combined), rownames(ph))
    if (length(common_ids) == 0) {
      warning("No matching sample IDs between dat and 'ph'; skipping tissue filtering.")
    } else {
      dat_combined <- dat_combined[common_ids, , drop = FALSE]
      ph <- ph[common_ids, , drop = FALSE]
    }
  }
  if ("GDSC Tissue descriptor 2" %in% colnames(ph)) {
    intestine_cells <- ph[ph$`GDSC Tissue descriptor 2` == "large_intestine", , drop = FALSE]
    CRC_cells <- rownames(intestine_cells)
    if (length(CRC_cells) > 0) {
      dat_combined <- dat_combined[intersect(rownames(dat_combined), CRC_cells), , drop = FALSE]
    } else {
      warning("No CRC samples found in phenotype table based on 'GDSC Tissue descriptor 2'. Using all samples.")
    }
  } else {
    warning("'GDSC Tissue descriptor 2' column not found in 'ph'; using all samples.")
  }
}

# Remove rows with any NA
dat_combined <- na.omit(dat_combined)

# Save dat for reproducibility
save(dat_combined, file = OUT_RDATA)
message("Saved combined data to: ", OUT_RDATA)

# ---------------------------
# Correlation analysis: gene vs each drug
# ---------------------------
# Compute Pearson correlation between the gene column and each other column
gene_vector <- dat_combined[[gene_symbol]]
results <- data.frame(
  column = colnames(dat_combined),
  correlation = NA_real_,
  p.value = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(dat_combined))) {
  varname <- colnames(dat_combined)[i]
  test_res <- cor.test(gene_vector, dat_combined[[varname]], method = "pearson")
  results$correlation[i] <- as.numeric(test_res$estimate)
  results$p.value[i] <- test_res$p.value
}

# Remove the self-correlation row (gene vs itself) if present
results <- results[results$column != gene_symbol, ]

# Save correlation results
write.csv(results, OUT_COR_RESULTS, row.names = FALSE)
message("Saved correlation results to: ", OUT_COR_RESULTS)

# ---------------------------
# Identify significant hits (example threshold)
# ---------------------------
sig_threshold_p <- 0.05
sig_threshold_cor <- 0.3
k <- results %>% filter(p.value < sig_threshold_p & abs(correlation) > sig_threshold_cor)

# If k is empty, inform user
if (nrow(k) == 0) {
  message("No significant drug correlations found with thresholds p<", sig_threshold_p, " and |cor|>", sig_threshold_cor)
} else {
  message("Identified ", nrow(k), " significant correlations.")
}

# ---------------------------
# Example visualization: scatter plot for a selected drug (if present)
# ---------------------------
# Choose a drug column for plotting (change as needed)
example_drug <- "PD173074"   # change to a drug present in your data
if (example_drug %in% colnames(dat_combined)) {
  pdat <- dat_combined[, c(example_drug, gene_symbol), drop = FALSE]
  colnames(pdat) <- c("DrugIC50", "GeneExpr")
  # Pearson correlation
  cor_test <- cor.test(pdat$DrugIC50, pdat$GeneExpr, method = "pearson")
  cor_value <- round(cor_test$estimate, 4)
  p_value <- format(cor_test$p.value, digits = 4, scientific = TRUE)
  
  p_scatter <- ggplot(pdat, aes(x = DrugIC50, y = GeneExpr)) +
    geom_point(color = "#2c7fb8") +
    geom_smooth(method = "lm", color = "#d95f0e", se = TRUE) +
    annotate("text", x = Inf, y = Inf, label = paste0("Pearson: ", cor_value, "\nP-value: ", p_value), hjust = 1.1, vjust = 1.1, size = 4) +
    labs(x = paste0(example_drug, " (LN IC50)"), y = paste0(gene_symbol, " expression (log2)"), title = paste("Correlation:", gene_symbol, "vs", example_drug)) +
    theme_minimal(base_size = 12) +
    theme(axis.line = element_line(color = "black"))
  print(p_scatter)
  ggsave(OUT_SCATTER_PNG, p_scatter, width = 6, height = 4, dpi = 300)
  message("Saved scatter plot to: ", OUT_SCATTER_PNG)
} else {
  message("Example drug '", example_drug, "' not found in data; skipping scatter plot.")
}

# ---------------------------
# Example group comparison: violin/box plot for a named drug vs gene-derived groups
# ---------------------------
# Create a grouping variable based on median gene expression
group_col <- paste0(gene_symbol, "_group")
dat_box <- dat_combined
dat_box[[group_col]] <- ifelse(dat_box[[gene_symbol]] > median(dat_box[[gene_symbol]], na.rm = TRUE), "high", "low")
# Select a drug for group comparison
example_drug2 <- "Tamoxifen"  # change to a drug present in your data
if (example_drug2 %in% colnames(dat_box)) {
  pdat_box <- dat_box %>% select(all_of(c(group_col, example_drug2))) %>% filter(!is.na(.data[[example_drug2]]))
  colnames(pdat_box) <- c("group", "IC50")
  pdat_box$group <- factor(pdat_box$group, levels = c("high", "low"))
  # t-test
  t_test_result <- t.test(IC50 ~ group, data = pdat_box)
  p_box <- ggplot(pdat_box, aes(x = group, y = IC50, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("high" = "#d73027", "low" = "#4575b4")) +
    stat_compare_means(method = "t.test", label = "p.format") +
    theme_minimal(base_size = 12) +
    labs(x = paste0(gene_symbol, " group"), y = paste0(example_drug2, " (LN IC50)"), title = paste(example_drug2, "by", gene_symbol, "group"))
  print(p_box)
  ggsave(OUT_BOX_PNG, p_box, width = 6, height = 4, dpi = 300)
  message("Saved group comparison plot to: ", OUT_BOX_PNG)
} else {
  message("Example drug '", example_drug2, "' not found in data; skipping group comparison plot.")
}

# End of script