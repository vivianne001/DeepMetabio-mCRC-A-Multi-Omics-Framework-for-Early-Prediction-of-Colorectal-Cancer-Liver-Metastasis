
#############################################################
# Batch Correlation Heatmap (Publication-ready)
# Purpose:
# - Load precomputed correlation results and sample data
# - Compute Pearson correlations (ACMSD vs drugs) and p-values
# - Plot upper-triangle correlation heatmap with significance annotations
# - Save figure to PNG (and optionally PDF) for manuscript use
#############################################################

rm(list = ls())

# Required libraries
library(dplyr)
library(corrplot)

# ---------------------------
# User settings
# ---------------------------
IO_DIR <- getwd()  # change to your data directory if needed
# Expected files (will try to load these RData files / objects):
# - "result of cor of GDSC_cor 0.3.rdata" (should define 'k' and optionally 'plot_df')
# - "dat of GDSC-CRC cells IC50.rdata" (should define 'dat')
RDATA_COR <- file.path(IO_DIR, "result of cor of GDSC_cor 0.3.rdata")
RDATA_DAT <- file.path(IO_DIR, "dat of GDSC-CRC cells IC50.rdata")

OUT_PLOT_PNG <- file.path(IO_DIR, "correlation_heatmap_ACMSD_upper_triangle.png")
OUT_PLOT_PDF <- file.path(IO_DIR, "correlation_heatmap_ACMSD_upper_triangle.pdf")  # optional

# ---------------------------
# Load data
# ---------------------------
if (file.exists(RDATA_COR)) {
  load(RDATA_COR)  # expects 'k' and possibly 'plot_df' etc.
  message("Loaded correlation results from: ", RDATA_COR)
} else if (exists("k")) {
  message("Using existing object 'k' in environment.")
} else {
  stop("Correlation results not found. Provide RData containing 'k' or place file: ", RDATA_COR)
}

if (file.exists(RDATA_DAT)) {
  load(RDATA_DAT)  # expects 'dat' (samples x variables including ACMSD and drug columns)
  message("Loaded sample data from: ", RDATA_DAT)
} else if (exists("dat")) {
  message("Using existing object 'dat' in environment.")
} else {
  stop("Sample data not found. Provide RData containing 'dat' or place file: ", RDATA_DAT)
}

# ---------------------------
# Subset dat to drug columns in k (if k defines relevant drug columns)
# ---------------------------
if (!exists("k") || is.null(k$column)) stop("'k' does not appear to contain a 'column' vector of drug names.")
drug_names <- k$column
# Keep only columns present in dat
drug_names_present <- intersect(drug_names, colnames(dat))
if (length(drug_names_present) == 0) stop("None of the drug names in 'k$column' found in dat columns.")
dat <- dat[, drug_names_present, drop = FALSE]

# Optional: if plot_df exists and has cpd_name, match order
if (exists("plot_df") && "cpd_name" %in% colnames(plot_df)) {
  ordered_cols <- intersect(plot_df$cpd_name, colnames(dat))
  if (length(ordered_cols) > 0) {
    dat <- dat[, ordered_cols, drop = FALSE]
  }
}

# ---------------------------
# Ensure ACMSD present and make it first column
# ---------------------------
if (!("ACMSD" %in% colnames(dat))) stop("Column 'ACMSD' not found in 'dat'.")
dat2 <- dat[, c("ACMSD", setdiff(colnames(dat), "ACMSD")), drop = FALSE]

# ---------------------------
# Compute Pearson correlation matrix and pairwise p-values
# ---------------------------
resmcor <- cor(dat2, method = "pearson", use = "pairwise.complete.obs")

# Function to compute pairwise p-values matrix
cor_p_matrix <- function(mat) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA_real_, n, n)
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # try-catch to avoid failures
      tt <- tryCatch(cor.test(mat[, i], mat[, j], method = "pearson"), error = function(e) NULL)
      if (!is.null(tt)) {
        p.mat[i, j] <- p.mat[j, i] <- tt$p.value
      } else {
        p.mat[i, j] <- p.mat[j, i] <- NA_real_
      }
    }
  }
  diag(p.mat) <- NA_real_
  return(p.mat)
}

p.mat <- cor_p_matrix(dat2)

# ---------------------------
# Plot heatmap (upper triangle) with corrplot
# ---------------------------
max_abs <- max(abs(resmcor), na.rm = TRUE)
cols <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(200)

# Adjust margins to make room for left labels
old_par <- par(no.readonly = TRUE)
par(mar = c(2, 8, 4, 2))

# Create plot and save to PNG (and PDF optionally)
png(OUT_PLOT_PNG, width = 2000, height = 1600, res = 300)
corrplot(resmcor,
         method = "color",
         type = "upper",
         order = "original",
         tl.cex = 0.7,
         tl.col = "black",
         col = cols,
         cl.lim = c(-max_abs, max_abs),
         p.mat = p.mat,
         sig.level = c(.001, .01, .05),
         insig = "label_sig",
         pch.cex = 0.8,
         pch.col = "white",
         na.label = " ",
         na.label.col = "grey95",
         addgrid.col = "grey90")
dev.off()
message("Saved correlation heatmap (PNG) to: ", OUT_PLOT_PNG)

# Optional: save to PDF as vector
pdf(OUT_PLOT_PDF, width = 10, height = 8)
corrplot(resmcor,
         method = "color",
         type = "upper",
         order = "original",
         tl.cex = 0.8,
         tl.col = "black",
         col = cols,
         cl.lim = c(-max_abs, max_abs),
         p.mat = p.mat,
         sig.level = c(.001, .01, .05),
         insig = "label_sig",
         pch.cex = 0.8,
         pch.col = "white",
         na.label = " ",
         na.label.col = "grey95",
         addgrid.col = "grey90")
dev.off()
message("Saved correlation heatmap (PDF) to: ", OUT_PLOT_PDF)

# Restore previous par settings
par(old_par)

# End