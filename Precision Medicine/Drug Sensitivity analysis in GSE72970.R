
###################### 
# Visualization of Drug Response Proportions by Treatment Group in GSE72970
# Purpose:
# - Visualize treatment response proportions by ACMSD expression groups (high vs low)
# - Separate plots for targeted therapy group and chemotherapy group
# - Save figures for manuscript use (modify output paths as needed)
#
# Instructions:
# - Place this script in your working directory or set IO_DIR to the folder that contains the processed data objects.
# - This script expects two data objects in the environment:
#     - ph: a data.frame of phenotype/sample annotation
#     - exp: a gene expression matrix or data.frame with genes in rows and sample IDs in columns
#   If you have an RData file that loads these objects, load it before running this script.
######################

# Clear workspace (optional)
rm(list = ls())

# Required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)

# ---------------------------
# User-editable settings
# ---------------------------
# Set IO_DIR to where your data / outputs should be read/written
# Example: IO_DIR <- "D:/your/path/to/data"
IO_DIR <- getwd()  # default to current working directory; change as needed

# Optionally set output figure filenames
OUT_PLOT_TARGET <- file.path(IO_DIR, "Targeted_Therapy_Response.png")
OUT_PLOT_CHEMO  <- file.path(IO_DIR, "Chemotherapy_Response.png")

# ---------------------------
# Pre-checks and data loading
# ---------------------------
# The script assumes 'ph' (phenotype) and 'exp' (expression) are already available in the environment.
# If you saved a processed workspace, you can load it here:
# load(file.path(IO_DIR, "GSE72970_processed.Rdata"))

if (!exists("ph") || !exists("exp")) {
  stop("Required objects 'ph' and 'exp' are not found in the environment. Please load them before running this script.")
}

# Inspect key columns used by this script
# For synchronous metastasis (filtering CRLM)
if (!("synchronous metastase:ch1" %in% colnames(ph))) {
  stop("Column 'synchronous metastase:ch1' not found in phenotype table 'ph'. Please check column names.")
}

# Subset to synchronous metastasis == 'Yes' (CRLM samples)
ph <- ph[ph$`synchronous metastase:ch1` == "Yes", ]

# Synchronize sample order between expression and phenotype
common_samples <- intersect(colnames(exp), rownames(ph))
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between 'exp' columns and 'ph' rownames. Please verify sample identifiers.")
}
dat <- exp[, common_samples, drop = FALSE]
ph  <- ph[common_samples, , drop = FALSE]
stopifnot(all(colnames(dat) == rownames(ph)))

# ---------------------------
# Group assignment by treatment regimen and response
# ---------------------------
# Create treatment group: 'target' vs 'chemotherapy'
ph$group <- ifelse(
  ph$characteristics_ch1.8 %in% c(
    "regimen: FOLFIRI+BEVACIZUMAB",
    "regimen: FOLFIRI+ERBITUX",
    "regimen: FOLFIRINOX+BEVACIZUMAB",
    "regimen: FOLFOX+BEVACIZUMAB"
  ),
  "target", "chemotherapy"
)

# Ensure response category column exists
if (!("response category:ch1" %in% colnames(ph))) {
  stop("Column 'response category:ch1' not found in phenotype table 'ph'.")
}

# Prepare data.frame for plotting: samples as rows
dat_df <- as.data.frame(t(dat))
dat_df$group <- ph$group
dat_df$outcome <- ph$`response category:ch1`

# ---------------------------
# Stratify by ACMSD expression (median split)
# ---------------------------
if (!("ACMSD" %in% colnames(dat_df))) {
  stop("Gene 'ACMSD' not found in expression matrix. Ensure gene names match and are present in 'exp'.")
}
median_score <- median(dat_df$ACMSD, na.rm = TRUE)
dat_df$ExprGroup <- ifelse(dat_df$ACMSD >= median_score, "high", "low")

# ---------------------------
# Separate by treatment group
# ---------------------------
dat_target <- dat_df[dat_df$group == "target", ]
dat_chemo  <- dat_df[dat_df$group == "chemotherapy", ]

# ---------------------------
# Function to prepare summary and plot stacked proportion bar
# ---------------------------
plot_response_proportions <- function(df, title, out_file = NULL, palette = NULL, response_levels = NULL) {
  if (nrow(df) == 0) {
    message("No samples for plotting group: ", title)
    return(NULL)
  }
  
  # Prepare summary counts and proportions
  summary_df <- df %>%
    group_by(ExprGroup, outcome) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(ExprGroup) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  # Ensure outcome is a factor with specified order if provided
  if (!is.null(response_levels)) {
    summary_df$outcome <- factor(summary_df$outcome, levels = response_levels)
  }
  
  # Prepare x-axis labels to include group sizes
  group_counts <- table(df$ExprGroup)
  x_labels <- paste(names(group_counts), "(n=", as.integer(group_counts), ")", sep = "")
  
  # Default palette
  if (is.null(palette)) {
    palette <- c("PD" = "#BA3E45", "SD" = "#20AEDD", "PR" = "#B6D7E9", "CR" = "#E1F3FB")
  }
  
  p <- ggplot(summary_df, aes(x = ExprGroup, y = prop, fill = outcome)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = palette, na.value = "grey80") +
    labs(x = "Expression Group", y = "Proportion", fill = "Treatment Response", title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.line = element_line(color = "grey50"),
      panel.grid.major.x = element_blank(),
      legend.position = "right"
    ) +
    scale_x_discrete(labels = x_labels)
  
  print(p)
  
  # Save plot if requested
  if (!is.null(out_file)) {
    ggsave(filename = out_file, plot = p, width = 6, height = 4.5, dpi = 300)
    message("Saved plot to: ", out_file)
  }
  
  return(p)
}

# ---------------------------
# Plot Targeted therapy response
# ---------------------------
target_levels <- c("PD", "SD", "PR")
p_target <- plot_response_proportions(dat_target, title = "Targeted Therapy Response", out_file = OUT_PLOT_TARGET, response_levels = target_levels)

# ---------------------------
# Plot Chemotherapy response
# ---------------------------
chemo_levels <- c("PD", "SD", "PR", "CR")
p_chemo <- plot_response_proportions(dat_chemo, title = "Chemotherapy Response", out_file = OUT_PLOT_CHEMO, response_levels = chemo_levels)

# End of script