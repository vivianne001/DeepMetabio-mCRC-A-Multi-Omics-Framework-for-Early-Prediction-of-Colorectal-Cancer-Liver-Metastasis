
############################################################
# Drug Target Attribute Classification and Visualization
# Purpose:
# - Map drugs to target categories using a provided 'edges' table
# - Merge drug-category mapping with drug correlation table 'k'
# - Produce a category-grouped correlation heatmap and save outputs
#
# Expected inputs (either as CSV files in working dir, or as R objects in the environment):
# - k (or k.csv)       : columns include 'column' (drug name), 'correlation', 'p.value'
# - drug2 (or drug2.csv): columns include 'Drug.name', 'Drug.target'
# - edges (or edges.csv): columns include 'Target', 'CategoryNode' (e.g., "Class:Tyrosine kinase")
#
# Set IO_DIR below to control output locations
############################################################

rm(list = ls())

# Required libraries
library(ggplot2)
library(dplyr)
library(scales)
library(readxl)
library(corrplot)

# ---------------------------
# User-editable settings
# ---------------------------
IO_DIR <- getwd()              # change to where your data and outputs should live
OUT_CSV <- file.path(IO_DIR, "drug_correlation_by_category.csv")
OUT_PLOT_PNG <- file.path(IO_DIR, "drug_correlation_by_category.png")
OUT_RDATA <- file.path(IO_DIR, "plot_df_of_drug_and_categories.rdata")

# ---------------------------
# Load data (CSV fallback or existing objects)
# ---------------------------
if (file.exists(file.path(IO_DIR, "k.csv"))) {
  k <- read.csv(file.path(IO_DIR, "k.csv"), stringsAsFactors = FALSE)
} else if (!exists("k")) {
  stop("k data not found. Provide 'k.csv' in IO_DIR or an R object named 'k' in the environment.")
}

if (file.exists(file.path(IO_DIR, "drug2.csv"))) {
  drug2 <- read.csv(file.path(IO_DIR, "drug2.csv"), stringsAsFactors = FALSE)
} else if (!exists("drug2")) {
  stop("drug2 data not found. Provide 'drug2.csv' in IO_DIR or an R object named 'drug2'.")
}

if (file.exists(file.path(IO_DIR, "edges.csv"))) {
  edges <- read.csv(file.path(IO_DIR, "edges.csv"), stringsAsFactors = FALSE)
} else if (!exists("edges")) {
  stop("edges data not found. Provide 'edges.csv' in IO_DIR or an R object named 'edges'.")
}

# ---------------------------
# Basic sanity checks
# ---------------------------
required_k_cols <- c("column", "correlation", "p.value")
if (!all(required_k_cols %in% colnames(k))) stop("k must contain columns: column, correlation, p.value")

if (!all(c("Drug.name", "Drug.target") %in% colnames(drug2))) stop("drug2 must contain columns: Drug.name, Drug.target")
if (!all(c("Target", "CategoryNode") %in% colnames(edges))) stop("edges must contain columns: Target, CategoryNode")

# ---------------------------
# Helper: split a target string into tokens
# ---------------------------
split_targets_string <- function(s) {
  if (is.na(s) || s == "") return(character(0))
  toks <- unlist(strsplit(as.character(s), "[,;/|]+"))
  toks <- trimws(toks)
  toks <- toks[toks != ""]
  return(toks)
}

# ---------------------------
# Map drug target tokens to category nodes using edges table
# ---------------------------
drug2_map <- drug2 %>%
  rowwise() %>%
  mutate(tokens = list(split_targets_string(Drug.target))) %>%
  ungroup()

map_tokens_to_categories <- function(tokens, edges_df) {
  if (length(tokens) == 0) return(NA_character_)
  cats <- character(0)
  for (tok in tokens) {
    # match any edges$Target pattern inside token (case-insensitive)
    matched <- edges_df$CategoryNode[
      vapply(edges_df$Target, function(t) {
        if (is.na(t) || t == "") return(FALSE)
        grepl(t, tok, ignore.case = TRUE)
      }, logical(1))
    ]
    if (length(matched) > 0) cats <- c(cats, matched)
  }
  cats <- unique(cats)
  if (length(cats) == 0) return(NA_character_)
  return(cats)
}

drug2_map$Categories <- lapply(drug2_map$tokens, map_tokens_to_categories, edges_df = edges)

# Choose a primary category per drug (join multiple matches with ';' if needed)
drug2_map <- drug2_map %>%
  mutate(PrimaryCategory = sapply(Categories, function(x) {
    if (all(is.na(x))) return(NA_character_)
    if (length(x) == 1) return(x)
    paste(x, collapse = ";")
  }))

drug_category <- drug2_map %>% transmute(Drug.name, PrimaryCategory = ifelse(is.na(PrimaryCategory), "Class:Unknown", PrimaryCategory))

# ---------------------------
# Merge correlation table 'k' with drug categories
# ---------------------------
plot_df <- k %>% rename(Drug.name = column, correlation = correlation, p.value = p.value) %>%
  left_join(drug_category, by = "Drug.name")

# For missing categories, try partial match on drug name
missing_idx <- which(is.na(plot_df$PrimaryCategory))
if (length(missing_idx) > 0) {
  for (i in missing_idx) {
    name <- plot_df$Drug.name[i]
    mm <- drug2$Drug.name[grepl(name, drug2$Drug.name, ignore.case = TRUE)]
    if (length(mm) >= 1) {
      plot_df$PrimaryCategory[i] <- drug_category$PrimaryCategory[match(mm[1], drug_category$Drug.name)]
    } else {
      plot_df$PrimaryCategory[i] <- "Class:Unknown"
    }
  }
}

# Create significance labels for display
plot_df <- plot_df %>%
  mutate(signif_label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

# Order drugs by PrimaryCategory then by absolute correlation descending
plot_df <- plot_df %>%
  group_by(PrimaryCategory) %>%
  arrange(PrimaryCategory, -abs(correlation)) %>%
  ungroup() %>%
  mutate(Drug = factor(Drug.name, levels = unique(Drug.name)))

# Save processed table for inspection
write.csv(plot_df, file = OUT_CSV, row.names = FALSE)
message("Saved merged table to: ", OUT_CSV)

# ---------------------------
# Plot: tile heatmap grouped by category
# ---------------------------
max_abs_cor <- max(abs(plot_df$correlation), na.rm = TRUE)
cor_limits <- c(-max_abs_cor, max_abs_cor)

p <- ggplot(plot_df, aes(x = PrimaryCategory, y = Drug, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif_label), color = "black", hjust = 1.05, size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = cor_limits, oob = scales::squish, name = "Correlation") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(title = "Drug correlation with ACMSD (grouped by category)")

print(p)

# Save figure as PNG for manuscript use
ggsave(filename = OUT_PLOT_PNG, plot = p, width = 10, height = 8, dpi = 300)
message("Saved plot to: ", OUT_PLOT_PNG)

# Optionally save R data for reproducibility
save(plot_df, file = OUT_RDATA)
message("Saved plot_df R data to: ", OUT_RDATA)

# ---------------------------
# Optional correlation matrix visualization (ACMSD vs drugs) using corrplot
# - Build an (n+1)x(n+1) matrix where first row/col is ACMSD and others are drugs
# ---------------------------
# Prepare plot_df ordering and drug list
plot_df$Drug <- as.character(ifelse(is.null(plot_df$Drug), plot_df$Drug.name, plot_df$Drug))
plot_df <- plot_df %>% arrange(PrimaryCategory, -abs(correlation))
drugs <- plot_df$Drug
n <- length(drugs)

if (n > 0) {
  mat <- matrix(NA_real_, nrow = n + 1, ncol = n + 1)
  rownames(mat) <- colnames(mat) <- c("ACMSD", drugs)
  diag(mat) <- 1
  mat["ACMSD", drugs] <- plot_df$correlation
  mat[drugs, "ACMSD"] <- plot_df$correlation
  
  p.mat <- matrix(NA_real_, nrow = n + 1, ncol = n + 1)
  rownames(p.mat) <- colnames(p.mat) <- c("ACMSD", drugs)
  diag(p.mat) <- 0
  p.mat["ACMSD", drugs] <- plot_df$p.value
  p.mat[drugs, "ACMSD"] <- plot_df$p.value
  
  # Color palette symmetric around zero
  col_palette <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(200)
  
  corrplot(mat,
           method = "color",
           type = "full",
           order = "original",
           tl.cex = 0.7,
           tl.col = "black",
           col = col_palette,
           cl.lim = cor_limits,
           p.mat = p.mat,
           sig.level = c(.001, .01, .05),
           insig = "label_sig",
           pch.cex = 0.8,
           pch.col = "white",
           na.label = " ",
           na.label.col = "grey95")
}

# End of script