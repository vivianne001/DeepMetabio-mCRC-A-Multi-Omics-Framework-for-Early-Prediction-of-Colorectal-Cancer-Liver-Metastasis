rm(list = ls())

###################### Multi-panel heatmap for NF-kB pathway / core TRS features #########################

library("pheatmap")
library("RColorBrewer")
library(dplyr)

## Load dataset
# load RData files that contain expression matrix `exp` and sample annotation `ph` (and grouplist)
load("GSE50760_Raw_Data.Rdata","group_metastasis_primary.Rdata")

# Ensure sample names in metadata match expression matrix column names
identical(rownames(ph), colnames(exp))

# Prepare column annotation data frame for heatmap
ann_col = ph[, c("geo_accession", "group")]
ann_col = as.data.frame(ann_col)

# Convert the 'group' column to a factor and set the desired level order
ann_col$group = factor(ann_col$group, levels = c("metastasis", "primary"))

# Sort ann_col by group so samples are grouped in the heatmap
ann_col = ann_col[order(ann_col$group), ]
cols = rownames(ann_col)

# Remove the geo_accession column and keep group only for annotation
ann_col = ann_col[, -1]
ann_col = as.data.frame(ann_col)
rownames(ann_col) = cols
colnames(ann_col) = "Group"

ann_col$Group = as.factor(ann_col$Group)

# Define annotation colors
ann_color = list(Group = c(metastasis = "#F79647", primary = "#93CCDD"))

# Reorder expression matrix columns to match annotation ordering
k = rownames(ann_col)
exp = exp[, match(k, colnames(exp))]

## Import feature subset to be plotted
core_genes = read.table("Core TRS Features.txt", header = TRUE, col.names = TRUE)
core_genes = as.data.frame(core_genes)

# Subset expression matrix for the core feature genes
exp = exp[core_genes$TRUE., ]
dat = exp

# Transform to log2 scale (log2(x + 1))
dat = as.data.frame(dat)
dat = log2(dat + 1)
range(dat)

#############################
# Heatmap: row-scaling, no clustering (rows and columns ordered as provided)
p = pheatmap(dat,
             scale = "row",
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             cutree_rows = NA,
             cutree_cols = NA,
             treeheight_col = 30,
             treeheight_row = 30,
             border_color = "grey60",
             display_numbers = FALSE,
             fontsize_number = 6,
             number_format = "%.2f",
             number_color = "grey30",
             fontsize = 10, fontsize_row = 6, fontsize_col = 10,
             show_colnames = FALSE, show_rownames = TRUE,
             main = "Core TRS Feature Expression between Groups",
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             angle_col = "45",
             gaps_row = NULL,
             gaps_col = c(18), # insert a column gap to visually separate the two groups
             annotation_col = ann_col, annotation_row = NA,
             annotation = NA, annotation_colors = ann_color,
             annotation_legend = TRUE,
             annotation_names_col = TRUE, annotation_names_row = TRUE)

print(p)

dev.off()

################### Boxplots for individual features ######################

# Transpose data for ggplot-friendly format: samples as rows, genes as columns
dat = t(dat)
dat = as.data.frame(dat)

# Attach group information from annotation
dat$group = ann_col$Group
dat$group = as.factor(dat$group)

# Example statistical tests for specific genes:
# Perform t-test and Wilcoxon test between groups for chosen genes
t_test_result <- t.test(dat$NADSYN1[dat$group == "metastasis"],
                        dat$NADSYN1[dat$group == "primary"])

wilcox_result <- wilcox.test(dat$KRAS[dat$group == "metastasis"],
                             dat$KRAS[dat$group == "primary"])

# Load plotting and significance annotation packages
library(ggsignif)
library(ggplot2)

# Map p-value to significance symbol helper
map_pvalue_to_signif_mark <- function(pvalue) {
  if (pvalue < 0.001) {
    return("***")
  } else if (pvalue < 0.01) {
    return("**")
  } else if (pvalue < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Example: get significance symbols for the performed tests
significance_symbols_ttest <- map_pvalue_to_signif_mark(t_test_result$p.value)
significance_symbols_wilcox <- map_pvalue_to_signif_mark(wilcox_result$p.value)

# Compute outliers for an example variable (placeholder)
# outliers <- boxplot.stats(testPtype.pr$BRD_K30748066)$out

########### Single-gene boxplot example #############
# Ensure the Group factor has the expected order
dat$group <- factor(dat$group, levels = c("primary", "metastasis"))

# Plot boxplot for a single gene (NADK) with significance from t-test
ggplot(dat, aes(x = group, y = NADK, fill = group)) +
  geom_boxplot(width = 0.75, outlier.shape = NA) + # hide outliers for cleaner appearance
  scale_fill_manual(values = c("metastasis" = "red", "primary" = "lightblue")) +
  geom_signif(comparisons = list(c("metastasis", "primary")),
              test = "t.test",
              map_signif_level = TRUE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14)) +
  ylab("NADK expression (log2 + 1)") +
  xlab("Group: ACMSD metastasis vs primary")

#########################
############## Multi-gene boxplots (long format) ###################

options(stringsAsFactors = FALSE)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(reshape)
library(reshape2)
library(IOBR)

# Convert data to long format for multiple genes
TME_NEW = melt(dat)
colnames(TME_NEW) = c("Group", "Gene", "Express")

# Ensure Express is numeric
TME_NEW$Express <- as.numeric(TME_NEW$Express)
head(TME_NEW)

# Define plotting order: median expression in metastasis group (descending)
plot_order = TME_NEW[TME_NEW$Group == "metastasis", ] %>%
  group_by(Gene) %>%
  summarise(m = median(Express, na.rm = TRUE)) %>%
  arrange(desc(m)) %>%
  pull(Gene)

# Ensure Group levels and labels for plotting
TME_NEW$Group <- factor(TME_NEW$Group, levels = c("primary", "metastasis"))

# Create boxplots for each gene and add t-test significance labels
p <- ggplot(TME_NEW, aes(x = Gene, y = Express, fill = Group)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, outlier.colour = "black") +
  scale_fill_manual(values = c("primary" = "lightblue", "metastasis" = "red"),
                    labels = c("primary" = "primary CRC", "metastasis" = "CRLM")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  stat_compare_means(aes(group = Group), method = "t.test", label = "p.signif", label.y = max(TME_NEW$Express, na.rm = TRUE))

print(p)