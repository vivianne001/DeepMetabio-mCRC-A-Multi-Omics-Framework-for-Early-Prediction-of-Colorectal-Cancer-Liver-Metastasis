############# ssGSEA Analysis of the Tumor Immune Microenvironment ############
################################################################################

# Section 1: Setup and Environment Preparation

rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
# install.packages(c("tidyverse", "GSVA", "pheatmap", "ggpubr", "paletteer", "corrplot", "rstatix", "readxl")) # Uncomment if needed
library(tidyverse)
library(readxl)
library(GSVA)
library(pheatmap)
library(ggpubr)
library(paletteer)
library(corrplot)
library(rstatix)

# Section 2: Expression Data Preprocessing (Counts -> TPM)

# load("GSE204805_processed.Rdata")
gene_length_data <- read.csv("All_hg19gene_len.csv")

count_df <- as.data.frame(exp)
count_df$Gene <- rownames(count_df)

merged_data <- left_join(count_df, gene_length_data, by = "Gene")
merged_data <- na.omit(merged_data)
rownames(merged_data) <- merged_data$Gene

kb <- merged_data$Length / 1000
count_matrix <- merged_data[, colnames(exp)]
rpk <- count_matrix / kb
per_million_scaling_factor <- colSums(rpk) / 1e6
tpm_matrix <- t(t(rpk) / per_million_scaling_factor)
tpm_matrix <- as.data.frame(tpm_matrix)
tpm_matrix <- tpm_matrix[!duplicated(rownames(tpm_matrix)), ]
exp <- log2(tpm_matrix + 1)

# Section 3: Immune Gene Set Preparation

geneset_wide <- read_xlsx("GeneList.xlsx", sheet = 2)
geneset_long <- geneset_wide %>%
  pivot_longer(cols = everything(), names_to = "Cell.type", values_to = "Metagene") %>%
  na.omit() %>%
  arrange(Cell.type, Metagene)
immune_gene_sets <- split(geneset_long$Metagene, geneset_long$Cell.type)

# Section 4: ssGSEA Score Calculation and Normalization

ssgsea_scores <- gsva(
  expr = as.matrix(exp),
  gset.idx.list = immune_gene_sets,
  method = "ssgsea",
  kcdf = "Gaussian",
  mx.diff = FALSE,
  verbose = FALSE
)

# Min-max normalization across samples for each cell type score
normalize_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
ssgsea_scores_normalized <- t(apply(ssgsea_scores, 1, normalize_min_max))

# Save ssGSEA outputs
save(ssgsea_scores, ssgsea_scores_normalized, file = "ssGSEA_results.Rdata")

# Section 5: Visualization - Heatmap of Immune Scores

annotation_col <- data.frame(
  Group = ifelse(grepl("metastasis", colnames(exp)), "Metastasis", "Primary"),
  row.names = colnames(exp)
)

pheatmap(
  ssgsea_scores_normalized,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  clustering_method = "complete",
  fontsize = 8
)

# Section 6: Visualization - Boxplots of Immune Cell Scores by Group

plot_data <- as.data.frame(t(ssgsea_scores_normalized)) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Cell_Type", values_to = "Score")
plot_data$Group <- annotation_col$Group[match(plot_data$Sample, rownames(annotation_col))]

stat_test <- plot_data %>%
  group_by(Cell_Type) %>%
  wilcox_test(Score ~ Group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Cell_Type", dodge = 0.8)

p_boxplot <- ggboxplot(
  plot_data, x = "Cell_Type", y = "Score",
  fill = "Group", color = "black",
  width = 0.7, alpha = 0.7,
  outlier.shape = 21, outlier.size = 1.5
) +
  scale_fill_manual(values = c("#4979b6", "#d9352a")) +
  labs(x = "Immune Cell Type", y = "Normalized Enrichment Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  stat_pvalue_manual(stat_test, label = 'p.adj.signif', tip.length = 0.01, bracket.size = 0.5)

print(p_boxplot)

# Section 7: Visualization - Correlation Heatmap of Immune Cells

cor_matrix <- cor(t(ssgsea_scores_normalized))
cor_test_results <- cor.mtest(cor_matrix, conf.level = 0.95)
p_matrix <- cor_test_results$p

corrplot(
  cor_matrix,
  method = "color",
  order = "hclust",
  tl.cex = 0.8,
  tl.col = "black",
  p.mat = p_matrix,
  sig.level = c(0.001, 0.01, 0.05),
  insig = "label_sig",
  pch.cex = 0.8,
  pch.col = "white"
)

# Section 8: Correlation between Immune Cell Scores and Target Genes

target_genes <- c("ACMSD", "SLC2A14", "MMP9", "TIMP1", "CD44", "CXCR1")
target_genes <- intersect(target_genes, rownames(exp))

exp_target_genes <- exp[target_genes, , drop = FALSE]
combined_matrix <- rbind(ssgsea_scores_normalized, exp_target_genes)

full_cor_matrix <- cor(t(combined_matrix))
full_cor_test <- cor.mtest(full_cor_matrix, conf.level = 0.95)
full_p_matrix <- full_cor_test$p

cross_cor_values <- full_cor_matrix[rownames(ssgsea_scores_normalized), target_genes]
cross_p_values <- full_p_matrix[rownames(ssgsea_scores_normalized), target_genes]

sig_markers <- case_when(
  cross_p_values < 0.001 ~ "***",
  cross_p_values < 0.01  ~ "**",
  cross_p_values < 0.05  ~ "*",
  TRUE ~ ""
)
dim(sig_markers) <- dim(cross_cor_values)

pheatmap(
  cross_cor_values,
  display_numbers = sig_markers,
  number_color = "black",
  fontsize_number = 12,
  color = colorRampPalette(c("#053061", "white", "#67001F"))(100),
  border_color = "grey60",
  angle_col = 45,
  main = "Correlation between Immune Cells and Target Genes"
)