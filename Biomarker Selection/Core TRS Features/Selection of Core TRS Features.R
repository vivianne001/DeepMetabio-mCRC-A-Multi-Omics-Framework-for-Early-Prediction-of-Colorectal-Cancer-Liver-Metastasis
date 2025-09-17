
rm(list = ls())
###################

# Load TRS correlation results (training/internal) and CRLM correlation table (external)
mCRC_genes = read.csv("trs_gene_correlations.csv")
mCRC_genes = as.data.frame(mCRC_genes)

CRLM_genes = read.csv("trs_gene_correlation in CRLM Samples.csv")
CRLM_genes = as.data.frame(CRLM_genes)

####################
# Select significant genes in mCRC dataset
sig_genes <- subset(mCRC_genes, as.numeric(P_value) < 0.05)

# Rank by absolute correlation and select top 50
top_genes <- head(sig_genes[order(-abs(as.numeric(sig_genes$Correlation))), ], 50)

# Save result
mCRC_top_genes <- top_genes
print(top_genes)

###################
# Select significant genes in CRLM dataset
sig_genes <- subset(CRLM_genes, as.numeric(P_Value) < 0.05)

# Rank by absolute correlation and select top 50
top_genes <- head(sig_genes[order(-abs(as.numeric(sig_genes$Correlation))), ], 50)

# Save result
CRLM_top_genes <- top_genes
print(CRLM_top_genes)

########
library(ggvenn)
library(ggplot2) # ggvenn is based on ggplot2

# Prepare gene lists for Venn diagram
# The list names will be used as legend labels
gene_lists <- list(
  `TOP 50 Risk Biomarkers in CRLM Samples` = CRLM_top_genes$Gene,
  `TOP 50 Risk Biomarkers in mCRC Samples` = mCRC_top_genes$Gene
)

# Draw Venn diagram with clear contrasting colors
ggvenn(
  gene_lists,
  fill_color = c("#6B98C4", "#F5867F"),
  stroke_size = 1.5,          # border line width
  set_name_size = 6,          # set name font size
  text_size = 17,             # inside-set count font size
  fill_alpha = 0.7            # fill transparency
) +
  labs(title = "Venn Diagram of TOP 50 Most Correlated Features with TRS") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

#############
# Write intersected core features to file
core_features <- intersect(mCRC_top_genes$Gene, CRLM_top_genes$Gene)
write.table(core_features, file = "Core TRS Features.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
