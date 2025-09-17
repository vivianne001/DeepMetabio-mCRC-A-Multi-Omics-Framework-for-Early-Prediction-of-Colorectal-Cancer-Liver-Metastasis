
######################## Bulk GSEA Analysis Using KEGG Database ########################
# GSE204805
# Clear workspace
rm(list = ls())

# Load DEG results
load("ACMSD group_high and low.rdata")
DEG = nrDEG

# Prepare ranked geneList for GSEA (named vector: logFC with ENTREZID names)
geneList = as.numeric(DEG$logFC)
names(geneList) = DEG$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

# Remove duplicated gene symbols if present
DEG$symbol = rownames(DEG)
DEG$gene = DEG$symbol
DEG = DEG[!duplicated(DEG$gene), ]
rownames(DEG) = DEG$gene

# Load required libraries for GSEA and annotation
library(clusterProfiler)
library(hugene10sttranscriptcluster.db)
ids = toTable(hugene10sttranscriptclusterSYMBOL)

library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

# Check package version (optional)
packageVersion('pheatmap')

# Trim ENSG IDs to 15 characters if they contain version suffixes
library(stringr)
rownames(DEG) <- str_sub(rownames(DEG), start = 1, end = 15)
DEG$symbol <- DEG$gene

# Map gene SYMBOLs to ENSEMBL and ENTREZID IDs
df <- bitr(DEG$symbol,
           fromType = "SYMBOL",
           toType = c("ENSEMBL", "ENTREZID"),
           OrgDb = org.Hs.eg.db)

# Merge ENTREZID into DEG dataframe
colnames(DEG)[7] = "SYMBOL"
DEG <- inner_join(DEG, df, by = "SYMBOL")

# Prepare geneList again using mapped ENTREZIDs and filter out zeros
geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[geneList != 0]

# Run GSEA using KEGG (gseKEGG)
kk_gse <- gseKEGG(
  geneList     = geneList,
  organism     = 'hsa',
  minGSSize    = 5,
  pvalueCutoff = 0.95,
  verbose      = FALSE
)
res = kk_gse@result

save(res, kk_gse, file = "KEGG_Bulk_GSEA.Rdata")

############################## Visualization ##############################

# 1. GSEA enrichment plots for selected pathways
library(enrichplot)
library(ggplot2)

# Specify target KEGG pathway IDs to visualize
index = c("hsa04350", "hsa04064")
selected_gene_sets <- kk_gse@result$ID %in% index

# Plot enrichment curves for selected gene sets
gseaplot2(
  kk_gse,
  geneSetID = kk_gse@result$ID[selected_gene_sets],
  ES_geom = 'line',
  pvalue_table = TRUE
)

# Example: GSEA plot for a single KEGG pathway
g = gseaplot2(
  kk_gse,
  geneSetID = "hsa04512",
  color = "red",
  pvalue_table = TRUE,
  title = "ECM-receptor interaction"
)
g

