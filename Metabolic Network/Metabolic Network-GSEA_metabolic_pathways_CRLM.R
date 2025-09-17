###################### Batch GSEA Using KEGG Database #######################
# GSEA analysis of differential expression
rm(list = ls())

# Load DEG results (Limma output)
nrDEG = read.csv(file = "nrDEG_grouping_results.csv")
DEG = nrDEG

# Prepare named numeric vector of log fold changes with ENTREZID names
geneList = as.numeric(DEG$logFC)
names(geneList) = DEG$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

# Remove duplicated SYMBOL entries if present
DEG$symbol = rownames(DEG)
DEG$gene = DEG$symbol
DEG = DEG[!duplicated(DEG$gene), ]
rownames(DEG) = DEG$gene

# Load required libraries for enrichment analysis
library(clusterProfiler)
library(hugene10sttranscriptcluster.db)
ids = toTable(hugene10sttranscriptclusterSYMBOL)

library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

# KEGG analysis: map SYMBOL to ENSEMBL/ENTREZ and prepare geneList
packageVersion('pheatmap')

# Trim ENSG IDs to remove version suffix when necessary
library(stringr)
rownames(DEG) <- str_sub(rownames(DEG), start = 1, end = 15)
DEG$symbol <- DEG$gene

# Map SYMBOL to ENSEMBL and ENTREZID
df <- bitr(DEG$symbol,
           fromType = "SYMBOL",
           toType = c("ENSEMBL", "ENTREZID"),
           OrgDb = org.Hs.eg.db)
# Note: some input gene IDs may fail to map

# Merge mapped ENTREZIDs into DEG table
colnames(DEG)[7] = "SYMBOL"
DEG <- inner_join(DEG, df, by = "SYMBOL")

# Prepare geneList for GSEA with ENTREZIDs
geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[geneList != 0]

class(geneList)

# Run GSEA using KEGG
kk_gse <- gseKEGG(
  geneList     = geneList,
  organism     = 'hsa',
  minGSSize    = 5,
  pvalueCutoff = 0.95,
  verbose      = FALSE
)
res = kk_gse@result

save(res, kk_gse, file = "KEGG_batch_GSEA.Rdata")

####################################
############################## Visualization ##################################

# 1. GSEA enrichment plots
library(clusterProfiler)
library(enrichplot)
library(cowplot) # or patchwork

# Select KEGG pathway IDs of interest
index <- c("hsa00830", "hsa00380")
selected_ids <- kk_gse@result$ID[kk_gse@result$ID %in% index]

# Generate GSEA plots for each selected pathway and assemble vertically
gsea_plots <- lapply(selected_ids, function(term_id) {
  gseaplot2(
    kk_gse,
    geneSetID = term_id,
    ES_geom = "line",
    pvalue_table = TRUE,
    title = kk_gse@result$Description[kk_gse@result$ID == term_id]
  )
})
plot_grid(plotlist = gsea_plots, ncol = 1, align = "v")

dev.off()

# Alternative single-call plotting for selected IDs
selected_gene_sets <- kk_gse@result$ID %in% index
gseaplot2(
  kk_gse,
  geneSetID = kk_gse@result$ID[selected_gene_sets],
  ES_geom = 'line',
  pvalue_table = TRUE
)

#############################
# GSEA plot for a single KEGG pathway (example)
g = gseaplot2(
  kk_gse,
  geneSetID = "hsa00380",
  color = "red",
  pvalue_table = TRUE,
  title = "Tryptophan Metabolism"
)
g

##############
