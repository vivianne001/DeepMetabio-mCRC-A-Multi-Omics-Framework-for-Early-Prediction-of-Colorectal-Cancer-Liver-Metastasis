
# 1. Gene list
gene_list <- c("CYP27A1","FMO3","CYP8B1","CYP2C8","ACMSD","ADH1A","SLC6A1",
               "ALDH8A1","CA5A","CYP1A2","CYP2A6","BCHE","AKR1C4","ATP2B2",
               "DHTKD1","CYP4A11","PC","DHODH","HSD17B6","AKR1C1","CYP2C9","RDH16")

# 2. Convert gene symbols to ENTREZ IDs
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

gene_df <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_df$ENTREZID

# 3. GO enrichment analysis (all ontologies)
ego <- enrichGO(gene         = entrez_ids,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "ALL",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)

# 4. KEGG enrichment analysis
ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05,
                    pAdjustMethod= "BH")

# 5. GO dotplot (adjust showCategory as needed)
dotplot(ego, showCategory = 15, x = "p.adjust") + ggtitle("GO Enrichment (adjusted p-value)")

# 6. KEGG dotplot / additional GO plotting
# Note: the following line refers to 'gene_diff' which should be defined if needed.
# If you intend to re-run GO with a different gene set, replace 'gene_diff' accordingly.
# ego = enrichGO(gene = gene_diff, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)

# Draw dotplot for KEGG enrichment (or GO if preferred)
dotplot(ekegg, font.size = 12, showCategory = 10) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 55)) +
  scale_fill_gradient(low = "red", high = "blue") +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))

# Save figure (adjust filename and format as required by the journal)
ggsave(filename = "GO_Enrichment_of_functional_DEGs.pdf", width = 8, height = 10)