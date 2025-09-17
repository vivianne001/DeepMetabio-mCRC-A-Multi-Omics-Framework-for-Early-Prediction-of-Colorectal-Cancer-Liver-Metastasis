
###################### Identification and Visualization of Functional DEGs #########################

rm(list = ls())

#############################################################################
## Functional DEGs identification

# Load DEG result and metabolic gene list
load("grouplist-volcano-metastasis-DATA.P.Value 0.05.Rdata")
metabolic_gene = read.csv("metabolic genes.csv", stringsAsFactors = FALSE)

# Extract upregulated and downregulated DEGs
DEGS = DEG[DEG$change %in% c("up", "down"), ]

# (Optional) save metabolic gene list as backup
write.csv(metabolic_gene, file = "metabolic_genes_backup.csv", row.names = FALSE)

# Intersect DEGs and metabolic genes to obtain functional DEGs
dd = intersect(DEGS$gene, metabolic_gene$geneSymbols)
function_genes = metabolic_gene[metabolic_gene$geneSymbols %in% dd, ]

# Save functional DEGs as CSV and RData
# write.csv(function_genes, file = "functional_dataset.csv", row.names = FALSE)
save(function_genes, file = "functional_dataset.Rdata")

##############################################################################
## Venn plot of Functional DEGs vs Metabolic Genes

# Prepare gene name vectors for Venn plot
DEGS_NAME = DEGS$gene_name
meta_genes = metabolic_gene$geneSymbols

# Draw Venn plot using ggvenn
library(ggvenn)

venn_list <- list(
  "Metastasis_DEGs" = DEGS_NAME,
  "Metabolic_Genes" = meta_genes
)

p_venn <- ggvenn(
  venn_list,
  c("Metastasis_DEGs", "Metabolic_Genes"),
  show_elements = FALSE,
  show_percentage = TRUE,
  digits = 3,
  fill_color = c("#EF8A43", "#4865A9"),
  fill_alpha = 0.5,
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "#75A4C9",
  set_name_size = 5,
  text_color = "black",
  text_size = 5,
  label_sep = ","
)

# Save Venn plot to files (PDF and PNG)
ggsave("VENN_DEGs_metabolic_genes.pdf", plot = p_venn, width = 6, height = 6, dpi = 300)
ggsave("VENN_DEGs_metabolic_genes.png", plot = p_venn, width = 6, height = 6, dpi = 300)