# ================ Load KEGG Pathway Gene Sets ===========================
# Read the TSV file for Tryptophan metabolism
data <- read.delim("KEGG_TRYPTOPHAN_METABOLISM.v2026.1.Hs.tsv", header = TRUE, sep = "\t")
class(data)
colnames(data) = data$STANDARD_NAME

# Inspect data structure and extract gene list
# Accessing the specific column and row containing the gene symbols
genes_tryptophan_raw = data$KEGG_TRYPTOPHAN_METABOLISM[17]
length(genes_tryptophan_raw)

# Split the comma-separated string into a vector of individual genes
genes_tryptophan <- strsplit(data$KEGG_TRYPTOPHAN_METABOLISM[17], ",")[[1]]
length(genes_tryptophan)

# -----------------------------------------------------------------
# Read the TSV file for Retinol metabolism
data2 <- read.delim("KEGG_RETINOL_METABOLISM.v2026.1.Hs.tsv", header = TRUE, sep = "\t")
class(data2)

# Inspect data structure and extract gene list
genes_retinol_raw = data2$KEGG_RETINOL_METABOLISM[17]
length(genes_retinol_raw)

# Split the comma-separated string into a vector of individual genes
genes_retinol <- strsplit(genes_retinol_raw, ",")[[1]]
length(genes_retinol)

# -----------------------------------------------------------------
# 4. Load the pre-defined list of 620 DEGs (Differentially Expressed Genes)
# Assuming the list is stored in the 'geneSymbols' column of the 'function_genes' object
deg_list_620 <- function_genes$geneSymbols

# 5. Perform intersection analysis to identify DEGs belonging to specific KEGG pathways
your_degs_in_tryptophan <- intersect(deg_list_620, genes_tryptophan)
your_degs_in_retinol <- intersect(deg_list_620, genes_retinol)

# 6. Display the results
cat("Number of DEGs identified in Tryptophan Metabolism: ", length(your_degs_in_tryptophan), "\n")
print(your_degs_in_tryptophan)

cat("\nNumber of DEGs identified in Retinol Metabolism: ", length(your_degs_in_retinol), "\n")
print(your_degs_in_retinol)

# 7. Save the resulting gene sets for downstream analysis
save(your_degs_in_retinol, your_degs_in_tryptophan, file = "tryptophan_and_retinol_genesets.rdata")
