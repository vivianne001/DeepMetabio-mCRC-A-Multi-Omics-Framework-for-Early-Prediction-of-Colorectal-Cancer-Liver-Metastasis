
rm(list = ls())
#############################################################################
## Import Core TRS Features KEGG enrichment results
KEGG_TRS = read.csv("ora_given_genes_kegg_results.csv")
KEGG_TRS = as.data.frame(KEGG_TRS)
KEGG_TRS = KEGG_TRS[KEGG_TRS$Adjusted.P.value < 0.05,]

#############################
## Import KEGG enrichment results from metabolomics (CRLM & ICRC)
library(readxl)

KEGG_metab = read_xlsx("dem_kegg_enrichment.xlsx", col_names = TRUE)
KEGG_metab = KEGG_metab[KEGG_metab$FDR < 0.05,]

####
# Compute intersection between transcriptional KEGG terms and metabolic KEGG terms
kegg_intersect <- intersect(KEGG_TRS$Term, KEGG_metab$name)

# --------------------------------------------------------------------------
# Draw Venn plot
library(ggvenn)

venn_list <- list(
  "Metabolites_Enrichment" = KEGG_metab$name,
  "Core_TRS_Enrichment" = KEGG_TRS$Term
)

# Basic Venn plot (default orientation)
p1 <- ggvenn(
  venn_list,
  c("Metabolites_Enrichment", "Core_TRS_Enrichment"),
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
print(p1)

# Flipped Venn plot (top/bottom circles inverted)
p2 <- ggvenn(
  venn_list,
  c("Metabolites_Enrichment", "Core_TRS_Enrichment"),
  show_elements = FALSE,
  show_percentage = TRUE,
  digits = 6,
  fill_color = c("#CC7979", "#BAD1E4"),
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
) + coord_flip()
print(p2)

# Save the flipped Venn plot to a PDF file in the working directory
# (Use ggsave to ensure the ggplot object is stored, rather than saving an example plot)
ggsave("VENN_PLOT_DEGs_and_metabolic_genes_flipped.pdf", plot = p2, width = 8, height = 6)