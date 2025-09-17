
#########################################################################
# Drug Target Network Visualization (GDSC / CTRP datasets)
# Purpose:
# - Extract drug targets from GDSC/CTRP tables, classify targets into categories
# - Build a bipartite network (Target <-> Category) and visualize with ggraph/igraph
# - Save processed data and figures for manuscript use
#
# Usage:
# - Ensure 'drug1' and 'drug2' data.frames are loaded in the environment:
#     - For 'drug1', expected columns: gene_symbol_of_protein_target, target_or_activity_of_compound (optional)
#     - For 'drug2', expected column: Drug.target (or similar)
# - Set IO_DIR to the directory where you want to save outputs
#########################################################################

# Libraries
library(igraph)
library(ggraph)
library(tidyverse)

# ---------------------------
# User settings
# ---------------------------
IO_DIR <- getwd()  # change to your desired output directory
if (!dir.exists(IO_DIR)) dir.create(IO_DIR, recursive = TRUE)

OUT_DATA_RDATA <- file.path(IO_DIR, "drug_target_network_data.rds")
OUT_PLOT_PNG <- file.path(IO_DIR, "drug_target_network.png")

# ---------------------------
# Keyword to category mapping
# - ordered patterns: the first matching pattern assigns the category
# ---------------------------
keyword_class_map <- list(
  "VEGFR|KDR|FLT"   = "Tyrosine kinase",
  "PDGFR"           = "Tyrosine kinase",
  "KIT"             = "Tyrosine kinase",
  "MET"             = "Tyrosine kinase",
  "RET"             = "Tyrosine kinase",
  "FGFR"            = "Tyrosine kinase",
  "EGFR|ERBB"       = "Tyrosine kinase",
  "NTRK"            = "Tyrosine kinase",
  "BTK"             = "Tyrosine kinase",
  "TIE2|TEK"        = "Tyrosine kinase",
  "RON"             = "Tyrosine kinase",
  "IR|IGF1R"        = "Tyrosine kinase",
  "JAK1|JAK2|JAK3"  = "Tyrosine kinase",
  "PTK2"            = "Tyrosine kinase",
  
  "MAPK|ERK|JNK|p38" = "Ser/Thr kinase",
  "GSK3"             = "Ser/Thr kinase",
  "PLK"              = "Ser/Thr kinase",
  "PERK"             = "Ser/Thr kinase",
  "FAK"              = "Non-receptor tyrosine kinase (FAK)",
  
  "HDAC"             = "Epigenetic regulator",
  "EHMT2|GLP"        = "Epigenetic regulator",
  "USP|deubiquitinase" = "Epigenetic regulator",
  "E3-ubiquitin ligase" = "Epigenetic regulator",
  "EZH2"             = "Epigenetic regulator",
  "EHMT1"            = "Epigenetic regulator",
  
  "RXR"              = "Nuclear receptor",
  "VDR"              = "Nuclear receptor",
  
  "dopamine receptor|DRD" = "GPCR",
  "PAR1"             = "GPCR",
  
  "PIK3CD|PIK3CG"    = "Lipid Kinase",
  
  "CLTA|CLTB|CLTC|CLTCL1" = "Endocytosis-related protein",
  
  "DNA alkylator|DNA damage|alkylator" = "DNA damaging agent",
  
  "phosphodiesterase|PDE" = "Metabolic enzyme",
  "MCL1"                 = "Apoptosis regulator"
)

# ---------------------------
# Helper functions
# ---------------------------
# Normalize and split strings of targets; treat separators ; , / | as delimiters
split_targets_from_column <- function(vec) {
  vec <- vec[!is.na(vec) & trimws(vec) != ""]
  if (length(vec) == 0) return(character(0))
  tokens <- unlist(strsplit(paste(vec, collapse = ";"), "[;,/|]+"))
  tokens <- trimws(tokens)
  tokens <- tokens[tokens != ""]
  tokens <- unique(tokens)
  # Collapse multiple spaces
  tokens <- gsub("\\s+", " ", tokens)
  return(tokens)
}

# Classify a token using the mapping (first match wins)
classify_token <- function(token, mapping) {
  for (pat in names(mapping)) {
    if (grepl(pat, token, ignore.case = TRUE)) {
      return(mapping[[pat]])
    }
  }
  return("Other")
}

# ---------------------------
# Data presence checks
# ---------------------------
if (!exists("drug1") && !exists("drug2")) {
  stop("Please load at least one of 'drug1' or 'drug2' data.frames into the environment before running this script.")
}

# ---------------------------
# Extract targets from drug1 (if present)
# ---------------------------
drug1_targets <- character(0)
if (exists("drug1")) {
  # Attempt to find a column for gene symbol or target
  possible_cols1 <- intersect(c("gene_symbol_of_protein_target","target_or_activity_of_compound","target"), colnames(drug1))
  if (length(possible_cols1) == 0) {
    # try more flexible matching
    possible_cols1 <- colnames(drug1)[grepl("gene|target|protein", colnames(drug1), ignore.case = TRUE)]
  }
  if (length(possible_cols1) > 0) {
    col1 <- possible_cols1[1]
    drug1_targets <- split_targets_from_column(drug1[[col1]])
    message("Extracted ", length(drug1_targets), " unique tokens from drug1 column: ", col1)
  } else {
    message("No suitable target column found in drug1; skipping drug1 extraction.")
  }
}

# ---------------------------
# Extract targets from drug2 (if present)
# ---------------------------
drug2_targets <- character(0)
if (exists("drug2")) {
  possible_cols2 <- intersect(c("Drug.target","drug_target","target"), colnames(drug2))
  if (length(possible_cols2) == 0) {
    possible_cols2 <- colnames(drug2)[grepl("drug|target", colnames(drug2), ignore.case = TRUE)]
  }
  if (length(possible_cols2) > 0) {
    col2 <- possible_cols2[1]
    drug2_targets <- split_targets_from_column(drug2[[col2]])
    message("Extracted ", length(drug2_targets), " unique tokens from drug2 column: ", col2)
  } else {
    message("No suitable target column found in drug2; skipping drug2 extraction.")
  }
}

# ---------------------------
# Combine and clean tokens
# ---------------------------
all_targets <- unique(c(drug1_targets, drug2_targets))
if (length(all_targets) == 0) {
  stop("No target tokens extracted from drug1 or drug2. Please verify input columns.")
}
# Further split any compound entries using ; , / or |
split_tokens <- unique(unlist(strsplit(all_targets, "[,;/|]+")))
split_tokens <- trimws(split_tokens)
split_tokens <- split_tokens[split_tokens != ""]
split_tokens <- gsub("\\s+", " ", split_tokens)

# ---------------------------
# Classify tokens into categories
# ---------------------------
categories <- vapply(split_tokens, classify_token, FUN.VALUE = character(1), mapping = keyword_class_map)
target_df <- tibble(Target = split_tokens, Category = categories)

# ---------------------------
# Build edges and nodes for visualization
# ---------------------------
edges <- target_df %>% mutate(CategoryNode = paste0("Class:", Category)) %>% select(Target, CategoryNode)
nodes <- tibble(name = unique(c(edges$Target, edges$CategoryNode))) %>%
  mutate(type = ifelse(startsWith(name, "Class:"), "Class", "Target"))

# Add an index to edges for coloring (optional)
edges <- edges %>% mutate(Index = row_number())

# Create igraph object
graph <- graph_from_data_frame(edges, vertices = nodes)

# ---------------------------
# Plot network with ggraph
# ---------------------------
set.seed(123)
p <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(color = Index), alpha = 0.6, show.legend = FALSE) +
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Class" = "#D73027", "Target" = "#4575B4")) +
  scale_edge_colour_gradient(low = "lightblue", high = "darkblue") +
  ggtitle("Drug Target - Category Network (GDSC / CTRP)") +
  theme_void()

print(p)

# Save plot to PNG for manuscript (300 dpi)
ggsave(OUT_PLOT_PNG, plot = p, width = 10, height = 8, dpi = 300)
message("Saved network plot to: ", OUT_PLOT_PNG)

# ---------------------------
# Save processed data for reproducibility
# ---------------------------
saveRDS(list(target_df = target_df, edges = edges, nodes = nodes), file = OUT_DATA_RDATA)
message("Saved processed network data to: ", OUT_DATA_RDATA)

# End of script